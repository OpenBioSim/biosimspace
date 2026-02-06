######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2025
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BioSimSpace is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""
Base class for wrapped Sire objects. This is an internal package and should
not be directly exposed to the user.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["SireWrapper"]


class SireWrapper:
    """A base class for wrapping Sire objects."""

    def __init__(self, object):
        """
        Constructor.

        Parameters
        ----------

        object : Sire.System.System, Sire.Mol.Molecule, Sire.Mol.Residue, Sire.Mol.Atom
            A Sire object.
        """

        # Store a deep copy of the Sire object.
        self._sire_object = object.__deepcopy__()

        # Initialise flags.
        self._is_multi_atom = False
        self._is_perturbable = False

    def __eq__(self, other):
        """Equals to operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._sire_object == other._sire_object
        else:
            return False

    def __ne__(self, other):
        """Not equals to operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._sire_object != other._sire_object
        else:
            return False

    def __hash__(self):
        """Hash operator."""
        return hash(self._sire_object)

    def copy(self):
        """
        Return a copy of this object. The return type is same as the object
        on which copy is called.

        Returns
        -------

        object : :class:`Atom <BioSimSpace._SireWrappers.Atom>`, \
                 :class:`Residue <BioSimSpace._SireWrappers.Residue>`, \
                 :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                 :class:`Molecules <BioSimSpace._SireWrappers.Molecules>`, \
                 :class:`System <BioSimSpace._SireWrappers.System>`
            A copy of the object.
        """
        return type(self)(self)

    def charge(self, property_map={}, is_lambda1=False):
        """
        Return the charge.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        is_lambda1 : bool
           Whether to use the charge at lambda = 1 if the object is peturbable.

        Returns
        -------

        charge : :class:`Charge <BioSimSpace.Types.Charge>`
            The charge.
        """
        from .. import Units as _Units

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'.")

        if not isinstance(is_lambda1, bool):
            raise TypeError("'is_lambda1' must be of type 'bool'.")

        # Copy the map.
        _property_map = property_map.copy()

        if property_map == {}:
            # This is a perturbable molecule.
            if self._is_perturbable:
                # Compute the charge for the chosen end state.
                if is_lambda1:
                    _property_map = {"charge": "charge1"}
                else:
                    _property_map = {"charge": "charge0"}

        # Calculate the charge.
        try:
            charge = self._sire_object.evaluate().charge(_property_map).value()
        except:
            charge = 0

        # Return the charge.
        return charge * _Units.Charge.electron_charge

    def smiles(self):
        """
        Return the SMILES string representation of this object.

        Returns
        -------

        smiles : str, [str]
            The SMILES string representation of the object(s).
        """

        try:
            return self._sire_object.smiles()
        except:
            return [obj.smiles() for obj in self]

    def translate(self, vector, property_map={}):
        """
        Translate the object.

        Parameters
        ----------

        vector : [:class:`Length <BioSimSpace.Types.Length>`]
            The translation vector in Angstroms.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """
        from .._Exceptions import IncompatibleError as _IncompatibleError
        from .. import _isVerbose
        from ..Types import Length as _Length
        from sire.legacy import Maths as _SireMaths

        # Convert tuple to a list.
        if isinstance(vector, tuple):
            vector = list(vector)

        # Validate input.
        if isinstance(vector, list):
            vec = []
            for x in vector:
                if type(x) is int:
                    vec.append(float(x))
                elif isinstance(x, float):
                    vec.append(x)
                elif isinstance(x, _Length):
                    vec.append(x.angstroms().value())
                else:
                    raise TypeError(
                        "'vector' must contain 'int', 'float', or "
                        "'BioSimSpace.Types.Length' types only!"
                    )
        else:
            raise TypeError("'vector' must be of type 'list' or 'tuple'")

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        try:
            # Make a local copy of the property map.
            _property_map = property_map.copy()

            # This is a perturbable molecule. First translate the lamba = 0 state.
            if self._is_perturbable:
                _property_map["coordinates"] = "coordinates0"

            # Perform the translation.
            self._sire_object = (
                self._sire_object.move()
                .translate(_SireMaths.Vector(vec), _property_map)
                .commit()
            )

            # This is a perturbable molecule. Now translate the lamba = 1 state.
            if self._is_perturbable:
                _property_map["coordinates"] = "coordinates1"

                # Perform the translation.
                self._sire_object = (
                    self._sire_object.move()
                    .translate(_SireMaths.Vector(vec), _property_map)
                    .commit()
                )

        except UserWarning as e:
            msg = (
                "Cannot compute axis-aligned bounding box "
                + "since the object has no 'coordinates' property."
            )
            if _isVerbose():
                raise _IncompatibleError(msg) from e
            else:
                raise _IncompatibleError(msg) from None

    def getAxisAlignedBoundingBox(self, property_map={}):
        """
        Get the axis-aligned bounding box enclosing the object.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        box_min : [:class:`Length <BioSimSpace.Types.Length>`]
            The minimum coordinates of the axis-aligned bounding box in
            each dimension.

        box_max : [:class:`Length <BioSimSpace.Types.Length>`]
            The minimum coordinates of the axis-aligned bounding box in
            each dimension.
        """
        from .. import Units as _Units

        aabox = self._getAABox(property_map)

        box_min = [x.value() * _Units.Length.angstrom for x in aabox.min_coords()]
        box_max = [x.value() * _Units.Length.angstrom for x in aabox.max_coords()]

        return box_min, box_max

    def _getCenterOfMass(self, space=None, property_map={}):
        """
        Get the center of mass for this object.

        Parameters
        ----------

        space : sire.legacy.Vol.PeriodicBox, sire.legacy.Vol.TriclinicBox
            The space associated with the object.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        com: [:class:`Length <BioSimSpace.Types.Length>`]
            The center of mass of the object.
        """
        from .. import Units as _Units
        from sire.legacy import Vol as _SireVol

        if space is None:
            space_prop = property_map.get("space", "space")
            try:
                space = self._sire_object.property("space")
            except:
                space = _SireVol.Cartesian()
        else:
            if not isinstance(
                space, (_SireVol.PeriodicBox, _SireVol.TriclinicBox, _SireVol.Cartesian)
            ):
                raise TypeError(
                    "'space' must be of type 'sire.legacy.Vol.PeriodicBox', "
                    "'sire.legacy.Vol.TriclinicBox', or 'sire.legacy.Vol.Cartesian'"
                )

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'.")

        # Whether this is an atom, molecule or residue, or system.
        is_atom = False
        is_mol_res = False
        is_system = False

        # Get the first atom in the object.
        try:
            atom = self[0].getAtoms()[0]
            is_system = True
        except:
            try:
                atom = self.getAtoms()[0]
                is_mol_res = True
            except:
                atom = self
                is_atom = True

        # Get the name of the required properties.
        mass_prop = property_map.get("mass", "mass")
        coord_prop = property_map.get("coordinates", "coordinates")

        # Intialise the total mass.
        try:
            total_mass = atom._sire_object.property(mass_prop).value()
        except:
            ValueError("Unable to compute center of mass. Missing 'mass' property!")

        # Store the reference coordinate.
        try:
            ref_coord = atom._sire_object.property(coord_prop)
        except:
            ValueError(
                "Unable to compute center of mass. Missing 'coordinates' property!"
            )

        # Intialise the center of mass.
        com = total_mass * ref_coord

        # If this is a single atom, then return immediately.
        if is_atom:
            return com

        def update_com(atoms, com, total_mass):
            """
            Helper function to compute the center of mass of a set of atoms.
            """
            from sire.legacy import Maths as _SireMaths

            for atom in atoms:
                # Update the total mass.
                try:
                    mass = atom._sire_object.property(mass_prop).value()
                    total_mass += mass
                except:
                    ValueError(
                        "Unable to compute center of mass. Missing 'mass' property!"
                    )

                # Get the coordinates and add to the reference using the
                # distance between them using the minimum image convention.
                try:
                    coord = atom._sire_object.property(coord_prop)
                    coord = ref_coord + _SireMaths.Vector(
                        space.calc_dist_vector(ref_coord, coord)
                    )
                except:
                    ValueError(
                        "Unable to compute center of mass. Missing 'coordinates' property!"
                    )

                # Update the center of mass.
                com += mass * coord

            return com, total_mass

        # Compute the center of mass for the atoms in all molecules in the system.
        if is_system:
            for mol in self:
                com, total_mass = update_com(mol.getAtoms(), com, total_mass)

        # Computer the center of mass for the atoms in this object.
        else:
            com, total_mass = update_com(self.getAtoms(), com, total_mass)

        # Normalise.
        com /= total_mass

        return [_Units.Length.angstrom * x.value() for x in com]

    def save(self, filebase):
        """
        Stream a wrapped Sire object to file.

        Parameters
        ----------

        sire_object : :class:`System <BioSimSpace._SireWrappers.SireWrapper>`
            A wrapped Sire object.

        filebase : str
            The base name of the binary output file.
        """
        from BioSimSpace.Stream import save as _save

        _save(self, filebase)

    def _getSireObject(self):
        """
        Return the underlying Sire object.

        Returns
        -------

        object : Sire.System.System, Sire.Mol.Molecule, Sire.Mol.Residue, Sire.Mol.Atom
            The Sire object that is being wrapped.
        """
        return self._sire_object

    def _getAABox(self, property_map={}):
        """
        Get the axis-aligned bounding box for the object.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        aabox : Sire.Vol.AABox
            The axis-aligned bounding box for the object.
        """
        from sire.legacy import Vol as _SireVol
        from .._Exceptions import IncompatibleError as _IncompatibleError
        from .. import _isVerbose

        # Initialise the coordinates vector.
        coord = []

        prop = property_map.get("coordinates", "coordinates")

        # Handle perturbable molecules.
        if self._is_perturbable:
            prop = "coordinates0"

        # Residues now have a coordinates property, but this is returned as a
        # Python list.
        try:
            c = self._sire_object.property(prop).to_vector()
        except:
            try:
                c = self._sire_object.property(prop)
            except Exception as e:
                msg = (
                    "Cannot compute axis-aligned bounding box "
                    + "since the object has no 'coordinates' property."
                )
                if _isVerbose():
                    print(msg)
                    raise _IncompatibleError(msg) from e
                else:
                    raise _IncompatibleError(msg) from None

        # We have a vector of coordinates. (Multiple atoms)
        if self._is_multi_atom:
            coord.extend(c)
        # Convert to a list.
        else:
            coord = [c]

        # Return the AABox for the coordinates.
        return _SireVol.AABox(coord)
