######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2023
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
A thin wrapper around Sire.Mol.Molecule. This is an internal package and should
not be directly exposed to the user.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Molecule"]

from math import isclose as _isclose
from warnings import warn as _warn

import numpy as _np
from sire.legacy import (
    Base as _SireBase,
    IO as _SireIO,
    MM as _SireMM,
    Maths as _SireMaths,
    Mol as _SireMol,
    System as _SireSystem,
    Units as _SireUnits,
)

from ._sire_wrapper import SireWrapper as _SireWrapper
from .. import _isVerbose
from ..Types import Coordinate as _Coordinate, Length as _Length, Vector as _BSSVector
from ..Units.Length import angstrom as _angstrom
from .._Exceptions import IncompatibleError as _IncompatibleError


class Molecule(_SireWrapper):
    """A container class for storing a molecule."""

    def __init__(self, molecule):
        """
        Constructor.

        Parameters
        ----------

        molecule : Sire.Mol.Molecule, :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            A Sire or BioSimSpace Molecule object.
        """

        # Set the force field variable. This records the force field with which
        # the molecule has been parameterised, i.e. by BSS.Parameters.
        self._forcefield = None

        # Set the components of the merged molecule to None.
        self._molecule0 = None
        self._molecule1 = None

        # Set the water model used to parameterise any structural ions
        # in the molecule.
        self._ion_water_model = None

        # Check that the molecule is valid.

        # A Sire Molecule object.
        if isinstance(molecule, _SireMol._Mol.Molecule):
            super().__init__(molecule)
            if self._sire_object.hasProperty("is_perturbable"):
                # Flag that the molecule is perturbable.
                self._is_perturbable = True

                # Extract the end states.
                if molecule.hasProperty("molecule0"):
                    self._molecule0 = Molecule(molecule.property("molecule0"))
                else:
                    self._molecule0, _ = self._extractMolecule()
                if molecule.hasProperty("molecule1"):
                    self._molecule1 = Molecule(molecule.property("molecule1"))
                else:
                    self._molecule1, _ = self._extractMolecule(is_lambda1=True)

        # Another BioSimSpace Molecule object.
        elif isinstance(molecule, Molecule):
            super().__init__(molecule._sire_object)
            if molecule._molecule0 is not None:
                self._molecule0 = Molecule(molecule._molecule0)
            if molecule._molecule1 is not None:
                self._molecule1 = Molecule(molecule._molecule1)
            self._forcefield = molecule._forcefield
            self._is_perturbable = molecule._is_perturbable

        # Invalid type.
        else:
            raise TypeError(
                f"'molecule' (type {type(molecule)}) must be of type 'Sire.Mol.Molecule' "
                "or 'BioSimSpace._SireWrappers.Molecule'."
            )

        # Flag that this object holds multiple atoms.
        self._is_multi_atom = True

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.Molecule: number=%d, nAtoms=%d, nResidues=%d>" % (
            self.number(),
            self.nAtoms(),
            self.nResidues(),
        )

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.Molecule: number=%d, nAtoms=%d, nResidues=%d>" % (
            self.number(),
            self.nAtoms(),
            self.nResidues(),
        )

    def __add__(self, other):
        """Addition operator."""

        # Convert tuple to a list.
        if isinstance(other, tuple):
            other = list(other)

        # Whether other is a container of molecules.
        is_sire_container = False

        # Validate the input.

        molecules = [self]

        # A System object.
        if isinstance(other, _System):
            system = _System(other)
            system.addMolecules(self)
            return system

        # A single Molecule object.
        elif isinstance(other, Molecule):
            molecules.append(other)

        # A Molecules object.
        elif isinstance(other, _Molecules):
            for molecule in other:
                molecules.append(molecule)

        # A list of Molecule objects.
        elif isinstance(other, list) and all(isinstance(x, Molecule) for x in other):
            for molecule in other:
                molecules.append(molecule)

        # Unsupported.
        else:
            raise TypeError(
                "'other' must be of type 'BioSimSpace._SireWrappers.System', "
                "'BioSimSpace._SireWrappers.Molecule', 'BioSimSpace._SireWrappers.Molecules' "
                "or a list of 'BioSimSpace._SireWrappers.Molecule' types"
            )

        # Create a new Molecules container.
        return _Molecules(molecules)

    def __contains__(self, other):
        """Return whether other is in self."""

        if not isinstance(other, (_Atom, _Residue)):
            raise TypeError(
                "'other' must be of type 'BioSimSpace._SireWrappers.Atom' "
                "or 'BioSimSpace._SireWrappers.Residue'."
            )

        # Return whether the object comes from this molecule.
        return self._sire_object.molecule() == other._sire_object.molecule()

    def copy(self):
        """
        Return a copy of this Molecule.

        Returns
        -------

        molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            A copy of the object.
        """
        # Copy the Sire object.
        mol = self._sire_object.__deepcopy__()

        # Give the molecule a unique number.
        mol = mol.edit().renumber(_SireMol.MolNum.getUniqueNumber()).commit().molecule()

        return Molecule(mol)

    def number(self):
        """
        Return the number of the molecule. Each molecule has a unique
        identification number.

        Returns
        -------

        mol_num : int
            The unique number of the molecule.
        """
        return self._sire_object.number().value()

    def coordinates(self, property_map={}):
        """
        Return the coordinates of the atoms in the molecule.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        [coordinates] : [class:`Coordinate <BioSimSpace.Types.Coordinate>`]
            The coordinates of the atoms in the molecule.
        """
        prop = property_map.get("coordinates", "coordinates")

        # Get the "coordinates" property from the molecule.
        try:
            sire_coord = self._sire_object.property(prop).toVector()
            coordinates = []
            for coord in sire_coord:
                coordinates.append(
                    _Coordinate(
                        _Length(coord[0], "Angstrom"),
                        _Length(coord[1], "Angstrom"),
                        _Length(coord[2], "Angstrom"),
                    )
                )
        except:
            return None

        # Return the coordinates.
        return coordinates

    def getResidues(self):
        """
        Return a list containing the residues in the molecule.

        Returns
        -------

        residues : [:class:`Residue <BioSimSpace._SireWrappers.Residue>`]
            The list of residues in the molecule.
        """
        residues = []
        for residue in self._sire_object.residues():
            residues.append(_Residue(residue))
        return residues

    def getAtoms(self):
        """
        Return a list containing the atoms in the molecule.

        Returns
        -------

        atoms : [:class:`Atom <BioSimSpace._SireWrappers.Atom>`]
            The list of atoms in the molecule.
        """
        atoms = []
        for atom in self._sire_object.atoms():
            atoms.append(_Atom(atom))
        return atoms

    def extract(self, indices, renumber=False, property_map={}):
        """
        Extract atoms at the specified indices from the molecule to create
        a new molecule.

        Parameters
        ----------

        indices : [ int ]
            The indices of the atoms to extract.

        renumber : bool
            Whether the returned molecule has the same number as the original.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            The extracted molecule.
        """

        # TODO: This method is slow for large molecules. Re-write in pure C++
        # and provide a suitable wrapper function.

        # Convert tuple to list.
        if isinstance(indices, tuple):
            indices = list(indices)

        if not isinstance(indices, list):
            raise TypeError("'indices' must be a list of 'int' types.")

        # A list to store the indices, mapped back to 0 --> nAtoms() - 1.
        indices_ = []

        # Store the maximum allowed index.
        max_index = self.nAtoms() - 1

        for x in indices:
            # Check type.
            if not type(x) is int:
                raise TypeError("'indices' must be a list of 'int' types.")

            # Map index back to range.
            if x < 0:
                x += self.nAtoms()

            # Check value.
            if x < 0 or x > max_index:
                raise TypeError(f"'indices' must be in range [0, {max_index}].")

            # Append to the updated atom index.
            indices_.append(_SireMol.AtomIdx(x))

        # Extract a partial molecule.
        try:
            # Create an empty atom selection for this molecule.
            selection = self._sire_object.selection()
            selection.selectNone()

            # Add the atom indices to the selection.
            for idx in indices_:
                selection.select(idx)

            partial_mol = (
                _SireMol.PartialMolecule(self._sire_object, selection)
                .extract()
                .molecule()
            )
        except Exception as e:
            msg = "Unable to create partial molecule!"
            if _isVerbose():
                raise _IncompatibleError(msg) from e
            else:
                raise _IncompatibleError(msg) from None

        # Get the "intrascale" property name.
        intrascale = property_map.get("intrascale", "intrascale")

        # Flag whether the molecule has an intrascale property.
        has_intrascale = self._sire_object.hasProperty(intrascale)

        # Remove the "intrascale" property, since this doesn't correspond to the
        # extracted molecule.
        if has_intrascale:
            partial_mol = (
                partial_mol.edit().removeProperty(intrascale).molecule().commit()
            )

            # Recreate the molecule.
            mol = Molecule(partial_mol)

            # Now parse the molecule as a GROMACS topology to recreate the
            # intrascale matrix.
            try:
                gro_top = _SireIO.GroTop(mol.toSystem()._sire_object)
            except Exception as e:
                msg = "Unable to recover non-bonded matrix for the extracted molecule!"
                if _isVerbose():
                    raise _IncompatibleError(msg) from e
                else:
                    raise _IncompatibleError(msg) from None

            # Convert back to a Sire system.
            gro_sys = gro_top.toSystem()

            # Add the intrascale property back into the molecule.
            edit_mol = mol._sire_object.edit()
            edit_mol.setProperty(
                intrascale, gro_sys[_SireMol.MolIdx(0)].property("intrascale")
            )

            # Recreate the molecule.
            mol = Molecule(edit_mol.commit())

        else:
            mol = Molecule(partial_mol)

        # Keep the same MolNum.
        if not renumber:
            mol._sire_object = (
                mol._sire_object.edit()
                .renumber(self._sire_object.number())
                .commit()
                .molecule()
            )

        return mol

    def molecule0(self):
        """
        Return the component of the merged molecule at lambda = 0.

        Returns
        -------

        molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            The component of the merged molecule at lambda = 0.
            Returns None if this isn't a merged molecule.
        """
        return self._molecule0

    def molecule1(self):
        """
        Return the component of the merged molecule at lambda = 1.

        Returns
        -------

        molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            The component of the merged molecule at lambda = 1.
            Returns None if this isn't a merged molecule.
        """
        return self._molecule1

    def nAtoms(self):
        """
        Return the number of atoms in the molecule.

        Returns
        -------

        num_atoms : int
            The number of atoms in the molecule.
        """
        return self._sire_object.nAtoms()

    def nResidues(self):
        """
        Return the number of residues in the molecule.

        Returns
        -------

        num_residues : int
            The number of residues in the molecule.
        """
        return self._sire_object.nResidues()

    def nChains(self):
        """
        Return the number of chains in the molecule.

        Returns
        -------

        num_chains : int
            The number of chains in the molecule.
        """
        return self._sire_object.nChains()

    def isPerturbable(self):
        """
        Whether this molecule is perturbable, i.e. it can be used in a
        free-energy perturbation simulation.

        Returns
        -------

        is_perturbable : bool
            Whether the molecule is perturbable.
        """
        return self._is_perturbable

    def isDecoupled(self):
        """
        Whether this molecule is decoupled, i.e. it can be used in a
        free-energy decoupling simulation.

        Returns
        -------

        is_decoupled : bool
            Whether the molecule is decoupled.
        """
        if self._sire_object.hasProperty("decouple"):
            return True
        else:
            return False

    def isML(self):
        """
        Whether this molecule is marked as ML molecule, i.e. it can be used in a
        MM to ML transformation simulation.

        Returns
        -------

        is_ML : bool
            Whether the molecule is marked as ML.
        """
        if self._sire_object.hasProperty("ML"):
            return True
        else:
            return False

    def isAlchemicalIon(self):
        """
        Whether this molecule is marked as Alchemical Ion.

        Returns
        -------

        isAlchemicalIon : bool
            Whether the molecule is marked as Alchemical Ion.
        """
        return self._sire_object.hasProperty("AlchemicalIon")

    def isWater(self, property_map={}):
        """
        Whether this is a water molecule.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        is_water : bool
            Whether this is a water molecule.
        """

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        return _SireIO.isWater(self._sire_object, property_map)

    def isAmberWater(self, property_map={}):
        """
        Whether this is an AMBER format water molecule.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        is_amber_water : bool
            Whether this molecule is an AMBER format water.
        """

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        return _SireIO.isAmberWater(self._sire_object, property_map)

    def isGromacsWater(self, property_map={}):
        """
        Whether this is a GROMACS format water molecule.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their

        Returns
        -------

        is_gromacs_water : bool
            Whether this molecule is a GROMACS format water.
        """

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        return _SireIO.isGromacsWater(self._sire_object, property_map)

    def toSystem(self):
        """
        Convert a single Molecule to a System.

        Returns
        -------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
        """
        return _System(self)

    def search(self, query, property_map={}):
        """
        Search the molecule for atoms and residues. Search results will be
        reduced to their minimal representation, i.e. a residue containing
        a single atom will be returned as a atom.

        Parameters
        ----------

        query : str
            The search query.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        results : [:class:`Atom <BioSimSpace._SireWrappers.Atom>`, \
                   :class:`Residue <BioSimSpace._SireWrappers.Residue>`, ...]
            A list of objects matching the search query.

        Examples
        --------

        Search for all residues named ALA.

        >>> result = molecule.search("resname ALA")

        Search for all oxygen or hydrogen atoms.

        >>> result = molecule.search("element oxygen or element hydrogen")

        Search for atom index 23.

        >>> result = molecule.search("atomidx 23")
        """

        if not isinstance(query, str):
            raise TypeError("'query' must be of type 'str'")

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        # Initialise a list to hold the search results.
        results = []

        try:
            query = _SireMol.Select(query)
        except Exception as e:
            msg = "'Invalid search query: %r (%s)" % (query, e)
            if _isVerbose():
                raise ValueError(msg) from e
            else:
                raise ValueError(msg) from None

        try:
            # Query the Sire molecule.
            search_result = query(self._sire_object, property_map)

        except Exception as e:
            msg = "'Invalid search query: %r : %s" % (query, e)
            if _isVerbose():
                raise ValueError(msg) from e
            else:
                raise ValueError(msg) from None

        return _SearchResult(search_result)

    def makeCompatibleWith(
        self,
        molecule,
        property_map={},
        overwrite=True,
        rename_atoms=False,
        verbose=False,
    ):
        """
        Make this molecule compatible with passed one, i.e. match atoms and
        add all additional properties from the passed molecule while preserving
        the topology and naming/numbering convention of this molecule.

        Parameters
        ----------

        molecule : Sire.Mol.Molecule, :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                   Sire.System.System, :class:`System <Sire._SireWrappers.System>`
            The molecule, or system of molecules, to match with.

        property_map : dict
            A map between property names and user supplied names.

        overwrite : bool
            Whether to overwrite any duplicate properties.

        rename_atoms : bool
            Whether to rename atoms if they have changed.

        verbose : bool
            Whether to report status updates to stdout.
        """

        # Validate input.

        is_system = False
        if isinstance(molecule, _SireMol.Molecule):
            mol1 = molecule
        elif isinstance(molecule, Molecule):
            mol1 = molecule._sire_object
        elif isinstance(molecule, _SireSystem.System):
            mol1 = molecule
            is_system = True
        elif isinstance(molecule, _System):
            mol1 = molecule._sire_object
            is_system = True
        else:
            raise TypeError(
                "'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule', or 'Sire.Mol.Molecule'"
            )

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        if not isinstance(overwrite, bool):
            raise TypeError("'overwrite' must be of type 'bool'")

        if not isinstance(rename_atoms, bool):
            raise TypeError("'rename_atoms' must be of type 'bool'")

        if not isinstance(verbose, bool):
            raise TypeError("'verbose' must be of type 'bool'")

        # Get the two Sire molecules.
        mol0 = self._sire_object

        # Store the number of atoms to match.
        num_atoms0 = mol0.nAtoms()

        # Work out the number of atoms in mol1.
        if is_system:
            num_atoms1 = _System(mol1).nAtoms()
        else:
            num_atoms1 = mol1.nAtoms()

        # The new molecule must have the same number of atoms.
        if num_atoms1 != num_atoms0:
            if is_system:
                raise _IncompatibleError(
                    "The passed system is incompatible with the original! "
                    "self.nAtoms() = %d, other.nAtoms() = %d" % (num_atoms0, num_atoms1)
                )
            else:
                raise _IncompatibleError(
                    "The passed molecule is incompatible with the original! "
                    "self.nAtoms() = %d, other.nAtoms() = %d" % (num_atoms0, num_atoms1)
                )

        # Whether the atoms have been renamed.
        is_renamed = False

        if not is_system:
            # Instantiate the default atom matcher (match by residue index and atom name).
            matcher = _SireMol.ResIdxAtomNameMatcher()

            # Match the atoms based on residue index and atom name.
            matches = matcher.match(mol0, mol1)

            # Have we matched all of the atoms?
            if len(matches) < num_atoms0:
                # Atom names might have changed. Try to match by residue index
                # and coordinates.
                matcher = _SireMol.ResIdxAtomCoordMatcher()
                matches = matcher.match(mol0, mol1)

                # We need to rename the atoms.
                is_renamed = True

                # Have we matched all of the atoms?
                if len(matches) < num_atoms0:
                    raise _IncompatibleError("Failed to match all atoms!")

                # Are the atoms in the same order?
                is_reordered = matcher.changesOrder(mol0, mol1)

            else:
                # Are the atoms in the same order?
                is_reordered = matcher.changesOrder(mol0, mol1)

            if verbose:
                print(
                    "\nAtom matching successful.\nAtom indices %s reordered."
                    % ("" if is_reordered else "not")
                )

            # Get a list of the property keys for each molecule.
            props0 = mol0.propertyKeys()
            props1 = mol1.propertyKeys()

            # Copy the property map.
            _property_map = property_map.copy()

            # See if any of the new properties are in the map, add them if not.
            for prop in props0:
                if not prop in _property_map:
                    _property_map[prop] = prop
            for prop in props1:
                if not prop in _property_map:
                    _property_map[prop] = prop

            # Make the molecule editable.
            edit_mol = mol0.edit()

            if "parameters" in _property_map:
                param = _property_map["parameters"]
            else:
                param = "parameters"

            # The atom order is the same, simply copy across properties as is.
            if not is_reordered:
                if verbose:
                    print("\nSetting properties...")
                # Loop over all of the keys in the new molecule.
                for prop in props1:
                    # Skip 'parameters' property, since it contains references to other parameters.
                    if prop != param:
                        # This is a new property, or we are allowed to overwrite.
                        if (not mol0.hasProperty(_property_map[prop])) or overwrite:
                            if verbose:
                                print("  %s" % _property_map[prop])
                            try:
                                edit_mol = edit_mol.setProperty(
                                    _property_map[prop], mol1.property(prop)
                                )
                            except Exception as e:
                                msg = (
                                    "Failed to set property '%s'" % _property_map[prop]
                                )
                                if _isVerbose():
                                    raise _IncompatibleError(msg) from e
                                else:
                                    raise _IncompatibleError(msg) from None

            # The atom order is different, we need to map the atoms when setting properties.
            else:
                # Create a dictionary to flag whether a property has been seen.
                seen_prop = {}
                for prop in props1:
                    seen_prop[prop] = False

                # First, set atom based properties.

                if verbose:
                    print("\nSetting atom properties...")

                # Loop over all of the keys in the new molecule.
                for prop in props1:
                    # This is a new property, or we are allowed to overwrite.
                    if (not mol0.hasProperty(_property_map[prop])) or overwrite:
                        # Loop over all of the atom mapping pairs and set the property.
                        for idx0, idx1 in matches.items():
                            # Does the atom have this property?
                            # If so, add it to the matching atom in this molecule.
                            if mol1.atom(idx1).hasProperty(prop):
                                if verbose:
                                    print(
                                        "  %-20s %s --> %s"
                                        % (_property_map[prop], idx1, idx0)
                                    )
                                try:
                                    edit_mol = (
                                        edit_mol.atom(idx0)
                                        .setProperty(
                                            _property_map[prop],
                                            mol1.atom(idx1).property(prop),
                                        )
                                        .molecule()
                                    )
                                    seen_prop[prop] = True
                                except Exception as e:
                                    msg = (
                                        "Failed to copy property '%s' from %s to %s."
                                        % (_property_map[prop], idx1, idx0)
                                    )
                                    if _isVerbose():
                                        raise _IncompatibleError(msg) from e
                                    else:
                                        raise _IncompatibleError(msg) from None

                # Now deal with all unseen properties. These will be non atom-based
                # properties, such as TwoAtomFunctions, StringProperty, etc.

                if verbose:
                    print("\nSetting molecule properties...")

                # Create a dictionary mapping the atom matches from mol1 to mol0.
                inv_matches = {}
                for key, value in matches.items():
                    inv_matches[value] = key

                # Loop over all of the unseen properties.
                for prop in seen_prop:
                    if not seen_prop[prop]:
                        # Skip 'parameters' property, since it contains references to other parameters.
                        if prop != "parameters":
                            # This is a new property, or we are allowed to overwrite.
                            if (not mol0.hasProperty(_property_map[prop])) or overwrite:
                                if verbose:
                                    print("  %s" % _property_map[prop])

                                # Get the property from the parameterised molecule.
                                propty = mol1.property(_property_map[prop])

                                # Try making it compatible with the original molecule.
                                if hasattr(propty, "makeCompatibleWith"):
                                    try:
                                        propty = propty.makeCompatibleWith(
                                            mol0, inv_matches
                                        )
                                    except Exception as e:
                                        msg = (
                                            "Incompatible property: %s"
                                            % _property_map[prop]
                                        )
                                        if _isVerbose():
                                            raise _IncompatibleError(msg) from e
                                        else:
                                            raise _IncompatibleError(msg) from None

                                # Now try to set the property.
                                edit_mol.setProperty(_property_map[prop], propty)

            # Finally, rename the atoms.

            if rename_atoms and is_renamed:
                if verbose:
                    print("\nRenaming atoms...")

                for idx0, idx1 in matches.items():
                    # Get the name of the atom in each molecule.
                    name0 = mol0.atom(idx0).name()
                    name1 = mol1.atom(idx1).name()

                    if verbose:
                        print("  %s --> %s" % (name0, name1))

                    # Try to rename the atom.
                    try:
                        edit_mol = (
                            edit_mol.atom(idx0)
                            .rename(mol1.atom(idx1).name())
                            .molecule()
                        )
                    except Exception as e:
                        msg = "Failed to rename atom: %s --> %s" % (name0, name1)
                        if _isVerbose():
                            raise _IncompatibleError(msg) from e
                        else:
                            raise _IncompatibleError(msg) from None

            # Commit the changes.
            self._sire_object = edit_mol.commit()

        # Atoms from molecule in the passed system (mol1) need to be matched against atoms
        # from the corresponding residue in mol0 and the properties aggregated.
        else:
            # Initialise a list to hold the matches for each molecule in mol1.
            matches = []

            # Tally counter for the total number of matches.
            num_matches = 0

            # Initialise the offset.
            offset = 0

            # Get the molecule numbers in the system.
            mol_nums = mol1.molNums()

            # Loop over all molecules in mol1.
            for num in mol_nums:
                # Extract the numbered molecule from the system mol1.
                mol = mol1[num]

                # Initialise the matcher.
                matcher = _SireMol.ResIdxAtomCoordMatcher(_SireMol.ResIdx(offset))

                # Get the matches for this molecule and append to the list.
                match = matcher.match(mol0, mol)
                matches.append(match)
                num_matches += len(match)

                # Increment the offset.
                offset += mol.nResidues()

            # Have we matched all of the atoms?
            if num_matches < num_atoms0:
                raise _IncompatibleError("Failed to match all atoms!")

            # Make the molecule editable.
            edit_mol = mol0.edit()

            # Create objects to hold all of the potential terms.
            bonds = _SireMM.TwoAtomFunctions(edit_mol.info())
            angles = _SireMM.ThreeAtomFunctions(edit_mol.info())
            dihedrals = _SireMM.FourAtomFunctions(edit_mol.info())
            impropers = _SireMM.FourAtomFunctions(edit_mol.info())

            # Next we need to work out what properties can be set at the atom level,
            # and which are associated with the molecule as a whole.
            mol_props = []
            atom_props = []

            # Each atom should have the same set of properties so we can check the
            # first atom in each molecule.
            for num in mol_nums:
                mol = mol1[num]

                # Get the molecule and atom properties.
                props_mol = mol.propertyKeys()
                props_atom = mol.atoms()[0].propertyKeys()

                # Check the atomic properties and add any new ones to the list.
                for prop in props_atom:
                    prop = property_map.get(prop, prop)
                    if prop not in atom_props:
                        atom_props.append(prop)

                # Check the molecular properties and add any new ones to the list.
                for prop in props_mol:
                    prop = property_map.get(prop, prop)
                    if prop not in mol_props and prop not in atom_props:
                        mol_props.append(prop)

            # Create a list of excluded molecular properties. These are ones that
            # must be re-mapped manually and set by hand.
            excluded_props = [
                property_map.get("bond", "bond"),
                property_map.get("angle", "angle"),
                property_map.get("dihedral", "dihedral"),
                property_map.get("improper", "improper"),
                property_map.get("connectivity", "connectivity"),
                property_map.get("intrascale", "intrascale"),
                property_map.get("parameters", "parameters"),
            ]

            # Loop over all atoms within each molecule. We set the allowed atomic
            # properties and build the molecular ones as we go.
            for idx, num in enumerate(mol_nums):
                # Extract the molecule and its associated info object.
                mol = mol1[num]
                info = mol.info()

                # An inverse mapping for indices in mol1 to those in mol0.
                inv_matches = {}

                # Loop over all matching atom pairs for this molecule.
                for idx0, idx1 in matches[idx].items():
                    # Add indices to the inverse mapping.
                    inv_matches[idx1] = idx0

                    # Check the atom names and see if they need updating.
                    if rename_atoms:
                        name0 = mol0.atom(idx0).name()
                        name1 = mol1.atom(idx1).name()

                        if verbose:
                            print("  %s --> %s" % (name0, name1))

                        # Try to rename the atom.
                        try:
                            edit_mol = (
                                edit_mol.atom(idx0)
                                .rename(mol1.atom(idx1).name())
                                .molecule()
                            )
                        except Exception as e:
                            msg = "Failed to rename atom: %s --> %s" % (name0, name1)
                            if _isVerbose():
                                raise _IncompatibleError(msg) from e
                            else:
                                raise _IncompatibleError(msg) from None

                    # Loop over all atom properties.
                    for prop in atom_props:
                        # This is a new property, or we're allowed to overwrite.
                        if (not mol0.hasProperty(prop)) or overwrite:
                            if verbose:
                                print("  %-20s %s --> %s" % (prop, idx1, idx0))
                            try:
                                edit_mol = (
                                    edit_mol.atom(idx0)
                                    .setProperty(prop, mol.atom(idx1).property(prop))
                                    .molecule()
                                )
                            except Exception as e:
                                msg = "Failed to copy property '%s' from %s to %s." % (
                                    prop,
                                    idx1,
                                    idx0,
                                )
                                if _isVerbose():
                                    raise _IncompatibleError(msg) from e
                                else:
                                    raise _IncompatibleError(msg) from None

                # Now deal with the molecular properties.
                for prop in mol_props:
                    # Get the property name from the user mapping.
                    prop = property_map.get(prop, prop)

                    # This is a new property, or we're allowed to overwrite, and it's not excluded.
                    if (
                        (not mol0.hasProperty(prop)) or overwrite
                    ) and prop not in excluded_props:
                        if verbose:
                            print("  %s" % prop)

                        # Get the property from the parameterised molecule.
                        propty = mol.property(prop)

                        # Try making the property compatible with the original molecule.
                        if hasattr(propty, "makeCompatibleWith"):
                            try:
                                propty = propty.makeCompatibleWith(mol0, inv_matches)
                            except Exception as e:
                                msg = "Incompatible property: %s" % prop
                                if _isVerbose():
                                    raise _IncompatibleError(msg) from e
                                else:
                                    raise _IncompatibleError(msg) from None

                        # Now try to set the property.
                        edit_mol.setProperty(prop, propty)

                # Now re-map and build the properties for each of the potential terms.

                # Bonds.
                prop = property_map.get("bond", "bond")
                if mol.hasProperty(prop):
                    if verbose:
                        print("  %s" % prop)
                    for bond in mol.property(prop).potentials():
                        # Extract the bond information.
                        atom0 = info.atomIdx(bond.atom0())
                        atom1 = info.atomIdx(bond.atom1())
                        exprn = bond.function()

                        # Map the atom indices to their position in the merged molecule.
                        atom0 = inv_matches[atom0]
                        atom1 = inv_matches[atom1]

                        # Set the new bond.
                        bonds.set(atom0, atom1, exprn)

                # Angles.
                prop = property_map.get("angle", "angle")
                if mol.hasProperty(prop):
                    if verbose:
                        print("  %s" % prop)
                    for angle in mol.property(prop).potentials():
                        # Extract the angle information.
                        atom0 = info.atomIdx(angle.atom0())
                        atom1 = info.atomIdx(angle.atom1())
                        atom2 = info.atomIdx(angle.atom2())
                        exprn = angle.function()

                        # Map the atom indices to their position in the merged molecule.
                        atom0 = inv_matches[atom0]
                        atom1 = inv_matches[atom1]
                        atom2 = inv_matches[atom2]

                        # Set the new angle.
                        angles.set(atom0, atom1, atom2, exprn)

                # Dihedrals.
                prop = property_map.get("dihedral", "dihedral")
                if mol.hasProperty(prop):
                    if verbose:
                        print("  %s" % prop)
                    for dihedral in mol.property(prop).potentials():
                        # Extract the dihedral information.
                        atom0 = info.atomIdx(dihedral.atom0())
                        atom1 = info.atomIdx(dihedral.atom1())
                        atom2 = info.atomIdx(dihedral.atom2())
                        atom3 = info.atomIdx(dihedral.atom3())
                        exprn = dihedral.function()

                        # Map the atom indices to their position in the merged molecule.
                        atom0 = inv_matches[atom0]
                        atom1 = inv_matches[atom1]
                        atom2 = inv_matches[atom2]
                        atom3 = inv_matches[atom3]

                        # Set the new dihedral.
                        dihedrals.set(atom0, atom1, atom2, atom3, exprn)

                # Impropers.
                prop = property_map.get("improper", "improper")
                if mol.hasProperty(prop):
                    if verbose:
                        print("  %s" % prop)
                    for improper in mol.property(prop).potentials():
                        # Extract the improper information.
                        atom0 = info.atomIdx(improper.atom0())
                        atom1 = info.atomIdx(improper.atom1())
                        atom2 = info.atomIdx(improper.atom2())
                        atom3 = info.atomIdx(improper.atom3())
                        exprn = improper.function()

                        # Map the atom indices to their position in the merged molecule.
                        atom0 = inv_matches[atom0]
                        atom1 = inv_matches[atom1]
                        atom2 = inv_matches[atom2]
                        atom3 = inv_matches[atom3]

                        # Set the new improper.
                        impropers.set(atom0, atom1, atom2, atom3, exprn)

            # Set properties for the molecular potential.

            # Bonds.
            if bonds.nFunctions() > 0:
                prop = property_map.get("bond", "bond")
                if verbose:
                    print("  %s" % prop)
                try:
                    edit_mol.setProperty(prop, bonds)
                except Exception as e:
                    msg = "Incompatible property: %s" % prop
                    if _isVerbose():
                        raise _IncompatibleError(msg) from e
                    else:
                        raise _IncompatibleError(msg) from None

            # Angles.
            if angles.nFunctions() > 0:
                prop = property_map.get("angle", "angle")
                if verbose:
                    print("  %s" % prop)
                try:
                    edit_mol.setProperty(prop, angles)
                except Exception as e:
                    msg = "Incompatible property: %s" % prop
                    if _isVerbose():
                        raise _IncompatibleError(msg) from e
                    else:
                        raise _IncompatibleError(msg) from None

            # Dihedrals.
            if dihedrals.nFunctions() > 0:
                prop = property_map.get("dihedral", "dihedral")
                if verbose:
                    print("  %s" % prop)
                try:
                    edit_mol.setProperty(prop, dihedrals)
                except Exception as e:
                    msg = "Incompatible property: %s" % prop
                    if _isVerbose():
                        raise _IncompatibleError(msg) from e
                    else:
                        raise _IncompatibleError(msg) from None

            # Impropers.
            if impropers.nFunctions() > 0:
                prop = property_map.get("improper", "improper")
                if verbose:
                    print("  %s" % prop)
                try:
                    edit_mol.setProperty(prop, impropers)
                except Exception as e:
                    msg = "Incompatible property: %s" % prop
                    if _isVerbose():
                        raise _IncompatibleError(msg) from e
                    else:
                        raise _IncompatibleError(msg) from None

            # Now generate the molecular connectivity.
            if bonds.nFunctions() > 0:
                prop = property_map.get("connectivity", "connectivity")
                if verbose:
                    print("  %s" % prop)
                conn = _SireMol.Connectivity(edit_mol.info()).edit()

                # Connect the bonded atoms. Connectivity is the same at lambda = 0
                # and lambda = 1.
                for bond in bonds.potentials():
                    conn.connect(bond.atom0(), bond.atom1())
                conn = conn.commit()

                try:
                    edit_mol.setProperty(prop, conn)
                except Exception as e:
                    msg = "Incompatible property: %s" % prop
                    if _isVerbose():
                        raise _IncompatibleError(msg) from e
                    else:
                        raise _IncompatibleError(msg) from None

            # Next we construct the intrascale matrix for the non-bonded
            # interactions. It is prohibitively slow to do this on-the-fly so
            # we use the GroTop parser to re-construct it for us, then copy it
            # back into the original system.

            prop = property_map.get("intrascale", "intrascale")
            if (not mol0.hasProperty(prop)) or overwrite:
                if verbose:
                    print("  %s" % prop)
                # Delete any existing intrascale property from the molecule.
                if mol0.hasProperty(prop):
                    edit_mol.removeProperty(prop)
                mol = edit_mol.commit()
                # Convert to a "GROMACS system" using the GroTop parser.
                gro_system = _SireIO.GroTop(
                    Molecule(mol).toSystem()._sire_object,
                    _SireBase.PropertyMap(property_map),
                ).toSystem()
                # Extract the only molecule in the system.
                gro_mol = gro_system[_SireMol.MolIdx(0)]
                edit_mol = mol.edit()
                try:
                    edit_mol.setProperty(prop, gro_mol.property(prop))
                except Exception as e:
                    msg = "Incompatible property: %s" % prop
                    if _isVerbose():
                        raise _IncompatibleError(msg) from e
                    else:
                        raise _IncompatibleError(msg) from None

            # Finally, commit the changes to the internal object.
            self._sire_object = edit_mol.commit()

    def translate(self, vector, property_map={}):
        """
        Translate each molecule in the container.

        Parameters
        ----------

        vector : [:class:`Length <BioSimSpace.Types.Length>`]
            The translation vector in Angstroms.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Convert tuple to a list.
        if isinstance(vector, tuple):
            vector = list(vector)

        # Validate input.
        if isinstance(vector, list):
            if len(vector) != 3:
                raise ValueError(
                    "'vector' must contain 3 items, i.e. x, y, z components!"
                )
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

        # Translate the molecule.
        # Copy the property map.
        _property_map = property_map.copy()

        # If this is a perturbable molecule then translate both lambda end states.
        if self._is_perturbable:
            # lambda = 0
            _property_map["coordinates"] = "coordinates0"
            mol = (
                self._sire_object.move()
                .translate(_SireMaths.Vector(vec), _property_map)
                .commit()
            )

            # lambda = 1
            _property_map["coordinates"] = "coordinates1"
            mol = mol.move().translate(_SireMaths.Vector(vec), _property_map).commit()

        else:
            mol = (
                self._sire_object.move()
                .translate(_SireMaths.Vector(vec), _property_map)
                .commit()
            )

        # Update the internal molecule object.
        self._sire_object = mol

    def repartitionHydrogenMass(
        self, factor=4, water="no", use_coordinates=False, property_map={}
    ):
        """
        Redistrubute mass of heavy atoms connected to bonded hydrogens into
        the hydrogen atoms. This allows the use of larger simulation
        integration time steps without encountering instabilities related
        to high-frequency hydrogen motion.

        Parameters
        ----------

        factor : float
            The repartioning scale factor. Hydrogen masses are scaled by this
            amount.

        water : str
            Whether to repartition masses for water molecules. Options are
            "yes", "no", and "exclusive", which can be used to repartition
            masses for water molecules only.

        use_coordinates : bool
            Whether to use the current molecular coordinates to work out
            the connectivity before repartitioning. If False, the information
            from the molecular topology, e.g. force field, will be used, if
            present.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Convert int to float.
        if type(factor) is int:
            factor = float(factor)

        # Check scale factor.
        if not isinstance(factor, float):
            raise TypeError("'factor' must be of type 'float'.")
        if factor <= 0:
            raise ValueError("'factor' must be positive!")

        # Check water handling.
        if not isinstance(water, str):
            raise TypeError("'water' must be of type 'str'.")

        # Strip whitespace and convert to lower case.
        water = water.replace(" ", "").lower()

        # Allowed options and mapping to Sire flag.
        water_options = {"no": 0, "yes": 1, "exclusive": 2}

        if water not in water_options:
            water_string = ", ".join(f"'{x}'" for x in water_options)
            raise ValueError(f"'water' must be one of: {water_string}")

        if not isinstance(use_coordinates, bool):
            raise TypeError("'use_coordinates' must be of type 'bool'.")

        # Check property map.
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'.")

        # Update the property map to indicate that coordinates will be used
        # to compute the connectivity.
        pmap = property_map.copy()
        if use_coordinates:
            pmap["use_coordinates"] = _SireBase.wrap(True)

        # Handle perturbable molecules separately.
        if self.isPerturbable():
            # Search for dummies in both end states.
            try:
                dummies0 = self.search(
                    "element Xx", property_map={"element": "element0"}
                )
            except:
                dummies0 = []
            try:
                dummies1 = self.search(
                    "element Xx", property_map={"element": "element1"}
                )
            except:
                dummies1 = []

            # Repartition masses for the lambda=0 state.
            pmap = {
                "mass": "mass0",
                "element": "element0",
                "coordinates": "coordinates0",
            }
            self._sire_object = _SireIO.repartitionHydrogenMass(
                self._sire_object, factor, water_options[water], pmap
            )

            # Repartition masses for the lambda=1 state.
            pmap = {
                "mass": "mass1",
                "element": "element1",
                "coordinates": "coordinates1",
            }
            self._sire_object = _SireIO.repartitionHydrogenMass(
                self._sire_object, factor, water_options[water], pmap
            )

            # Now replace atom dummy atom masses with the reparitioned mass
            # from the opposite end state.

            edit_mol = self._sire_object.edit()

            for dummy in dummies0:
                idx = dummy._sire_object.index()
                mass1 = self._sire_object.atom(idx).property("mass1")
                edit_mol = edit_mol.atom(idx).setProperty("mass0", mass1).molecule()
            for dummy in dummies1:
                idx = dummy._sire_object.index()
                mass0 = self._sire_object.atom(idx).property("mass0")
                edit_mol = edit_mol.atom(idx).setProperty("mass1", mass0).molecule()

            self._sire_object = edit_mol.commit()

        else:
            self._sire_object = _SireIO.repartitionHydrogenMass(
                self._sire_object, factor, water_options[water], pmap
            )

    def _getPropertyMap0(self):
        """Generate a property map for the lambda = 0 state of the merged molecule."""

        property_map = {}

        if self._is_perturbable:
            for prop in self._sire_object.propertyKeys():
                if prop[-1] == "0":
                    property_map[prop[:-1]] = prop

        return property_map

    def _getPropertyMap1(self):
        """Generate a property map for the lambda = 1 state of the merged molecule."""

        property_map = {}

        if self._is_perturbable:
            for prop in self._sire_object.propertyKeys():
                if prop[-1] == "1":
                    property_map[prop[:-1]] = prop

        return property_map

    def _fixCharge(self, property_map={}):
        """
        Make the molecular charge an integer value.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps "properties" in this molecule to their
            user defined values. This allows the user to refer to properties
            with their own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Get the user defined charge property.
        prop = property_map.get("charge", "charge")

        if not self._sire_object.hasProperty(prop):
            raise _IncompatibleError(
                "Molecule does not have charge property: '%s'." % prop
            )

        # Calculate the charge.
        charge = self.charge(property_map=property_map).value()

        # Calculate the difference from the nearest integer value.
        delta = round(charge) - charge

        # The difference is too small to care about.
        if _isclose(charge + delta, charge, rel_tol=1e-6):
            return

        # Normalise by the number of atoms.
        delta /= self.nAtoms()

        # Make the molecule editable.
        edit_mol = self._sire_object.edit()

        # Shift the charge of each atom in the molecule by delta.
        # Make sure to invert the sign of the charge since it is in
        # units of -1 |e|.
        for atom in edit_mol.atoms():
            charge = edit_mol.atom(atom.index()).property(prop).value()
            charge = -(charge + delta)
            edit_mol = (
                edit_mol.atom(atom.index())
                .setProperty(prop, charge * _SireUnits.e_charge)
                .molecule()
            )

        # Update the Sire molecule.
        self._sire_object = edit_mol.commit()

    def _toRegularMolecule(
        self,
        property_map={},
        is_lambda1=False,
        convert_amber_dummies=False,
        generate_intrascale=False,
    ):
        """
        Internal function to convert a merged molecule to a regular molecule.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        is_lambda1 : bool
            Whether to use the molecule at the lambda = 1 end state.
            By default, the state at lambda = 0 is used.

        convert_amber_dummies : bool
            Whether to convert dummies to the correct AMBER formatting for
            non-FEP simulations. This will replace the "du" ambertype
            and "Xx" element with the properties from the other end state.

        generate_intrascale : bool
            Whether to regenerate the intrascale matrix.

        Returns
        -------

        molecule : BioSimSpace._SireWrappers.Molecule
            The molecule at the chosen end state.
        """

        if not isinstance(is_lambda1, bool):
            raise TypeError("'is_lambda1' must be of type 'bool'")

        if not isinstance(convert_amber_dummies, bool):
            raise TypeError("'convert_amber_dummies' must be of type 'bool'")

        if is_lambda1:
            lam = "1"
        else:
            lam = "0"

        if not self._is_perturbable:
            return Molecule(self._sire_object)

        # Extract and copy the Sire molecule.
        mol = self._sire_object.__deepcopy__()

        # Make the molecule editable.
        mol = mol.edit()

        # Remove the perturbable molecule flag.
        mol = mol.removeProperty("is_perturbable").molecule()

        # Flag that the molecule was perturbable, so that dummies should be
        # treated as "normal" atoms.
        mol = mol.setProperty("was_perturbable", _SireBase.wrap(True)).molecule()

        # Rename all properties in the molecule for the corresponding end state,
        # e.g.: "prop0" --> "prop". Then delete all properties named "prop0"
        # and "prop1".
        for prop in mol.propertyKeys():
            if prop[-1] == lam:
                # See if this property exists in the user map.
                new_prop = property_map.get(prop[:-1], prop[:-1])

                # Copy the property using the updated name.
                mol = mol.setProperty(new_prop, mol.property(prop)).molecule()

                # Store the amber types in the opposite end state.
                if prop[:-1] == "ambertype":
                    if lam == "0":
                        amber_types = mol.property("ambertype1").toVector()
                    else:
                        amber_types = mol.property("ambertype0").toVector()

                elif prop[:-1] == "element":
                    if lam == "0":
                        elements = mol.property("element1").toVector()
                    else:
                        elements = mol.property("element0").toVector()

                else:
                    # Delete redundant properties.
                    mol = mol.removeProperty(prop[:-1] + "1").molecule()
                    mol = mol.removeProperty(prop[:-1] + "0").molecule()

        # Convert ambertype and element property of dummies to those of the
        # other end state.
        if convert_amber_dummies:
            amber_type = property_map.get("ambertype", "ambertype")
            element = property_map.get("element", "element")
            if mol.hasProperty(amber_type) and mol.hasProperty(element):
                # Search for any dummy atoms.
                try:
                    search = mol.atoms("element Xx")
                except KeyError:
                    search = []

                # Replace the ambertype.
                for dummy in search:
                    index = dummy.index()
                    amber_type_value = amber_types[index.value()]
                    element_value = elements[index.value()]

                    # We have a dummy in both endstates so we try to infer the element.
                    # This is not general so it is only suitable for some common cases.
                    if element_value.symbol() == "Xx":
                        name = property_map.get("name", "name")
                        element_symbol = dummy.property(name)[0].upper()
                        element_value = _SireMol.Element(element_symbol)

                    mol = (
                        mol.atom(index)
                        .setProperty(amber_type, amber_type_value)
                        .molecule()
                    )
                    mol = mol.atom(index).setProperty(element, element_value).molecule()

                # Delete redundant properties.
                mol = mol.removeProperty("ambertype0").molecule()
                mol = mol.removeProperty("ambertype1").molecule()
                mol = mol.removeProperty("element0").molecule()
                mol = mol.removeProperty("element1").molecule()

        if generate_intrascale:
            # First we regenerate the connectivity based on the bonds.
            conn = _SireMol.Connectivity(mol.info()).edit()
            for bond in mol.property("bond").potentials():
                conn.connect(bond.atom0(), bond.atom1())
            mol.setProperty("connectivity", conn.commit())

            # Now we have the correct connectivity, we can regenerate the exclusions.
            gro_sys = _SireIO.GroTop(_System(mol)._sire_object).toSystem()
            mol.setProperty("intrascale", gro_sys[0].property("intrascale"))

        # Return the updated molecule.
        return Molecule(mol.commit())

    def _extractMolecule(self, property_map={}, is_lambda1=False):
        """
        Internal function to extract an "original" molecule from a merged
        molecule, i.e. one of the original molecules that was used to
        create the merge.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        is_lambda1 : bool
            Whether to use the molecule at the lambda = 1 end state.
            By default, the state at lambda = 0 is used.

        Returns
        -------

        molecule : BioSimSpace._SireWrappers.Molecule
            The molecule at the chosen end state, with dummy atoms removed.

        dummy_indices : [ int ]
            The indices of any dummy atoms in the original molecule.
        """

        if not isinstance(is_lambda1, bool):
            raise TypeError("'is_lambda1' must be of type 'bool'")

        if not self._is_perturbable:
            return Molecule(self._sire_object)

        # First extract the required end state.
        mol = self._toRegularMolecule(property_map=property_map, is_lambda1=is_lambda1)

        # Now find any non-dummy atoms in the molecule.
        query = "not element Xx"
        try:
            search_result = mol.search(query, property_map)
        except:
            msg = "All atoms in the selection are dummies. Unable to extract."
            if _isVerbose():
                raise _IncompatibleError(msg) from e
            else:
                raise _IncompatibleError(msg) from None

        # If there are no dummies, then simply return this molecule.
        if len(search_result) == mol.nAtoms():
            return mol, []

        else:
            # Store the indices of the non-dummy atoms.
            non_dummies = []
            for atom in search_result:
                non_dummies.append(atom.index())

            # Now search for the dummy atoms.
            query = "element Xx"
            try:
                search_result = mol.search(query, property_map)
            except:
                search_result = []

            # Store the indices of the dummy atoms.
            dummies = []
            for atom in search_result:
                dummies.append(atom.index())

            # Extract the non-dummy atoms from the molecule and return, along
            # with the indices of the dummy atoms.
            return mol.extract(non_dummies), dummies

    def _getPerturbationIndices(self):
        """
        Return the indices of the atoms that are perturbed, i.e. those
        that change one of the following properties: "ambertype", "LJ",
        or "charge".

        Returns
        -------

        idxs : [int]
            The indices of the atoms that are perturbed.
        """

        idxs = []

        if not self._is_perturbable:
            _warn(
                "You are trying to get the perturbation indices for a "
                "molecule that isn't perturbable!"
            )
            return idxs

        for idx, atom in enumerate(self.getAtoms()):
            atom = atom._sire_object
            if (
                atom.property("ambertype0") != atom.property("ambertype1")
                or atom.property("LJ0") != atom.property("LJ1")
                or atom.property("charge0") != atom.property("charge1")
            ):
                idxs.append(idx)

        return idxs

    def getCOMIdx(self):
        """Get the index of the atom that closest to the center of mass."""
        if self.isPerturbable():
            property_map = {"coordinates": "coordinates0", "mass": "mass0"}
        else:
            property_map = {"coordinates": "coordinates", "mass": "mass"}
        coords = self.coordinates(property_map=property_map)
        com = self._getCenterOfMass(property_map=property_map)
        com = _BSSVector(*[e / _angstrom for e in com])
        diffs = [coord.toVector() - com for coord in coords]
        sq_distances = [diff.dot(diff) for diff in diffs]
        idx = int(_np.argmin(sq_distances))
        return idx


# Import at bottom of module to avoid circular dependency.
from ._atom import Atom as _Atom
from ._molecules import Molecules as _Molecules
from ._residue import Residue as _Residue
from ._search_result import SearchResult as _SearchResult
from ._system import System as _System
