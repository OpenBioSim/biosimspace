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
A thin wrapper around Sire.System.System. This is an internal package and should
not be directly exposed to the user.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["System"]

import warnings as _warnings

from sire.legacy import IO as _SireIO
from sire.legacy import Maths as _SireMaths
from sire.legacy import Mol as _SireMol
from sire.legacy import System as _SireSystem
from sire.legacy import Vol as _SireVol
from sire.legacy import Units as _SireUnits

from .. import _isVerbose
from .._Exceptions import IncompatibleError as _IncompatibleError
from ..Types import Angle as _Angle
from ..Types import Length as _Length
from .. import Units as _Units

from ._sire_wrapper import SireWrapper as _SireWrapper

from sire.mol import Select as _Select


class System(_SireWrapper):
    """A container class for storing molecular systems."""

    def __init__(self, system):
        """
        Constructor.

        Parameters
        ----------

        system : Sire.System.System, :class:`System <BioSimSpace._SireWrappers.System>`, \
                 Sire.Mol._Mol.Molecule, :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                 :class:`Molecules <BioSimSpace._SireWrappers.Molecules>` \
                 [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
            A Sire or BioSimSpace System object, a Sire or BioSimSpace Molecule object,
            a BioSimSpace Molecules object, or a list of BioSimSpace molecule objects.
        """

        # Check that the system is valid.

        # Convert tuple to a list.
        if isinstance(system, tuple):
            system = list(system)

        # A Sire System object.
        if isinstance(system, _SireSystem.System):
            super().__init__(system)

        # Another BioSimSpace System object.
        elif isinstance(system, System):
            super().__init__(system._sire_object)

        # A Sire Molecule object.
        elif isinstance(system, _SireMol.Molecule):
            sire_object = _SireSystem.System("BioSimSpace_System.")
            super().__init__(sire_object)
            self.addMolecules(_Molecule(system))

        # A BioSimSpace Molecule object.
        elif isinstance(system, _Molecule):
            sire_object = _SireSystem.System("BioSimSpace_System.")
            super().__init__(sire_object)
            self.addMolecules(system)

        # A BioSimSpace Molecules object.
        elif isinstance(system, _Molecules):
            sire_object = _SireSystem.System("BioSimSpace_System.")
            super().__init__(sire_object)
            self.addMolecules(system)

        # A list of BioSimSpace Molecule objects.
        elif isinstance(system, list):
            if not all(isinstance(x, _Molecule) for x in system):
                raise TypeError(
                    "'system' must contain a list of 'BioSimSpace._SireWrappers.Molecule' types."
                )
            else:
                sire_object = _SireSystem.System("BioSimSpace_System.")
                super().__init__(sire_object)
                self.addMolecules(system)

        # Invalid type.
        else:
            raise TypeError(
                "'system' must be of type 'Sire.System.System', 'BioSimSpace._SireWrappers.System', "
                " Sire.Mol.Molecule', 'BioSimSpace._SireWrappers.Molecule', "
                "or a list of 'BioSimSpace._SireWrappers.Molecule' types. "
                f"The type {type(system)} is not supported."
            )

        # Flag that this object holds multiple atoms.
        self._is_multi_atom = True

        # Initialise dictionaries mapping MolNums to the number cumulative
        # number of atoms/residues in the system, up to that MolNum.
        self._atom_index_tally = {}
        self._residue_index_tally = {}

        # Initialise dictionary to map MolNum to MolIdx.
        self._molecule_index = {}

        # Store the molecule numbers.
        self._mol_nums = self._sire_object.molNums()

        # Initialise the iterator counter.
        self._iter_count = 0

        # Copy any fileformat property to each molecule.
        if "fileformat" in self._sire_object.propertyKeys():
            fileformat = self._sire_object.property("fileformat")
            for num in self._mol_nums:
                edit_mol = self._sire_object[num].edit()
                edit_mol = edit_mol.setProperty("fileformat", fileformat)
                self._sire_object.update(edit_mol.commit())
        else:
            # If a molecule has a fileformat property, use the first
            # that we find.
            for mol in self:
                if mol._sire_object.hasProperty("fileformat"):
                    self._sire_object.setProperty(
                        "fileformat", mol._sire_object.property("fileformat")
                    )
                    break

        # Reset the index mappings.
        self._reset_mappings()

        # Update the molecule numbers.
        self._mol_nums = self._sire_object.molNums()

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.System: nMolecules=%d>" % self.nMolecules()

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.System: nMolecules=%d>" % self.nMolecules()

    def __add__(self, other):
        """Addition operator."""

        # Create a copy of the current system.
        system = System(self._sire_object.__deepcopy__())

        # Add the new molecules.
        system.addMolecules(other)

        # Return the combined system.
        return system

    def __sub__(self, other):
        """Subtraction operator."""

        # Create a copy of the current system.
        system = System(self._sire_object.__deepcopy__())

        # Remove the molecules from the other system.
        if isinstance(other, System):
            system.removeMolecules(other.getMolecules())
        else:
            system.removeMolecules(other)

        # Return the new system.
        return system

    def __contains__(self, other):
        """Return whether other is in self."""

        if not isinstance(other, (_Molecule, _Atom, _Residue)):
            raise TypeError(
                "'other' must be of type 'BioSimSpace._SireWrappers.Molecule' "
                "'BioSimSpace._SireWrappers.Atom', or 'BioSimSpace._SireWrappers.Residue'."
            )

        # Return whether the object is in the system.
        return self._sire_object.contains(other._sire_object.molecule())

    def __getitem__(self, key):
        """Get a molecule from the system."""

        # Slice.
        if isinstance(key, slice):
            # Create a list to hold the molecules.
            molecules = []

            # Iterate over the slice.
            for x in range(*key.indices(self.nMolecules())):
                molecules.append(self[x])

            # Convert to a Molecules container and return.
            return _Molecules(molecules)

        # Index.
        else:
            try:
                key = int(key)
            except:
                raise TypeError("'key' must be of type 'int'")

            if key < -self.nMolecules() or key > self.nMolecules() - 1:
                raise IndexError("Molecules index is out of range.")

            if key < 0:
                key = key + self.nMolecules()

            # Extract and return the corresponding molecule.
            return _Molecule(self._sire_object[_SireMol.MolIdx(key)])

    def __setitem__(self, key, value):
        """Set a molecule in the container."""
        raise TypeError("'System' object does not support assignment.")

    def __iter__(self):
        """An iterator for the object."""
        # Reset the iterator counter and return the object.
        self._iter_count = 0
        return self

    def __next__(self):
        """An iterator for the object."""

        # Stop if we've reached the end of the container.
        if self._iter_count == self.nMolecules():
            raise StopIteration

        # Extract the next molecule in the container.
        molecule = self[self._iter_count]

        # Update the iterator counter.
        self._iter_count += 1

        # Return the molecule.
        return molecule

    def __len__(self):
        """Return the number of molecules in the system."""
        return self.nMolecules()

    def nMolecules(self):
        """
        Return the number of molecules in the system.

        Returns
        -------

        num_molecules : int
            The number of molecules in the system.
        """
        return self._sire_object.nMolecules()

    def nResidues(self):
        """
        Return the number of residues in the system.

        Returns
        -------

        num_residues : int
            The number of residues in the system.
        """

        tally = 0

        for n in self._sire_object.molNums():
            tally += self._sire_object[n].nResidues()

        return tally

    def nChains(self):
        """
        Return the number of chains in the system.

        Returns
        -------

        num_chains : int
            The number of chains in the system.
        """

        tally = 0

        for n in self._sire_object.molNums():
            tally += self._sire_object[n].nChains()

        return tally

    def nAtoms(self):
        """
        Return the number of atoms in the system.

        Returns
        -------

        num_atoms : int
            The number of atoms in the system.
        """

        tally = 0

        for n in self._sire_object.molNums():
            tally += self._sire_object[n].nAtoms()

        return tally

    def charge(self, property_map={}, is_lambda1=False):
        """
        Return the total molecular charge.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        is_lambda1 : bool
           Whether to use the charge at lambda = 1 if the molecule is merged.

        Returns
        -------

        charge : :class:`Charge <BioSimSpace.Types.Charge>`
            The molecular charge.
        """

        # Zero the charge.
        charge = 0 * _Units.Charge.electron_charge

        # Loop over all molecules and add the charge.
        for mol in self.getMolecules():
            charge += mol.charge(property_map, is_lambda1)

        # Return the total charge.
        return charge

    def fileFormat(self, property_map={}):
        """
        Return the file formats associated with the system.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        format : str
           The file formats associated with the system.
        """

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        prop = property_map.get("fileformat", "fileformat")

        try:
            return self._sire_object.property(prop).value()
        except:
            return None

    def isSame(
        self,
        other,
        excluded_properties=[],
        property_map0={},
        property_map1={},
        skip_water=True,
    ):
        """
        Check whether "other" is the same as this system.

        Parameters
        ----------

        other : :class:`System <BioSimSpace._SireWrappers.System>`
            Another BioSimSpace system.

        excluded_properties : [str]
            A list of properties to exclude when comparing systems. This allows
            you to check whether a subset of the system is the same, e.g.
            checking whether the topology is the same, but allowing different
            coordinates.

        property_map0 : dict
            A dictionary that maps "properties" in this system to their user
            defined values. This allows the user to refer to properties with
            their own naming scheme, e.g. { "charge" : "my-charge" }

        property_map1 : dict
            A dictionary that maps "properties" in "other" to their user
            defined values.

        skip_water : bool
            Whether to skip water molecules when comparing systems.

        Returns
        -------

        is_same : bool
            Whether the systems are the same.
        """

        # Validate input.

        if not isinstance(other, System):
            raise TypeError(
                "'other' must be of type 'BioSimSpace._SireWrappers.System'."
            )

        if not isinstance(excluded_properties, (list, tuple)):
            raise TypeError("'excluded_properties' must be a list of 'str' types.")

        if not all(isinstance(x, str) for x in excluded_properties):
            raise TypeError("'excluded_properties' must be a list of 'str' types.")

        if not isinstance(property_map0, dict):
            raise TypeError("'property_map0' must be of type 'dict'.")

        if not isinstance(property_map1, dict):
            raise TypeError("'property_map1' must be of type 'dict'.")

        if not isinstance(skip_water, bool):
            raise TypeError("'skip_water' must be of type 'bool'.")

        # If the system UIDs differ, then they are definitely different.
        if self._sire_object.uid() != other._sire_object.uid():
            return False

        # Return False if the systems have a different number of molecules,
        # atoms, or residues.
        if (
            self.nMolecules() != other.nMolecules()
            or self.nResidues() != other.nResidues()
            or self.nAtoms() != other.nAtoms()
        ):
            return False

        # Invert the property maps.
        inv_prop_map0 = {v: k for k, v in property_map0.items()}
        inv_prop_map1 = {v: k for k, v in property_map1.items()}

        # Add some additional properties to the excluded list. These are
        # used for internal metadata to aid object recovery.
        _excluded_properties = excluded_properties.copy()
        _excluded_properties.extend(["fileformat", "is_perturbable", "was_perturbable"])

        def _object_compare(object0, object1):
            """Helper function to check whether two Sire objects are the same."""

            # Store the two sets of properties.
            props0 = object0.propertyKeys()
            props1 = object1.propertyKeys()

            # Loop over all properties of object0.
            for p0 in props0:
                # Get the actual property name.
                name0 = inv_prop_map0.get(p0, p0)

                # Skip if excluded.
                if not name0 in _excluded_properties:
                    # Get the property name in other.
                    name1 = property_map1.get(name0, name0)

                    # Does object1 have this property?
                    if name1 in object1.propertyKeys():
                        # Do the property versions match?
                        try:
                            if object0.version(name0) != object1.version(name1):
                                return False
                        except:
                            # Get the properties in objects.
                            prop0 = object0.property(name0)
                            prop1 = object1.property(name1)
                            # Are the values the same?
                            try:
                                if prop0.value() != prop1.value():
                                    return False
                            except:
                                # Are they equal?
                                try:
                                    if prop0 != prop1:
                                        return False
                                except:
                                    return False

                    # Property is missing, so the objects differ.
                    else:
                        return False

            # Now check that there aren't any additional properties in object1
            # that aren't excluded.
            for p1 in props1:
                name1 = inv_prop_map1.get(p1, p1)
                # This is a property unique to object1, so they differ.
                if not name1 in _excluded_properties and name1 not in props0:
                    return False

            # If we get this far, then the objects are the same.
            return True

        # First compare system properties.
        is_same = _object_compare(self._sire_object, other._sire_object)

        if not is_same:
            return False

        if skip_water:
            molecules0 = self.search("not water", property_map0).molecules()
            molecules1 = other.search("not water", property_map1).molecules()
        else:
            molecules0 = self.getMolecules()
            molecules1 = other.getMolecules()

        # Now compare molecules.
        for molecule0, molecule1 in zip(molecules0, molecules1):
            if not _object_compare(molecule0._sire_object, molecule1._sire_object):
                return False

        # If we made it this far, then the systems are the same.
        return True

    def addMolecules(self, molecules):
        """
        Add a molecule, or list of molecules to the system.

        Parameters
        ----------

        molecules : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                    :class:`Molecules <BioSimSpace._SireWrappers.Molecules>`, \
                    [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`], \
                    :class:`System <BioSimSpace._SireWrappers.System>`
            A Molecule, Molecules object, a list of Molecule objects, or a System containing molecules.
        """

        # Whether the molecules are in a Sire container.
        is_sire_container = False

        # Convert tuple to a list.
        if isinstance(molecules, tuple):
            molecules = list(molecules)

        # A Molecule object.
        if isinstance(molecules, _Molecule):
            molecules = _Molecules([molecules])

        # A Molecules object.
        if isinstance(molecules, _Molecules):
            pass

        # A System object.
        elif isinstance(molecules, System):
            molecules = molecules.getMolecules()

        # A list of Molecule objects.
        elif isinstance(molecules, list) and all(
            isinstance(x, _Molecule) for x in molecules
        ):
            molecules = _Molecules(molecules)

        # Invalid argument.
        else:
            raise TypeError(
                "'molecules' must be of type 'BioSimSpace._SireWrappers.Molecule', "
                ", 'BioSimSpace._SireWrappers.System', or a list of "
                "'BioSimSpace._SireWrappers.Molecule' types."
            )

        # Store the existing number of molecules.
        num_mols = self._sire_object.nMolecules()

        # The system is empty: create an empty Sire system.
        if num_mols == 0:
            self._sire_object = self._createSireSystem()
            self._mol_nums = []

        # Only continue if there are molecules to add.
        if len(molecules) > 0:
            # Extract the molecule numbers for the current system and
            # the molecules to add.
            mol_nums0 = self._mol_nums
            mol_nums1 = molecules._sire_object.molNums()

            # There are molecule numbers in both sets, or the molecules
            # to add contains duplicates.
            if (set(mol_nums0) & set(mol_nums1)) or (
                len(mol_nums1) != len(set(mol_nums1))
            ):
                raise ValueError(
                    "'BioSimSpace._SireWrappers.System' can only "
                    "contain unique molecules. Use the 'copy' method "
                    "of 'BioSimSpace._SireWrappers' objects to "
                    "create a new version of them."
                )

            # Add the molecules to the system.
            self._sire_object.add(molecules._sire_object, _SireMol.MGName("all"))

            # Reset the index mappings.
            self._reset_mappings()

            # Update the molecule numbers.
            self._mol_nums = self._sire_object.molNums()

        # Remove velocities if any molecules are missing them.
        if self.nMolecules() > 1:
            # Search for water molecules in the system.
            try:
                mols_with_velocities = self.search(
                    f"mols with property velocity"
                ).molecules()
                num_vels = len(mols_with_velocities)
            except:
                num_vels = 0

            # Not all molecules have velocities.
            if num_vels > 0 and num_vels != self.nMolecules():
                _warnings.warn(
                    "Not all molecules have velocities. The 'velocity' property will be removed."
                )
                try:
                    self._sire_object = _SireIO.removeProperty(
                        self._sire_object, "velocity"
                    )
                except:
                    _warnings.warn(
                        "Failed to remove 'velocity' property from all molecules!"
                    )

    def removeMolecules(self, molecules):
        """
        Remove a molecule, or list of molecules from the system.

        Parameters
        ----------

        molecules : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                    :class:`Molecules <BioSimSpace._SireWrappers.Molecules>`, \
                    [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
            A Molecule, Molecules object, or list of Molecule objects.
        """

        # Whether the molecules are in a Sire container.
        is_sire_container = False

        # Convert tuple to a list.
        if isinstance(molecules, tuple):
            molecules = list(molecules)

        # A Molecule object.
        if isinstance(molecules, _Molecule):
            molecules = [molecules]

        # A Molecules object.
        if isinstance(molecules, _Molecules):
            is_sire_container = True

        # A list of Molecule objects.
        elif isinstance(molecules, list) and all(
            isinstance(x, _Molecule) for x in molecules
        ):
            pass

        # Invalid argument.
        else:
            raise TypeError(
                "'molecules' must be of type 'BioSimSpace._SireWrappers.Molecule' "
                "or a list of 'BioSimSpace._SireWrappers.Molecule' types."
            )

        # Remove the molecules in the system.
        if is_sire_container:
            self._sire_object.remove(molecules._sire_object, _SireMol.MGName("all"))
        else:
            for mol in molecules:
                self._sire_object.remove(mol._sire_object.number())

        # Reset the index mappings.
        self._reset_mappings()

        # Update the molecule numbers.
        self._mol_nums = self._sire_object.molNums()

    def removeWaterMolecules(self):
        """Remove all of the water molecules from the system."""

        # Get the list of water molecules.
        waters = self.getWaterMolecules()

        # Remove the molecules in the system.
        self._sire_object.remove(waters._sire_object, _SireMol.MGName("all"))

        # Reset the index mappings.
        self._reset_mappings()

        # Update the molecule numbers.
        self._mol_nums = self._sire_object.molNums()

    def updateMolecule(self, index, molecule):
        """
        Update the molecule at the given index.

        Parameters
        ----------

        index : int
            The index of the molecule.

        molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            The updated (or replacement) molecule.
        """

        if type(index) is not int:
            raise TypeError("'index' must be of type 'int'")

        if index < -self.nMolecules() or index >= self.nMolecules():
            raise IndexError("The molecule 'index' is out of range.")

        if not isinstance(molecule, _Molecule):
            raise TypeError(
                "'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'."
            )

        # The molecule numbers match.
        if self[index].number() == molecule.number():
            self.updateMolecules(molecule)

        # The molecule numbers don't match.
        else:
            # Create a copy of the system.
            system = self.copy()._sire_object

            # Update the molecule in the system, preserving the
            # original molecular ordering.
            system = _SireIO.updateAndPreserveOrder(
                system, molecule._sire_object, index
            )

            # Update the Sire object.
            self._sire_object = system

            # Reset the index mappings.
            self._reset_mappings()

            # Update the molecule numbers.
            self._mol_nums = self._sire_object.molNums()

    def updateMolecules(self, molecules):
        """
        Update a molecule, or list of molecules in the system.

        Parameters
        ----------

        molecules : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                    [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
            A Molecule, or list of Molecule objects.
        """

        # Convert tuple to a list.
        if isinstance(molecules, tuple):
            molecules = list(molecules)

        # A Molecule object.
        if isinstance(molecules, _Molecule):
            molecules = [molecules]

        # A list of Molecule objects.
        elif isinstance(molecules, list) and all(
            isinstance(x, _Molecule) for x in molecules
        ):
            pass

        # Invalid argument.
        else:
            raise TypeError(
                "'molecules' must be of type 'BioSimSpace._SireWrappers.Molecule' "
                "or a list of 'BioSimSpace._SireWrappers.Molecule' types."
            )

        # Update each of the molecules.
        # TODO: Currently the Sire.System.update method doesn't work correctly
        # for certain changes to the Molecule molInfo object. As such, we remove
        # the old molecule from the system, then add the new one in. Make sure to
        # operate on a copy of the system since the original will be destroyed if
        # an exception is thrown.
        for mol in molecules:
            # Only try to update the molecule if it exists in the system.
            if _SireMol.MolNum(mol.number()) in self._mol_nums:
                try:
                    system = self.copy()._sire_object
                    system.update(mol._sire_object)
                except:
                    # The UUID of the molecule has changed. Normally we would
                    # need to remove and re-add the molecule, but this would
                    # change the molecular ordering. Instead we create a new
                    # system by re-adding all molecules in the original order.

                    system = self.copy()._sire_object

                    # Get the index of the molecule in the existing molNums
                    # list.
                    idx = self._mol_nums.index(_SireMol.MolNum(mol.number()))

                    # Update the molecule in the system, preserving the
                    # original molecular ordering.
                    system = _SireIO.updateAndPreserveOrder(
                        system, mol._sire_object, idx
                    )
            else:
                raise ValueError(f"System doesn't contain molecule: {mol}")

        # Update the Sire object.
        self._sire_object = system

        # Reset the index mappings.
        self._reset_mappings()

        # Update the molecule numbers.
        self._mol_nums = self._sire_object.molNums()

    def getMolecule(self, index):
        """
        Return the molecule at the given index.

        Parameters
        ----------

        index : int
            The index of the molecule.

        Returns
        -------

        molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            The requested molecule.
        """
        return self[index]

    def getMolecules(self, group="all"):
        """
        Return a list containing all of the molecules in the specified group.

        Parameters
        ----------

        group : str
            The name of the molecule group.

        Returns
        -------

        molecules : [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
            The list of molecules in the group.
        """

        if not isinstance(group, str):
            raise TypeError("'group' must be of type 'str'")

        # Try to extract the molecule group.
        try:
            molgrp = self._sire_object.group(_SireMol.MGName(group))
        except:
            raise ValueError("No molecules in group '%s'" % group)

        # Return a molecules container.
        return _Molecules(molgrp)

    def getAtoms(self):
        """
        Return a list containing the atoms in the system.

        Returns
        -------

        atoms : [:class:`Atom <BioSimSpace._SireWrappers.Atom>`]
            The list of atoms in the system.
        """
        atoms = []
        for mol in self:
            atoms.extend(mol.getAtoms())
        return atoms

    def getAtom(self, index):
        """
        Return the atom specified by the absolute index.

        Parameters
        ----------

        index : int
            The absolute index of the atom in the system.

        Returns
        -------

        atom : :class:`Atom <BioSimSpace._SireWrappers.Atom>`
            The atom at the specified index.
        """
        if not type(index) is int:
            raise TypeError("'index' must be of type 'int'.")

        # Set the MolNum to atom index mapping.
        self._set_atom_index_tally()

        # Work out the total number of atoms.
        total_atoms = list(self._atom_index_tally.values())[-1] + self[-1].nAtoms()

        if index < 0 or index >= total_atoms:
            raise ValueError(f"'index' must be in range (0-{total_atoms}].")

        # Loop backwards over the molecules until we find the one
        # containing the index.
        for mol_idx in range(self.nMolecules() - 1, -1, -1):
            mol_num = self._mol_nums[mol_idx]
            if index >= self._atom_index_tally[mol_num]:
                break

        # Get the relative index of the atom in the molecule.
        rel_idx = index - self._atom_index_tally[mol_num]

        # Return the atom.
        return self[mol_idx].getAtoms()[rel_idx]

    def getResidues(self):
        """
        Return a list containing the residues in the system.

        Returns
        -------

        residues : [:class:`Residue <BioSimSpace._SireWrappers.Residue>`]
            The list of residues in the system.
        """
        residues = []
        for mol in self:
            residues.extend(mol.getResidues())
        return residues

    def getResidue(self, index):
        """
        Return the residue specified by the absolute index.

        Parameters
        ----------

        index : int
            The absolute index of the residue in the system.

        Returns
        -------

        residue : :class:`Residue <BioSimSpace._SireWrappers.Residue>`
            The residue at the specified index.
        """
        if not type(index) is int:
            raise TypeError("'index' must be of type 'int'.")

        # Set the MolNum to residue index mapping.
        self._set_residue_index_tally()

        # Work out the total number of residues.
        total_residues = (
            list(self._residue_index_tally.values())[-1] + self[-1].nResidues()
        )

        if index < 0 or index >= total_residues:
            raise ValueError(f"'index' must be in range (0-{total_residues}].")

        # Loop backwards over the molecules until we find the one
        # containing the index.
        for mol_idx in range(self.nMolecules() - 1, -1, -1):
            mol_num = self._mol_nums[mol_idx]
            if index >= self._residue_index_tally[mol_num]:
                break

        # Get the relative index of the residue in the molecule.
        rel_idx = index - self._residue_index_tally[mol_num]

        # Return the atom.
        return self[mol_idx].getResidues()[rel_idx]

    def getWaterMolecules(self, property_map={}):
        """
        Return a list containing all of the water molecules in the system.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        molecules : :class:`Molecules <BioSimSpace._SireWrappers.Molecules>`
            A container of water molecule objects. The container will be
            empty if no water molecules are present.
        """
        return _Molecules(self._sire_object.search("water").toGroup())

    def nWaterMolecules(self):
        """
        Return the number of water molecules in the system.

        Returns
        -------

        num_waters : int
            The number of water molecules in the system.
        """
        return len(self.getWaterMolecules())

    def getPerturbableMolecules(self):
        """
        Return a list containing all of the perturbable molecules in the system.

        Returns
        -------

        molecules : :class:`Molecules <BioSimSpace._SireWrappers.Molecules>`
            A container of perturbable molecule objects. The container will
            be empty if no perturbable molecules are present.
        """
        return _Molecules(
            self._sire_object.search("molecules with property is_perturbable").toGroup()
        )

    def nPerturbableMolecules(self):
        """
        Return the number of perturbable molecules in the system.

        Returns
        -------

        num_perturbable : int
            The number of perturbable molecules in the system.
        """
        return len(self.getPerturbableMolecules())

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

        # Repartion hydrogen masses for all molecules in this system.
        self._sire_object = _SireIO.repartitionHydrogenMass(
            self._sire_object, factor, water_options[water], pmap
        )

    def search(self, query, property_map={}):
        """
        Search the system for atoms, residues, and molecules. Search results
        will be reduced to their minimal representation, i.e. a molecule
        containing a single residue will be returned as a residue.

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
                   :class:`Residue <BioSimSpace._SireWrappers.Residue>`, \
                   :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, ...]
            A list of objects matching the search query.

        Examples
        --------

        Search for residues names ALA.

        >>> result = system.search("resname ALA")

        Search for residues names ALA using a case insensitive search.

        >>> result = system.search("resname /ala/i")

        Search for all carbon atoms in residues named ALA.

        >>> result = system.search("resname ALA and atomname C")

        Search for all oxygen or hydrogen atoms.

        >>> result = system.search("element oxygen or element hydrogen")

        Search for atom index 23 in molecule index 10.

        >>> result = system.search("molidx 10 and atomidx 23")
        """

        if not isinstance(query, str):
            raise TypeError("'query' must be of type 'str'")

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        # Initialise a list to hold the search results.
        results = []

        try:
            # Query the Sire system.
            search_result = _Select(query)(self._sire_object, property_map)

        except Exception as e:
            msg = "'Invalid search query: %r : %s" % (query, e)
            if _isVerbose():
                raise ValueError(msg) from e
            else:
                raise ValueError(msg) from None

        return _SearchResult(search_result)

    def getIndex(self, item):
        """
        Convert indices of atoms and residues to their absolute values in
        this system.

        Parameters
        ----------

        item : :class:`Atom <BioSimSpace._SireWrappers.Atom>`, \
               :class:`Residue <BioSimSpace._SireWrappers.Residue>` \
               :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            An Atom, Residue, or Molecule object from the System, or a
            list containing objects of these types.

        Returns
        -------

        index : int, [int]
            The absolute index of the atom/residue/molecule in the system.
        """

        # Convert single object to list.
        if not isinstance(item, (tuple, list)):
            item = [item]

        # Create a list to hold the indices.
        indices = []

        # Loop over all of the items.
        for x in item:
            if x is None:
                continue
            # Atom.
            if isinstance(x, _Atom):
                # Create the MolNum to atom index mapping dictionary.
                self._set_atom_index_tally()

                # Get the AtomIdx and MolNum from the Atom.
                index = x.index()
                mol_num = x._sire_object.molecule().number()

                try:
                    index += self._atom_index_tally[mol_num]
                except KeyError:
                    raise KeyError(
                        "The atom belongs to molecule '%s' that is not part of "
                        "this system!" % mol_num
                    )

                indices.append(index)

            # Residue.
            elif isinstance(x, _Residue):
                # Create the MolNum to residue index mapping dictionary.
                self._set_residue_index_tally()

                # Get the ResIdx and MolNum from the Residue.
                index = x.index()
                mol_num = x._sire_object.molecule().number()

                try:
                    index += self._residue_index_tally[mol_num]
                except KeyError:
                    raise KeyError(
                        "The residue belongs to molecule '%s' that is not part of "
                        "this system!" % mol_num
                    )

                indices.append(index)

            # Molecule.
            elif isinstance(x, _Molecule):
                # Create the MolNum to molecule index mapping dictionary.
                self._set_molecule_index_tally()

                # Get the MolNum from the molecule.
                mol_num = x._sire_object.molecule().number()

                try:
                    index = self._molecule_index[mol_num]
                except KeyError:
                    raise KeyError(
                        "The molecule '%s' is not part of this system!" % mol_num
                    )

                indices.append(index)

            # Unsupported.
            else:
                raise TypeError(
                    f"'item' (class {type(x)}) must be of type 'BioSimSpace._SireWrappers.Atom' "
                    "or 'BioSimSpace._SireWrappers.Residue'"
                )

        # If this was a single object, then return a single index.
        if len(indices) == 1:
            return indices[0]
        # Return the list of indices.
        else:
            return indices

    def setBox(self, box, angles=3 * [_Angle(90, "degree")], property_map={}):
        """
        Set the size of the periodic simulation box.

        Parameters
        ----------

        box : [:class:`Length <BioSimSpace.Types.Length>`]
            The box vector magnitudes in each dimension.

        angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
            The box vector angles: yz, xz, and xy.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Convert tuples to lists.
        if isinstance(box, tuple):
            box = list(box)
        if isinstance(angles, tuple):
            angles = list(angles)

        # Validate input.
        if not isinstance(box, list) or not all(isinstance(x, _Length) for x in box):
            raise TypeError(
                "'box' must be a list of 'BioSimSpace.Types.Length' objects."
            )

        if not isinstance(angles, list) or not all(
            isinstance(x, _Angle) for x in angles
        ):
            raise TypeError(
                "'angles' must be a list of 'BioSimSpace.Types.Angle' objects."
            )

        if len(box) != 3:
            raise ValueError("'box' must contain three items.")

        if len(box) != 3:
            raise ValueError("'angles' must contain three items.")

        # Convert sizes to Anstrom.
        vec = [x.angstroms().value() for x in box]

        # Whether the box is triclinic.
        is_triclinic = False

        # If any angle isn't 90 degrees, then the box is triclinic.
        for x in angles:
            if abs(x.degrees().value() - 90) > 1e-6:
                is_triclinic = True
                break

        if is_triclinic:
            # Convert angles to Sire units.
            ang = [x.degrees().value() * _SireUnits.degree for x in angles]

            # Create a triclinic box object.
            space = _SireVol.TriclinicBox(
                vec[0], vec[1], vec[2], ang[0], ang[1], ang[2]
            )

        else:
            # Create a periodic box object.
            space = _SireVol.PeriodicBox(_SireMaths.Vector(vec))

        # Set the "space" property.
        self._sire_object.setProperty(property_map.get("space", "space"), space)

    def getBox(self, property_map={}):
        """
        Get the size of the periodic simulation box.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        box_size : [:class:`Length <BioSimSpace.Types.Length>`]
            The box vector magnitudes in each dimension.

        angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
            The box vector angles: yz, xz, and xy.
        """

        # Get the "space" property and convert to a list of BioSimSpace.Type.Length
        # objects.
        try:
            space = self._sire_object.property(property_map.get("space", "space"))

            # Periodic box.
            if isinstance(space, _SireVol.PeriodicBox):
                box = [_Length(x, "Angstrom") for x in space.dimensions()]
                angles = 3 * [_Angle(90, "degrees")]

            # TriclinicBox box.
            elif isinstance(space, _SireVol.TriclinicBox):
                box = [
                    _Length(space.vector0().magnitude(), "Angstrom"),
                    _Length(space.vector1().magnitude(), "Angstrom"),
                    _Length(space.vector2().magnitude(), "Angstrom"),
                ]
                angles = [
                    _Angle(space.alpha(), "degree"),
                    _Angle(space.beta(), "degree"),
                    _Angle(space.gamma(), "degree"),
                ]

            else:
                raise TypeError(f"Unsupported box type: {space} - {type(space)}")
        except:
            box = None
            angles = None

        return box, angles

    def translate(self, vector, property_map={}):
        """
        Translate the system.

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

        # Translate each of the molecules in the system.
        for mol in self.getMolecules():
            mol.translate(vector, property_map)
            self._sire_object.update(mol._sire_object)

    def getRestraintAtoms(
        self,
        restraint,
        mol_index=None,
        is_absolute=True,
        allow_zero_matches=False,
        property_map={},
    ):
        """
        Get the indices of atoms involved in a restraint.

        Parameters
        ----------

        restraint : str
            The type of restraint.

        mol_index : int
            The index of the molecule of interest. If None, then the
            entire system is searched.

        is_absolute : bool
            Whether the indices are absolute, i.e. indexed within the
            entire system. If False, then indices are relative to the
            molecule in which the atom is found.

        allow_zero_matches : bool
            Whether to raise an exception when no atoms match the restraint.
            Setting to false is useful if you need to determine restraints
            on a per-molecule basis, i.e. so molecules won't contain atoms
            matching the restraint, but the entire system will.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their

        Returns
        -------

        indices : [int]
            A list of the backbone atom indices.
        """

        if not isinstance(restraint, str):
            raise TypeError("'restraint' must be of type 'str'.")

        # Allowed keyword options.
        allowed = ["backbone", "heavy", "all"]

        # Convert to lower case and strip whitespace.
        restraint = restraint.lower().replace(" ", "")
        if restraint not in allowed:
            raise ValueError(f"'restraint' must be one of: {allowed}")

        if mol_index is not None and not type(mol_index) is int:
            raise TypeError("'mol_index' must be of type 'int'.")

        if not isinstance(is_absolute, bool):
            raise TypeError("'is_absolute' must be of type 'bool'.")

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'.")

        # Initialise the list of indices.
        indices = []

        # A set of protein residues. Taken from MDAnalysis.
        prot_res = {
            # CHARMM top_all27_prot_lipid.rtf
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HSD",
            "HSE",
            "HSP",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
            "ALAD",
            ## 'CHO','EAM', # -- special formyl and ethanolamine termini of gramicidin
            # PDB
            "HIS",
            "MSE",
            # from Gromacs 4.5.3 oplsaa.ff/aminoacids.rtp
            "ARGN",
            "ASPH",
            "CYS2",
            "CYSH",
            "QLN",
            "PGLU",
            "GLUH",
            "HIS1",
            "HISD",
            "HISE",
            "HISH",
            "LYSH",
            # from Gromacs 4.5.3 gromos53a6.ff/aminoacids.rtp
            "ASN1",
            "CYS1",
            "HISA",
            "HISB",
            "HIS2",
            # from Gromacs 4.5.3 amber03.ff/aminoacids.rtp
            "HID",
            "HIE",
            "HIP",
            "ORN",
            "DAB",
            "LYN",
            "HYP",
            "CYM",
            "CYX",
            "ASH",
            "GLH",
            "ACE",
            "NME",
            # from Gromacs 2016.3 amber99sb-star-ildn.ff/aminoacids.rtp
            "NALA",
            "NGLY",
            "NSER",
            "NTHR",
            "NLEU",
            "NILE",
            "NVAL",
            "NASN",
            "NGLN",
            "NARG",
            "NHID",
            "NHIE",
            "NHIP",
            "NTRP",
            "NPHE",
            "NTYR",
            "NGLU",
            "NASP",
            "NLYS",
            "NPRO",
            "NCYS",
            "NCYX",
            "NMET",
            "CALA",
            "CGLY",
            "CSER",
            "CTHR",
            "CLEU",
            "CILE",
            "CVAL",
            "CASF",
            "CASN",
            "CGLN",
            "CARG",
            "CHID",
            "CHIE",
            "CHIP",
            "CTRP",
            "CPHE",
            "CTYR",
            "CGLU",
            "CASP",
            "CLYS",
            "CPRO",
            "CCYS",
            "CCYX",
            "CMET",
            "CME",
            "ASF",
        }

        # A list of ion elements.
        ions = [
            "F",
            "Cl",
            "Br",
            "I",
            "Li",
            "Na",
            "K",
            "Rb",
            "Cs",
            "Mg",
            "Tl",
            "Cu",
            "Ag",
            "Be",
            "Cu",
            "Ni",
            "Pt",
            "Zn",
            "Co",
            "Pd",
            "Ag",
            "Cr",
            "Fe",
            "Mg",
            "V",
            "Mn",
            "Hg",
            "Cd",
            "Yb",
            "Ca",
            "Sn",
            "Pb",
            "Eu",
            "Sr",
            "Sm",
            "Ba",
            "Ra",
            "Al",
            "Fe",
            "Cr",
            "In",
            "Tl",
            "Y",
            "La",
            "Ce",
            "Pr",
            "Nd",
            "Sm",
            "Eu",
            "Gd",
            "Tb",
            "Dy",
            "Er",
            "Tm",
            "Lu",
            "Hf",
            "Zr",
            "Ce",
            "U",
            "Pu",
            "Th",
        ]

        # Whether we've searched a perturbable system. If so, then we need to process
        # the search results differently. (There will be one for each molecule.)
        is_perturbable_system = False

        # Search the entire system.
        if mol_index is None:
            # Only search the system directly if there are no perturbable molecules.
            if self.nPerturbableMolecules() == 0:
                # Backbone restraints.
                if restraint == "backbone":
                    # Find all N, CA, C, and O atoms in protein residues.
                    string = (
                        "(not water) and (resname "
                        + ",".join(prot_res)
                        + ") and (atomname N,CA,C,O)"
                    )
                    try:
                        search = self.search(string, property_map)
                    except:
                        search = []

                elif restraint == "heavy":
                    # Convert to a formatted string for the search.
                    ion_string = ",".join(ions + ["H,Xx"])
                    # Find all non-water, non-hydrogen, non-ion elements.
                    string = f"(not water) and (not element {ion_string})"
                    try:
                        search = self.search(string, property_map)
                    except:
                        search = []

                elif restraint == "all":
                    # Convert to a formatted string for the search.
                    ion_string = ",".join(ions)
                    # Find all non-water, non-ion elements.
                    string = f"(not water) and (not element {ion_string})"
                    try:
                        search = self.search(string, property_map)
                    except:
                        search = []

            # Search each molecule individually, using the property map to specify the
            # correct name for the "element" property in any perturble molecules.
            else:
                is_perturbable_system = True

                # Initialise a list to hold all of the search results.
                results = []

                for mol in self.getMolecules():
                    # Reset the local property map.
                    _property_map = property_map.copy()

                    # Initialise an empty search.
                    search = []

                    # If perturbable, use the 'element0' property for the search
                    # unless it's already been set to the lambda = 1 property.
                    if mol.isPerturbable():
                        if _property_map.get("element") != "element1":
                            _property_map["element"] = "element0"

                    if restraint == "backbone":
                        if not mol.isWater():
                            # Find all N, CA, C, and O atoms in protein residues.
                            string = (
                                "(resname "
                                + ",".join(prot_res)
                                + ") and (atomname N,CA,C,O)"
                            )
                            try:
                                search = mol.search(string, _property_map)
                            except:
                                search = []

                    elif restraint == "heavy":
                        if not mol.isWater():
                            # Convert to a formatted string for the search.
                            ion_string = ",".join(ions + ["H,Xx"])
                            # Find all non-water, non-hydrogen, non-ion elements.
                            string = f"not element {ion_string}"
                            try:
                                search = mol.search(string, _property_map)
                            except:
                                search = []

                    elif restraint == "all":
                        if not mol.isWater():
                            # Convert to a formatted string for the search.
                            ion_string = ",".join(ions)
                            # Find all non-water, non-ion elements.
                            string = f"not element {ion_string}"
                            try:
                                search = mol.search(string, _property_map)
                            except:
                                search = []

                    # Append the search result for this molecule.
                    if len(search) > 0:
                        results.append(search)

        # Search the chosen molecule.
        else:
            # Create a local copy of the property map.
            _property_map = property_map.copy()

            # Initialise an empty search.
            search = []

            # Extract the molecule.
            mol = self[mol_index]

            # If perturbable, use the 'element0' property for the search
            # unless it's already been set to the lambda = 1 property.
            if mol.isPerturbable():
                if _property_map.get("element") != "element1":
                    _property_map["element"] = "element0"

            if restraint == "backbone":
                if not mol.isWater():
                    # Find all N, CA, C, and O atoms in protein residues.
                    string = (
                        "(resname " + ",".join(prot_res) + ") and (atomname N,CA,C,O)"
                    )
                    try:
                        search = mol.search(string, _property_map)
                    except:
                        search = []

            elif restraint == "heavy":
                if not mol.isWater():
                    # Convert to a formatted string for the search.
                    ion_string = ",".join(ions + ["H,Xx"])
                    # Find all non-water, non-hydrogen, non-ion elements.
                    string = f"not element {ion_string}"
                    try:
                        search = mol.search(string, _property_map)
                    except:
                        search = []

            elif restraint == "all":
                if not mol.isWater():
                    # Convert to a formatted string for the search.
                    ion_string = ",".join(ions)
                    # Find all non-water, non-ion elements.
                    string = f"not element {ion_string}"
                    try:
                        search = mol.search(string, _property_map)
                    except:
                        search = []

        if is_perturbable_system:
            # Raise an exception if no atoms match the restraint.
            if not allow_zero_matches:
                num_results = 0
                for search in results:
                    num_results += len(search)
                if num_results == 0:
                    msg = "No atoms matched the restraint!"
                    if restraint == "backbone":
                        msg += " Backbone restraints only apply to atoms in protein residues."
                    raise _IncompatibleError(msg)

            # Loop over the searches for each molecule.
            for search in results:
                # Now loop over all matching atoms and get their indices.
                for atom in search:
                    if is_absolute:
                        indices.append(self.getIndex(atom))
                    else:
                        indices.append(atom.index())

        else:
            # Raise an exception if no atoms match the restraint.
            if len(search) == 0 and not allow_zero_matches:
                msg = "No atoms matched the restraint!"
                if restraint == "backbone":
                    msg += (
                        " Backbone restraints only apply to atoms in protein residues."
                    )
                raise _IncompatibleError(msg)

            # Now loop over all matching atoms and get their indices.
            for atom in search:
                if is_absolute:
                    indices.append(self.getIndex(atom))
                else:
                    indices.append(atom.index())

        # The indices should be sorted, but do so just in case.
        indices.sort()

        return indices

    def _isParameterised(self, property_map={}):
        """
        Whether the system is parameterised, i.e. can we run a simulation
        using this system. Essentially we check whether every molecule in
        the system has "bond" and "LJ" properties and assume that is
        parameterised if so.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their

        Returns
        -------

        is_parameterised : bool
            Whether the system is parameterised.
        """

        # Get the "bond" property.
        bond = property_map.get("bond", "bond")
        # Test for perturbable molecules too.
        bond0 = bond + "0"

        # Get the "LJ" property.
        LJ = property_map.get("LJ", "LJ")
        # Test for perturbable molecules too.
        LJ0 = LJ + "0"

        # Check each molecule for "bond" and "LJ" properties.
        for mol in self.getMolecules():
            props = mol._sire_object.propertyKeys()
            if (bond not in props and bond0 not in props) or (
                LJ not in props and LJ0 not in props
            ):
                return False

        # If we get this far, then all molecules are okay.
        return True

    def _getAABox(self, property_map={}):
        """
        Get the axis-aligned bounding box for the molecular system.

        Parameters
        ----------

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        aabox : Sire.Vol.AABox
            The axis-aligned bounding box for the molecule.
        """

        # Initialise the coordinates vector.
        coord = []

        # Get all of the molecules in the system.
        mols = self.getMolecules()

        # Loop over all of the molecules.
        for idx, mol in enumerate(mols):
            # Extract the atomic coordinates and append them to the vector.
            try:
                if "coordinates" in property_map:
                    prop = property_map["coordinates"]
                else:
                    if mol.isPerturbable():
                        prop = "coordinates0"
                    else:
                        prop = "coordinates"
                coord.extend(mol._sire_object.property(prop).toVector())

            except UserWarning as e:
                msg = (
                    "Unable to compute the axis-aligned bounding "
                    + "box since a molecule has no 'coordinates' property."
                )
                if _isVerbose():
                    raise _IncompatibleError(msg) from e
                else:
                    raise _IncompatibleError(msg) from None

        # Return the AABox for the coordinates.
        return _SireVol.AABox(coord)

    def _renumberMolecules(self, molecules, is_rebuild=False):
        """
        Helper function to renumber the molecules to be consistent with the
        system.

        Parameters
        ----------

        molecules : [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
            A list of molecule objects.

        Returns
        -------

        molecules : [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
            The renumber list of molecule objects.
        """

        # Renumber everything.
        if is_rebuild:
            num_molecules = 0
            num_residues = 0
            num_atoms = 0

        # Get the current number of molecules, residues, and atoms.
        else:
            num_molecules = self.nMolecules()
            num_residues = self.nResidues()
            num_atoms = self.nAtoms()

        # Create a list to hold the modified molecules.
        new_molecules = []

        # Loop over all of the molecules.
        for mol in molecules:
            # Create a copy of the molecule.
            new_mol = _Molecule(mol)

            # Get the Sire molecule and make it editable.
            edit_mol = new_mol._sire_object.edit()

            # Renumber the molecule.
            edit_mol = edit_mol.renumber(_SireMol.MolNum(num_molecules + 1)).molecule()
            num_molecules += 1

            # A hash mapping between old and new numbers.
            num_hash = {}

            # Loop over all residues and add them to the hash.
            for res in edit_mol.residues():
                num_hash[res.number()] = _SireMol.ResNum(num_residues + 1)
                num_residues += 1

            # Renumber the residues.
            edit_mol = edit_mol.renumber(num_hash).molecule()

            # Clear the hash.
            num_hash = {}

            # Loop over all of the atoms and add them to the hash.
            for atom in edit_mol.atoms():
                num_hash[atom.number()] = _SireMol.AtomNum(num_atoms + 1)
                num_atoms += 1

            # Renumber the atoms.
            edit_mol = edit_mol.renumber(num_hash).molecule()

            # Commit the changes and replace the molecule.
            new_mol._sire_object = edit_mol.commit()

            # Append to the list of molecules.
            new_molecules.append(new_mol)

        # Return the renumbered molecules.
        return new_molecules

    def _createSireSystem(self):
        """
        Create an empty Sire system with a single molecule group called "all".

        Returns
        -------

        system : Sire.System.System
            A Sire system object.
        """

        # Create an empty Sire System.
        system = _SireSystem.System("BioSimSpace_System")

        # Create a new "all" molecule group.
        molgrp = _SireMol.MoleculeGroup("all")

        # Add the molecule group to the system.
        system.add(molgrp)

        # Copy any existing system properties.
        for prop in self._sire_object.propertyKeys():
            system.setProperty(prop, self._sire_object.property(prop))

        return system

    def _getRelativeIndices(self, abs_index):
        """
        Given an absolute index, get the number of the molecule to which it
        belongs and the relative index within that molecule.

        Parameters
        ----------

        abs_index : int
            The absolute index of the atom in the system.

        Returns
        -------

        mol_index : int
            The molecule index to which the atom belongs.

        rel_index : int
            The relative index of the atom in the molecule to which it
            belongs.
        """
        # Make sure the atom index tally has been created.
        self.getIndex(self[0].getAtoms()[0])

        # Set tally counters for the total number of atoms to date
        # and the total up to the previous molecule.
        tally = 0
        tally_last = 0

        # Set the current and previous MolNum.
        mol_num = self._mol_nums[0]
        mol_num_last = self._mol_nums[0]

        # Loop over each molecule until the total atom number is equal
        # to or greater than the absolute index. When it is, return the
        # absolute index minus the tally up to the previous molecule.
        for num, num_atoms in self._atom_index_tally.items():
            if abs_index >= tally:
                tally_last = tally
                tally = num_atoms
                mol_num_last = mol_num
                mol_num = num
            else:
                mol_idx = self._molecule_index[mol_num_last]
                return mol_idx, abs_index - tally_last

        raise ValueError("'abs_index' exceeded system atom tally!")

    def _reset_mappings(self):
        """Internal function to reset index mapping dictionaries."""

        # Clear dictionaries.
        self._molecule_index = {}
        self._atom_index_tally = {}
        self._residue_index_tally = {}

        # Rebuild the MolNum to index mapping.
        for idx in range(0, self.nMolecules()):
            self._molecule_index[self._sire_object[_SireMol.MolIdx(idx)].number()] = idx

    def _set_water_topology(self, format, property_map={}):
        """
        Internal function to swap the water topology to AMBER or GROMACS format.

        Parameters
        ----------

        format : string
            The format to convert to: either "AMBER" or "GROMACS".

        property_map : dict
           A dictionary that maps system "properties" to their user defined
           values. This allows the user to refer to properties with their
           own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Validate input.

        if not isinstance(format, str):
            raise TypeError("'format' must be of type 'str'")

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        # Strip whitespace and convert to upper case.
        format = format.replace(" ", "").upper()

        # We allow conversion to AMBER or GROMACS format water toplogies. While
        # AMBER requires a specific topology, GROMACS can handle whatever. We
        # leave this format as an option in case we want to allow the user to
        # swap the topology in future.
        if format not in ["AMBER", "GROMACS"]:
            raise ValueError("'format' must be 'AMBER' or 'GROMACS'.")

        # Get the water molecules.
        waters = self.getWaterMolecules()

        if len(waters) > 0:
            # Don't perform conversion if the topology already matches the
            # the template for the desired format.

            if format == "AMBER":
                if waters[0].isAmberWater():
                    return

            elif format == "GROMACS":
                if waters[0].isGromacsWater():
                    return

            # There will be a "water_model" system property if this object was
            # solvated by BioSimSpace.
            if "water_model" in self._sire_object.propertyKeys():
                water_model = self._sire_object.property("water_model").toString()

            # Otherwise, convert to an appropriate topology.
            else:
                num_point = waters[0].nAtoms()

                if num_point == 3:
                    # TODO: Assume TIP3P. Not sure how to detect SPC/E at present.
                    water_model = "TIP3P"
                elif num_point == 4:
                    water_model = "TIP4P"
                elif num_point == 5:
                    water_model = "TIP5P"

            if format == "AMBER":
                self._sire_object = _SireIO.setAmberWater(
                    self._sire_object, water_model, property_map
                )
            else:
                self._sire_object = _SireIO.setGromacsWater(
                    self._sire_object, water_model, property_map
                )

    def _set_atom_index_tally(self):
        """
        Internal helper function to create a dictionary mapping molecule
        numbers to the cumulative total of atoms in the system.
        """
        # Only compute the atom index mapping if it hasn't already
        # been created.
        if len(self._atom_index_tally) == 0:
            # Loop over all molecules in the system and keep track
            # of the cumulative number of atoms.
            num_atoms = 0
            for num in self._sire_object.molNums():
                self._atom_index_tally[num] = num_atoms
                num_atoms += self._sire_object.molecule(num).nAtoms()

    def _set_residue_index_tally(self):
        """
        Internal helper function to create a dictionary mapping molecule
        numbers to the cumulative total of atoms in the system.
        """
        # Only compute the residue index mapping if it hasn't already
        # been created.
        if len(self._residue_index_tally) == 0:
            # Loop over all molecules in the system and keep track
            # of the cumulative number of residues.
            num_residues = 0
            for num in self._sire_object.molNums():
                self._residue_index_tally[num] = num_residues
                num_residues += self._sire_object.molecule(num).nResidues()

    def _set_molecule_index_tally(self):
        """
        Internal helper function to create a dictionary mapping molecule
        numbers to the cumulative total of atoms in the system.
        """

        # Only compute the molecule index mapping if it hasn't already
        # been created.
        if len(self._molecule_index) == 0:
            # Loop over all molecules in the system by index and map
            # the MolNum to index.
            for idx in range(0, self.nMolecules()):
                self._molecule_index[
                    self._sire_object[_SireMol.MolIdx(idx)].number()
                ] = idx


# Import at bottom of module to avoid circular dependency.
from ._atom import Atom as _Atom
from ._molecule import Molecule as _Molecule
from ._molecules import Molecules as _Molecules
from ._residue import Residue as _Residue
from ._search_result import SearchResult as _SearchResult
