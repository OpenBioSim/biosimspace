######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2024
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

"""Functionality for a root-mean-square deviation collective variable."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["RMSD"]

from math import ceil as _ceil
from math import sqrt as _sqrt

from sire.legacy import IO as _SireIO
from sire.legacy import Mol as _SireMol

from sire.mol import selection_to_atoms as _selection_to_atoms

from ... import _isVerbose
from ..._Exceptions import IncompatibleError as _IncompatibleError
from ..._SireWrappers import Atom as _Atom
from ..._SireWrappers import Molecule as _Molecule
from ..._SireWrappers import System as _System
from ...Align import rmsdAlign as _rmsdAlign

from ._collective_variable import CollectiveVariable as _CollectiveVariable
from .._bound import Bound as _Bound
from .._grid import Grid as _Grid
from ...Types import Length as _Length


class RMSD(_CollectiveVariable):
    """A class for a root-mean-square deviation (RMSD) collective variable."""

    def __init__(
        self,
        system,
        reference,
        align_selection,
        rmsd_selection,
        reference_mapping=None,
        hill_width=_Length(0.1, "nanometer"),
        lower_bound=None,
        upper_bound=None,
        grid=None,
        alignment_type="optimal",
        pbc=True,
        property_map={},
    ):
        """
        Constructor.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system of interest.

        reference : :class:`System <BioSimSpace._SireWrappers.System>`
            The reference system, against which the RMSD will be measured.

        align_selection : str
            A Sire selection string that defines the atoms to be used
            when aligning the two structures. If None, then the RMSD
            will be calculated without alignment.

        rmsd_selection : str
            A Sire selection string that defines the atoms to be used
            when calculating the RMSD.

        reference_mapping : dict
            A dictionary mapping molecule indices in the reference to
            those in the system. This must be used when the reference
            represents a sub-set of the system.

        hill_width : :class:`Length <BioSimSpace.Types.Length>`
            The width of the Gaussian hill used to sample this variable.

        lower_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
            A lower bound on the value of the collective variable.

        upper_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
            An upper bound on the value of the collective variable.

        grid : :class:`Grid <BioSimSpace.Metadynamics.Grid>`
            The grid on which the collective variable will be sampled.
            This can help speed up long metadynamics simulations where
            the number of Gaussian kernels can become prohibitive.

        alignment_type : str
            The mannier in which RMSD alignment is performed. Options are
            "optimal" or "simple".

        pbc : bool
            Whether to use periodic boundary conditions when computing the
            collective variable.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__()

        # Set the types associated with this collective variable.
        self._types = [_Length]

        # Validate input.

        if not isinstance(system, _System):
            raise TypeError(
                "'system' must be of type 'BioSimSpace._SireWrappers.System'."
            )
        self._system = system.copy()

        if not isinstance(reference, _System):
            raise TypeError(
                "'reference' must be of type 'BioSimSpace._SireWrappers.System'."
            )
        self._reference = reference.copy()

        from sire.system import System

        sire_reference = System(reference._sire_object)

        if reference_mapping is None:
            # Make sure that the reference system is compatible with the system.
            if system.nMolecules() != reference.nMolecules():
                raise _IncompatibleError(
                    "The number of molecules in 'system' and 'reference' must match."
                )
            if system.nResidues() != reference.nResidues():
                raise _IncompatibleError(
                    "The number of residues in 'system' and 'reference' must match."
                )
            if system.nAtoms() != reference.nAtoms():
                raise _IncompatibleError(
                    "The number of atoms in 'system' and 'reference' must match."
                )

            self._reference_mapping = {i: i for i in range(system.nMolecules())}
        else:
            if not isinstance(reference_mapping, dict):
                raise TypeError("'reference_mapping' must be of type 'dict'")

            # Validate the mapping.
            for k, v in reference_mapping.items():
                if not isinstance(k, int):
                    raise TypeError("Keys of 'reference_mapping' must be of type 'int'")
                if not isinstance(v, int):
                    raise TypeError(
                        "Values of 'reference_mapping' must be of type 'int'"
                    )

                # Make sure the mapped molecules have the same number of atoms.
                if reference[k].nAtoms() != system[v].nAtoms():
                    raise _IncompatibleError(
                        "Molecules mapped via the 'reference_mapping' must have the same number of atoms."
                    )

            self._reference_mapping = reference_mapping

        # Validate alignment selection string.
        if not isinstance(align_selection, str):
            raise TypeError("'align_selection' must be of type 'str'")
        self._align_selection = align_selection

        try:
            self._aligment_atoms = _selection_to_atoms(sire_reference, align_selection)
        except Exception as e:
            msg = "Invalid 'align_selection' string."
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Validate RMSD selection string.
        if not isinstance(rmsd_selection, str):
            raise TypeError("'rmsd_selection' must be of type 'str'")
        self._rmsd_selection = rmsd_selection

        try:
            self._rmsd_atoms = _selection_to_atoms(sire_reference, rmsd_selection)
        except Exception as e:
            msg = "Invalid 'rmsd_selection' string."
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Check that the map is valid.
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        # Work out the unique molecule numbers used for alignment and RMSD.

        mol_nums = set()

        for atom in self._aligment_atoms:
            if atom.molecule().number() not in mol_nums:
                mol_nums.add(atom.molecule().number())

        for atom in self._rmsd_atoms:
            if atom.molecule().number() not in mol_nums:
                mol_nums.add(atom.molecule().number())

        # List to store the molecule indices in the system.
        self._molecule_indices = []

        # List to store index pairs for the mapped molecule (system, reference).
        molecule_pairs = []

        # List to store the absolute atom indices.
        abs_atom_indices = []

        # Dictionary to store alignment and RMSD indices for each molecule.
        align_indices = {}
        rmsd_indices = {}

        # Create mappings for the alignment and RMSD.

        # Loop over the molecules numbers.
        for num in mol_nums:
            # Extract the molecule from the reference system.
            molecule = reference._sire_object[num]

            # Work out the index of the molecule in the reference.
            index = reference.getIndex(_Molecule(molecule))

            # Map the index to the system.
            self._molecule_indices.append(self._reference_mapping[index])

            # Store the index pair.
            molecule_pairs.append((self._reference_mapping[index], index))

            # Set of atoms to select.
            selected = set()

            # Create a cursor for editing the molecule.
            cursor = molecule.cursor()

            # Loop over the atoms.
            for i, atom in enumerate(molecule.atoms()):
                is_align = False
                is_rmsd = False
                # This atom is used for alignment.
                if atom in self._aligment_atoms:
                    is_align = True
                    try:
                        align_indices[num].append(atom.index())
                    except:
                        align_indices[num] = [atom.index()]
                # This atom is used for RMSD.
                if atom in self._rmsd_atoms:
                    is_rmsd = True
                    try:
                        rmsd_indices[num].append(atom.index())
                    except:
                        rmsd_indices[num] = [atom.index()]

                if is_align or is_rmsd:
                    # Append to the list of atoms to select.
                    selected.add(atom.index())

                    # Add the absolute atom index.
                    abs_atom_indices.append(1 + reference.getIndex(_Atom(atom)))

                    # Set occupancy and beta factor.
                    if is_align:
                        cursor.atom(i)["occupancy"] = 1.0
                    else:
                        cursor.atom(i)["occupancy"] = 0.0
                    if is_rmsd:
                        cursor.atom(i)["beta_factor"] = 1.0
                    else:
                        cursor.atom(i)["beta_factor"] = 0.0

            # Commit the changes.
            new_molecule = cursor.commit()

            # Create an AtomSelection.
            selection = new_molecule.selection()

            # Unselect all of the atoms.
            selection.selectNone()

            # Now add all of the atoms that appear in the reference.
            for idx in selected:
                selection.select(idx)

            # Create a partial molecule and extract the atoms.
            partial_molecule = (
                _SireMol.PartialMolecule(new_molecule, selection).extract().molecule()
            )

            # Update the new system.
            try:
                new_system += _Molecule(partial_molecule)
            except:
                new_system = _Molecule(partial_molecule).toSystem()

        # Parse as a PDB file and store the lines.
        pdb = _SireIO.PDB2(new_system._sire_object)
        lines = pdb.toLines()

        # Format for PLUMED, making sure to use the same indices as in the system.
        # Also strip any TER records.
        self._reference_pdb = []
        for line, idx in zip(lines[1:-2], abs_atom_indices):
            if not "TER" in line:
                self._reference_pdb.append(line[:6] + str(idx).rjust(5) + line[11:])
        self._reference_pdb.append(lines[-1])

        # Store the initial value of the RMSD. This is useful to use as a starting
        # point for the restraint when performing steered molecular dynamics.
        self._initial_value = self._compute_initial_rmsd(
            system,
            reference,
            molecule_pairs,
            align_indices,
            rmsd_indices,
            property_map,
        )

        # Set the "settable" parameters.
        self.setHillWidth(hill_width)
        self.setAlignmentType(alignment_type)
        self.setPeriodicBoundaries(pbc)

        # Set defaults for optional values.
        self._lower_bound = None
        self._upper_bound = None
        self._grid = None

        # Set the optional parameters.
        if lower_bound is not None:
            self.setLowerBound(lower_bound)
        if upper_bound is not None:
            self.setUpperBound(upper_bound)
        if grid is not None:
            self.setGrid(grid)

        # Validate that the state is self-consistent.
        self._validate()

        # Flag that the object has been instantiated, i.e. it is no longer "new".
        self._is_new_object = False

    def __str__(self):
        """Return a human readable string representation of the object."""
        string = "<BioSimSpace.Metadynamics.CollectiveVariable.RMSD: "
        string += " align_selection=%s" % self._align_selection
        string += ", rmsd_selection=%s" % self._rmsd_selection
        string += ", hill_width=%s" % self._hill_width
        if self._lower_bound is not None:
            string += ", lower_bound=%s" % self._lower_bound
        if self._upper_bound is not None:
            string += ", upper_bound=%s" % self._upper_bound
        if self._grid is not None:
            string += ", grid=%s" % self._grid
        string += ", alignment_type=%s" % self._alignment_type
        string += ", pbc=%s" % self._pbc
        string += ">"
        return string

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return self.__str__()

    def __eq__(self, other):
        """Equality operator."""
        return (
            self._reference == other._reference
            and self._hill_width == other._hill_width
            and self._lower_bound == other._lower_bound
            and self._upper_bound == other._upper_bound
            and self._grid == other._grid
            and self._alignment_type == other._alignment_type
            and self._pbc == other._pbc
        )

    def getReferencePDB(self):
        """
        Return the reference PDB file as a list of strings.

        Returns
        -------

        pdb : [str]
            The reference PDB file as list of strings.
        """
        return self._reference_pdb

    def getInitialValue(self):
        """
        Return the initial value of the collective variable.

        Returns
        -------

        rmsd : :class:`Length <BioSimSpace.Types.Length>`
            The initial value of the collective variable.
        """
        return self._initial_value

    def setHillWidth(self, hill_width):
        """
        Set the width of the Gaussian hills used to bias this collective
        variable.

        hill_width : :class:`Length <BioSimSpace.Types.Length>`
            The width of the Gaussian hill.
        """
        if not isinstance(hill_width, _Length):
            raise TypeError("'hill_width' must be of type 'BioSimSpace.Types.Length'")

        if hill_width.value() < 0:
            raise ValueError("'hill_width' must have a value of > 0")

        # Convert to the internal unit.
        self._hill_width = hill_width.nanometers()

    def getHillWidth(self):
        """
        Return the width of the Gaussian hill used to bias this collective
        variable.

        Returns
        -------

        hill_width : :class:`Length <BioSimSpace.Types.Length>`
            The width of the Gaussian hill.
        """
        return self._hill_width

    def setAlignmentType(self, alignment_type):
        """
        Set the RMSD alignment type. Options are "optimal" or "simple".

        Parameters
        ----------

        alignment_type : str
            The RMSD alignment type.
        """
        if not isinstance(alignment_type, str):
            raise TypeError("'alignment_type' must be of type 'str'")

        alignment_type = alignment_type.lower().replace(" ", "")
        if alignment_type not in ["optimal", "simple"]:
            raise ValueError("'alignment_type' must be 'optimal' or 'simple'.")

        self._alignment_type = alignment_type

    def getAlignmentType(self):
        """
        Return the RMSD alignment type.

        Returns
        -------

        alignment_type : str
            The RMSD alignment type.
        """
        return self._alignment_type

    def getMoleculeIndices(self):
        """
        Return the indices of molecules involved in the collective variable.

        Returns
        -------

        molecule_indices : int
            The indices of molecules involved in the collective variable.
        """
        return self._molecule_indices

    def setPeriodicBoundaries(self, pbc):
        """
        Set whether to use periodic_boundaries when calculating the
        collective variable.

        Parameters
        ----------

        pbc : bool
            Whether to use periodic boundaries conditions.
        """
        if not isinstance(pbc, bool):
            raise TypeError("'pbc' must be of type 'bool'")
        self._pbc = pbc

    def getPeriodicBoundaries(self):
        """
        Return whether to take account of periodic boundary conditions
        when computing the collective variable.

        Returns
        -------

        pbc : bool
            Whether to use periodic boundaries conditions.
        """
        return self._pbc

    def _compute_initial_rmsd(
        self,
        system,
        reference,
        molecule_pairs,
        align_indices,
        rmsd_indices,
        property_map={},
    ):
        """
        Compute the initial value of the RMSD collective variable.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system of interest.

        reference : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            The reference molecule, against which the RMSD will be measured.
            This molecule should match with a single molecule from the
            system, i.e. contain the same residues as the matching molecule
            in the same order.

        molecule_pairs : [(int, int), ...]
            The indices of molecules in the system and reference that contain
            atoms involved in alignment and RMSD.

        align_indices : {Sire.Mol.MolNum: [Sire.Mol.AtomIdx, ...]}
            A dictionary mapping molecules to the indices of atoms that will
            be used for alignment.

        rmsd_indices : {Sire.Mol.MolNum: [Sire.Mol.AtomIdx, ...]}
            A dictionary mapping molecules to the indices of atoms that will
            be used for the RMSD.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

         Returns
         -------

         rmsd : :class:`Length <BioSimSpace.Types.Length>`
             The initial value of the RMSD.
        """

        # Note that we need to do this manually, since Sire.Mol.Evaluator doesn't
        # work correctly for molecules with different numbers of coordinate groups.

        if not isinstance(system, _System):
            raise TypeError(
                "'system' must be of type 'BioSimSpace._SireWrappers.System'."
            )

        if not isinstance(reference, _System):
            raise TypeError(
                "'reference' must be of type 'BioSimSpace._SireWrappers.System'."
            )

        if not isinstance(molecule_pairs, list):
            raise TypeError("'molecule_pairs' must be a list of integer tuples.")
        for pair in molecule_pairs:
            if not isinstance(pair, tuple):
                raise TypeError("'molecule_pairs' must be a list of integer tuples.")
            if len(pair) != 2:
                raise ValueError("'molecule_pairs' must be a list of integer tuples.")
            if not isinstance(pair[0], int):
                raise TypeError("'molecule_pairs' must be a list of integer tuples.")
            if not isinstance(pair[1], int):
                raise TypeError("'molecule_pairs' must be a list of integer tuples.")

        if not isinstance(align_indices, dict):
            raise TypeError("'align_indices' must be a dictionary.")
        for key, value in align_indices.items():
            if not isinstance(key, _SireMol.MolNum):
                raise TypeError(
                    "Keys of 'align_indices' must be of type 'sire.legacy.Mol.MolMolNum'."
                )
            if not isinstance(value, list):
                raise TypeError(
                    "Values of 'align_indices' must be lists of 'sire.legacy.Mol.AtomIdx' types."
                )
            for idx in value:
                if not isinstance(idx, _SireMol.AtomIdx):
                    raise TypeError(
                        "Values of 'align_indices' must be lists of 'sire.legacy.Mol.AtomIdx' types."
                    )

        if not isinstance(rmsd_indices, dict):
            raise TypeError("'rmsd_indices' must be a dictionary.")
        for key, value in rmsd_indices.items():
            if not isinstance(key, _SireMol.MolNum):
                raise TypeError(
                    "Keys of 'rmsd_indices' must be of type 'sire.legacy.Mol.MolNum'."
                )
            if not isinstance(value, list):
                raise TypeError(
                    "Values of 'rmsd_indices' must be lists of 'sire.legacy.Mol.AtomIdx' types."
                )
            for idx in value:
                if not isinstance(idx, _SireMol.AtomIdx):
                    raise TypeError(
                        "Values of 'rmsd_indices' must be lists of 'sire.legacy.Mol.AtomIdx' types."
                    )

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        # Get the 'space' property from the system.
        try:
            space_prop = property_map.get("space", "space")
            space = system._sire_object.property(space_prop)
        except:
            raise ValueError(
                f"'system' has no '{space_prop}' property. Unable to compute RMSD!"
            )

        # Set the user-define coordinates property.
        coord_prop = property_map.get("coordinates", "coordinates")

        # Total squared distance.
        dist2 = 0

        # Total number of RMSD atoms.
        num_rmsd = 0

        # Loop over the molecules.
        for idx_system, idx_ref in molecule_pairs:
            mol = system[idx_system]
            ref = reference[idx_ref]

            align_mapping = {}

            try:
                align_mapping = {
                    i.value(): i.value()
                    for i in align_indices[ref._sire_object.number()]
                }
            except Exception as e:
                pass

            if len(align_mapping) > 0:
                try:
                    new_mol = _rmsdAlign(
                        mol,
                        ref,
                        align_mapping,
                        property_map0=property_map,
                        property_map1=property_map,
                    )
                except:
                    ValueError(
                        "Unable to align 'molecule' to 'reference' based on 'align_mapping'."
                    )

            rmsd_mapping = {}

            try:
                rmsd_mapping = {i: i for i in rmsd_indices[ref._sire_object.number()]}
            except Exception as e:
                pass

            if len(rmsd_mapping) > 0:
                # Loop over all atom matches and compute the squared distance.
                for idx0, idx1 in rmsd_mapping.items():
                    try:
                        coord0 = new_mol._sire_object.atom(idx0).property(coord_prop)
                        coord1 = ref._sire_object.atom(idx1).property(coord_prop)
                    except:
                        raise ValueError(
                            "Could not calculate initial RMSD due to missing coordinates!"
                        )
                    dist2 += space.calcDist2(coord0, coord1)
                    num_rmsd += 1

        # Compute the RMSD.
        dist2 /= num_rmsd
        rmsd = _sqrt(dist2)

        return _Length(rmsd, "Angstrom")

    def _validate(self):
        """Internal function to check that the object is in a consistent state."""

        if self._lower_bound is not None:
            if type(self._lower_bound.getValue()) not in self._types:
                raise TypeError(
                    "'lower_bound' must be of type 'BioSimSpace.Types.Length'"
                )
        if self._upper_bound is not None:
            if type(self._upper_bound.getValue()) not in self._types:
                raise TypeError(
                    "'upper_bound' must be of type 'BioSimSpace.Types.Length'"
                )
        if self._lower_bound is not None and self._upper_bound is not None:
            if self._lower_bound.getValue() >= self._upper_bound.getValue():
                raise TypeError("'lower_bound' must less than 'upper_bound'")

        if self._grid is not None:
            if type(self._grid.getMinimum()) not in self._types:
                raise TypeError(
                    "'grid' minimum must be of type 'BioSimSpace.Types.Length'"
                )
            if type(self._grid.getMaximum()) not in self._types:
                raise TypeError(
                    "Grid 'maximum' must be of type 'BioSimSpace.Types.Length'"
                )
            if (
                self._lower_bound is not None
                and self._grid.getMinimum() > self._lower_bound.getValue()
            ):
                raise ValueError("'lower_bound' is less than 'grid' minimum.")
            if (
                self._upper_bound is not None
                and self._grid.getMaximum() < self._upper_bound.getValue()
            ):
                raise ValueError("'upper_bound' is greater than 'grid' maximum.")

            # If the number of bins isn't specified, estimate it out from the hill width.
            if self._grid.getBins() is None:
                grid_range = (self._grid.getMaximum() - self._grid.getMinimum()).value()
                num_bins = _ceil(5.0 * (grid_range / self._hill_width.value()))
                self._grid.setBins(num_bins)
