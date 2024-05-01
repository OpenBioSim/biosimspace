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

from ..._Exceptions import IncompatibleError as _IncompatibleError
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
        rmsd_indices,
        reference_index=None,
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

        reference : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            The reference molecule, against which the RMSD will be measured.
            This molecule should match with a single molecule from the
            system, i.e. contain the same residues as the matching molecule
            in the same order.

        rmsd_indices : [int]
            The indices of the atoms within the reference that should be
            used when calculating the RMSD. Atoms from the reference that
            are not in rmsd_indices will instead be used for alignment.

        reference_index : int
            The index of the molecule in the system that matches the
            reference. If none, then we will assume that the molecule with
            the closest number of residues is the match. If not all atoms
            from the reference are matched in the system, then an exception
            will be thrown.

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

        if not isinstance(reference, _Molecule):
            raise TypeError(
                "'reference' must be of type 'BioSimSpace._SireWrappers.Molecule'."
            )
        self._reference = reference.copy()

        # Check that the map is valid.
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        # Extract the corresponding molecule from the system by index.
        if reference_index is not None:
            if not type(reference_index) is int:
                raise TypeError("'reference_index' must be of type 'int'.")
            self._reference_index = reference_index
            molecule = system[reference_index]

        # Match to the molecule with the closest number of residues.
        else:
            # Store the number of residues in the target.
            num_res_ref = reference.nResidues()

            # Create a list to hold the number of residues for each molecule in the system.
            num_res_system = []
            for mol in system:
                num_res_system.append(mol.nResidues())

            # Find the index of the molecule with the closest number of residues.
            reference_index = min(
                enumerate(num_res_system), key=lambda x: abs(x[1] - num_res_ref)
            )[0]

            # Extract the molecule from the system.
            self._reference_index = reference_index
            molecule = system[reference_index]

        if isinstance(rmsd_indices, (list, tuple)):
            # Make sure the indices are unique.
            rmsd_indices = list(set(rmsd_indices))

            if len(rmsd_indices) == 0:
                raise ValueError("'rmsd_indices' contains no items?")
            elif len(rmsd_indices) >= reference.nAtoms():
                raise ValueError(
                    "'rmsd_indices' must refer to a subset of the "
                    "atoms from the reference molecule."
                )

            for idx in rmsd_indices:
                if not type(idx) is int:
                    raise TypeError("'rmsd_indices' must be a list of 'int' types.")
                if idx < 0:
                    idx += reference.nAtoms()
                if idx < 0 or idx >= reference.nAtoms():
                    raise ValueError(
                        "'rmsd_indices' is out of range: "
                        f"{-reference.nAtoms()}:{reference.nAtoms()-1}"
                    )
        else:
            raise TypeError("'rmsd_indices' must be a list of 'int' types.")

        # Create the matcher.
        matcher = _SireMol.ResNumAtomNameMatcher()

        # Match the atoms by residue number and name.
        matches = matcher.match(reference._sire_object, molecule._sire_object)

        # We need to match all of the atoms in the reference.
        if len(matches) < reference.nAtoms():
            raise _IncompatibleError(
                "Didn't match all of the atoms in the reference molecule: "
                f"Found {len(matches)}, expected {reference.nAtoms()}."
            )

        # Invert the matches so that we map from the system to the reference.
        matches = {v: k for k, v in matches.items()}

        # Get the indices of the matching atoms.
        idx_matches = matches.keys()

        # Extract the molecule of interest and make it editable.
        edit_mol = molecule._sire_object.edit()

        # A list of atoms to select.
        selected = []

        # Create a mapping dictionary for the atoms involved in the RMSD and alignment.
        rmsd_mapping = {}
        align_mapping = {}

        # Get the coordinates property from the user mapping.
        coord_prop = property_map.get("coordinates", "coordinates")

        # Loop over all atoms in the molecule.
        for x in range(0, molecule.nAtoms()):
            # Convert to an AtomIdx.
            idx = _SireMol.AtomIdx(x)
            # This atom was matched to an atom in the reference.
            if idx in idx_matches:
                # Append to the list of atoms to select.
                selected.append(idx)
                # Copy across the coordinates.
                edit_mol = (
                    edit_mol.atom(idx)
                    .setProperty(
                        coord_prop,
                        reference._sire_object.atom(matches[idx]).property(
                            "coordinates"
                        ),
                    )
                    .molecule()
                )
                # Set occupancy and beta factor.
                if matches[idx].value() in rmsd_indices:
                    # This atom is used to compute the RMSD.
                    edit_mol = (
                        edit_mol.atom(idx).setProperty("occupancy", 0.0).molecule()
                    )
                    edit_mol = (
                        edit_mol.atom(idx).setProperty("beta_factor", 1.0).molecule()
                    )
                    rmsd_mapping[idx] = matches[idx]
                else:
                    # This atom is used for alignment.
                    edit_mol = (
                        edit_mol.atom(idx).setProperty("occupancy", 1.0).molecule()
                    )
                    edit_mol = (
                        edit_mol.atom(idx).setProperty("beta_factor", 0.0).molecule()
                    )
                    align_mapping[idx.value()] = matches[idx].value()

        # Store the initial value of the RMSD. This is useful to use as a starting
        # point for the restraint when performing steered molecular dynamics.
        self._initial_value = self._compute_initial_rmsd(
            system,
            reference,
            reference_index,
            rmsd_mapping,
            align_mapping,
            property_map,
        )

        # Commit the changes.
        new_molecule = edit_mol.commit()

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

        # Save the changes to the molecule.
        molecule = _Molecule(partial_molecule)

        # Parse as a PDB file and store the lines.
        pdb = _SireIO.PDB2(molecule.toSystem()._sire_object)
        lines = pdb.toLines()

        # Format for PLUMED, making sure to use the same indices as in the system.
        new_indices = [idx.value() + 1 for idx in selected]
        self._reference_pdb = [
            line[:6] + str(idx).rjust(5) + line[11:]
            for line, idx in zip(lines[1:-2], new_indices)
        ]
        self._reference_pdb.append(lines[-1])

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
        string += " reference=%s" % self._reference
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

    def getReferenceIndex(self):
        """
        Return the index of the molecule in the system that corresponds
        to the reference.

        Returns
        -------

        reference_index : int
            The index of the molecule in the system that corresponds to
            the reference.
        """
        return self._reference_index

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
        reference_index,
        rmsd_mapping,
        align_mapping,
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

        reference_index : int
            The index of the molecule in the system that matches the
            reference.

        rmsd_mapping : {Sire.Mol.AtomIdx : Sire.Mol.AtomIdx, ...}
             A dictionary mapping atoms in the system to those in the
             reference molecule.

        align_mapping : { int: int, ...}
             A dictionary mapping atoms in the system to those in the
             reference molecule that will be used for alignment prior
             to computing the RMSD.

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

        if not isinstance(reference, _Molecule):
            raise TypeError(
                "'reference' must be of type 'BioSimSpace._SireWrappers.Molecule'."
            )

        if not type(reference_index) is int:
            raise TypeError("'reference_index' must be of type 'int'.")

        if not isinstance(rmsd_mapping, dict):
            raise TypeError("'rmsd_mapping' must be of type 'dict'.")

        for k, v in rmsd_mapping.items():
            if not isinstance(k, _SireMol.AtomIdx) or not isinstance(
                v, _SireMol.AtomIdx
            ):
                raise TypeError(
                    "'mapping' must contain 'Sire.Mol.AtomIdx' key-value pairs!"
                )

        if not isinstance(align_mapping, dict):
            raise TypeError("'align_mapping' must be of type 'dict'.")

        for k, v in align_mapping.items():
            if not type(k) is int or not type(v) is int:
                raise TypeError("'align_mapping' must contain 'int' key-value pairs!")

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        try:
            molecule = system[reference_index]
        except:
            raise ValueError("Molecule at 'reference_index' not found in 'system'!")

        # Get the 'space' property from the system.
        try:
            space_prop = property_map.get("space", "space")
            space = system._sire_object.property(space_prop)
        except:
            raise ValueError(
                f"'system' has no '{space_prop}' property. Unable to compute RMSD!"
            )

        # Align the molecule to the reference using the alignment mapping.
        if len(align_mapping) > 0:
            try:
                molecule = _rmsdAlign(
                    molecule,
                    reference,
                    align_mapping,
                    property_map0=property_map,
                    property_map1=property_map,
                )
            except:
                ValueError(
                    "Unable to align 'molecule' to 'reference' based on 'align_mapping'."
                )

        # Set the user-define coordinates property.
        coord_prop = property_map.get("coordinates", "coordinates")

        # Loop over all atom matches and compute the squared distance.
        dist2 = 0
        for idx0, idx1 in rmsd_mapping.items():
            try:
                coord0 = molecule._sire_object.atom(idx0).property(coord_prop)
                coord1 = reference._sire_object.atom(idx1).property(coord_prop)
            except:
                raise ValueError(
                    "Could not calculate initial RMSD due to missing coordinates!"
                )
            dist2 += space.calcDist2(coord0, coord1)

        # Compute the RMSD.
        dist2 /= len(rmsd_mapping)
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
