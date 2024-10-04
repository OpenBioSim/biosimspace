######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2024
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#          Matthew Burman <matthew@openbiosim.org>
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
######################################################################

# Functionality for creating and viewing systems for Atomic transfer.

__all__ = ["AToM"]

import copy as _copy
import json as _json
import os as _os
import shutil as _shutil
import warnings as _warnings

from sire.legacy import IO as _SireIO

from .._SireWrappers import Molecule as _Molecule
from .._SireWrappers import System as _System
from .. import _Utils
from ..Types import Length as _Length
from ..Types import Vector as _Vector
from ..Types import Coordinate as _Coordinate
from ..Align import matchAtoms as _matchAtoms
from ..Align import rmsdAlign as _rmsdAlign
from ..Notebook import View as _View
from .. import _isVerbose
from ..Process import OpenMM as _OpenMM
from ..Process import ProcessRunner as _ProcessRunner


class AToM:
    """A class for creating, setting up, running and analysing AToM
    simulations."""

    def __init__(
        self,
        system=None,
        receptor=None,
        ligand_bound=None,
        ligand_free=None,
        protein_index=0,
        ligand_bound_index=1,
        ligand_free_index=2,
    ):
        """Constructor for the AToM class.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            A pre-prepared AToM system containing protein and ligands placed in their correct positions.
            If provided takes precedence over protein, ligand_bound and ligand_free.

        receptor : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            A receptor molecule. Will be used along with ligand_bound and ligand_free to create a system.

        ligand_bound : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            The bound ligand. Will be used along with protein and ligand_free to create a system.

        ligand_free : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            The free ligand. Will be used along with protein and ligand_bound to create a system.

        protein_index : int, [int]
            If passing a pre-prepared system, the index (or indices) of the protein molecule in the system (Default 0).

        ligand_bound_index : int
            If passing a pre-prepared system, the index of the bound ligand molecule in the system (Default 1).

        ligand_free_index : int
            If passing a pre-prepared system, the index of the free ligand molecule in the system (Default 2).
        """
        # make sure that either system or protein, ligand_bound and ligand_free are given
        if system is None and not all(
            x is not None for x in [receptor, ligand_bound, ligand_free]
        ):
            raise ValueError(
                "Either a pre-prepared system or protein, bound ligand and free ligand must be given."
            )
        # check that the system is a BioSimSpace system
        # or the other inputs are BioSimSpace molecules
        if system is not None and not isinstance(system, _System):
            raise ValueError("The system must be a BioSimSpace System object.")
        elif not all(
            isinstance(x, _Molecule)
            for x in [receptor, ligand_bound, ligand_free]
            if x is not None
        ):
            raise ValueError(
                "The protein, bound ligand and free ligand must be BioSimSpace Molecule objects."
            )
        self._is_prepared = False
        self._setSystem(system)
        if not self._is_prepared:
            self._setProtein(receptor)
            self._setLigandBound(ligand_bound)
            self._setLigandFree(ligand_free)
        else:
            self._setProteinIndex(protein_index)
            self._setLigandBoundIndex(ligand_bound_index)
            self._setLigandFreeIndex(ligand_free_index)

    def _setSystem(self, system, is_prepped=True):
        """Set the system for the AToM simulation.

        Parameters
        ----------

        system : BioSimSpace._SireWrappers.System
            The system for the AToM simulation.
        """
        if system is not None:
            if not isinstance(system, _System):
                raise ValueError(
                    f"The system must be a BioSimSpace System object. It is currently {type(system)}."
                )
            elif len(system.getMolecules()) < 3:
                raise ValueError(
                    "The system must contain at least three molecules (a protein and two ligands)."
                )
            else:
                self._system = system
                self._is_prepared = is_prepped
        else:
            self._system = None
            self._is_prepared = False

    def _getSystem(self):
        """Get the system for the AToM simulation.

        Returns
        -------
        BioSimSpace._SireWrappers.System
            The system for the AToM simulation.
        """
        return self._system

    def _setProtein(self, protein):
        """Set the protein for the AToM simulation.

        Parameters
        ----------

        protein : BioSimSpace._SireWrappers.Molecule
            The protein for the AToM simulation.
        """
        if protein is not None:
            if not isinstance(protein, _Molecule):
                raise ValueError("The protein must be a BioSimSpace Molecule object.")
            else:
                self._protein = protein
        else:
            self._protein = None

    def _getProtein(self):
        """Get the protein for the AToM simulation.

        Returns
        -------
        BioSimSpace._SireWrappers.Molecule
            The protein for the AToM simulation.
        """
        return self._protein

    def _setLigandBound(self, ligand_bound):
        """Set the bound ligand for the AToM simulation.

        Parameters
        ----------

        ligand_bound : BioSimSpace._SireWrappers.Molecule
            The bound ligand for the AToM simulation.
        """
        if ligand_bound is not None:
            if not isinstance(ligand_bound, _Molecule):
                raise ValueError(
                    "The bound ligand must be a BioSimSpace Molecule object."
                )
            else:
                self._ligand_bound = ligand_bound
        else:
            self._ligand_bound = None

    def _getLigandBound(self):
        """Get the bound ligand for the AToM simulation.

        Returns
        -------
        BioSimSpace._SireWrappers.Molecule
            The bound ligand for the AToM simulation.
        """
        return self._ligand_bound

    def _setLigandFree(self, ligand_free):
        """Set the free ligand for the AToM simulation.

        Parameters
        ----------

        ligand_free : BioSimSpace._SireWrappers.Molecule
            The free ligand for the AToM simulation.
        """
        if ligand_free is not None:
            if not isinstance(ligand_free, _Molecule):
                raise ValueError(
                    "The free ligand must be a BioSimSpace Molecule object."
                )
            else:
                self._ligand_free = ligand_free
        else:
            self._ligand_free = None

    def _getLigandFree(self):
        """Get the free ligand for the AToM simulation.

        Returns
        -------
        BioSimSpace._SireWrappers.Molecule
            The free ligand for the AToM simulation.
        """
        return self._ligand_free

    def _setDisplacement(self, displacement):
        """Set the displacement of the free ligand along the normal vector."""
        if isinstance(displacement, str):
            try:
                self.displacement = _Length(displacement)
            except Exception as e:
                raise ValueError(
                    f"Could not convert {displacement} to a BSS length, due to the following error: {e}"
                )
        elif isinstance(displacement, _Length):
            self.displacement = displacement
        elif isinstance(displacement, list):
            if len(displacement) != 3:
                raise ValueError("displacement must have length 3")
            if all(isinstance(x, (float, int)) for x in displacement):
                self.displacement = _Vector(*displacement)
            elif all(isinstance(x, _Length) for x in displacement):
                self.displacement = _Vector([x.value() for x in displacement])
            else:
                raise TypeError("displacement must be a list of floats or BSS lengths")
        elif isinstance(displacement, _Vector):
            self.displacement = displacement
        else:
            raise TypeError(
                f"displacement must be a string, BSS length or list. It is currently {type(displacement)}."
            )
        if self._is_prepared:
            if not isinstance(self.displacement, _Vector):
                raise ValueError(
                    "Displacement must be a vector or list if a pre-prepared system is given"
                )

    def _getDisplacement(self):
        """Get the displacement of the free ligand along the normal vector.

        Returns
        -------
        BioSimSpace.Types.Length
            The displacement of the free ligand along the normal vector.
        """
        return self.displacement

    def _setLigandBoundRigidCore(self, ligand_bound_rigid_core):
        """Set the indices for the rigid core atoms of ligand 1.

        Parameters
        ----------

        ligand_bound_rigid_core : BioSimSpace._SireWrappers.Molecule
            The rigid core of the bound ligand for the AToM simulation.
        """
        if ligand_bound_rigid_core is None:
            self.ligand_bound_rigid_core = None
        else:
            if not isinstance(ligand_bound_rigid_core, list):
                raise TypeError("ligand_bound_rigid_core must be a list")
            if len(ligand_bound_rigid_core) != 3:
                raise ValueError("ligand_bound_rigid_core must have length 3")
            # make sure all indices are ints
            if not all(isinstance(x, int) for x in ligand_bound_rigid_core):
                raise TypeError("ligand_bound_rigid_core must contain only integers")
            if any(x >= self.ligand_bound_atomcount for x in ligand_bound_rigid_core):
                raise ValueError(
                    "ligand_bound_rigid_core contains an index that is greater than the number of atoms in the ligand"
                )
            self.ligand_bound_rigid_core = ligand_bound_rigid_core

    def _getLigandBoundRigidCore(self):
        """Get the indices for the rigid core atoms of ligand 1.

        Returns
        -------
        list
            The indices for the rigid core atoms of ligand 1.
        """
        return self.ligand_bound_rigid_core

    def _setLigandFreeRigidCore(self, ligand_free_rigid_core):
        """Set the indices for the rigid core atoms of ligand 2.

        Parameters
        ----------

        ligand_free_rigid_core : BioSimSpace._SireWrappers.Molecule
            The rigid core of the free ligand for the AToM simulation.
        """
        if ligand_free_rigid_core is None:
            self.ligand_free_rigid_core = None
        else:
            if not isinstance(ligand_free_rigid_core, list):
                raise TypeError("ligand_free_rigid_core must be a list")
            if len(ligand_free_rigid_core) != 3:
                raise ValueError("ligand_free_rigid_core must have length 3")
            # make sure all indices are ints
            if not all(isinstance(x, int) for x in ligand_free_rigid_core):
                raise TypeError("ligand_free_rigid_core must contain only integers")
            if any(x >= self.ligand_free_atomcount for x in ligand_free_rigid_core):
                raise ValueError(
                    "ligand_free_rigid_core contains an index that is greater than the number of atoms in the ligand"
                )
            self.ligand_free_rigid_core = ligand_free_rigid_core

    def _getLigandFreeRigidCore(self):
        """Get the indices for the rigid core atoms of ligand 2.

        Returns
        -------
        list
            The indices for the rigid core atoms of ligand 2.
        """
        return self.ligand_free_rigid_core

    def _setProteinIndex(self, protein_index):
        """
        Set the index of the protein in the system

        Parameters
        ----------

        protein_index : list
            The index or indices of the protein in the system.
        """
        if isinstance(protein_index, list):
            # check that all elements are ints
            if not all(isinstance(x, int) for x in protein_index):
                raise TypeError("protein_index must be a list of ints or a single int")
            for p in protein_index:
                if p < 0:
                    raise ValueError("protein_index must be a positive integer")
                if self._system[p].isWater():
                    _warnings.warn(
                        f"The molecule at index {p} is a water molecule, check your protein_index list."
                    )
            self.protein_index = protein_index
        elif isinstance(protein_index, int):
            self.protein_index = [protein_index]
        else:
            raise TypeError("protein_index must be an int or a list of ints")

    def _getProteinIndex(self):
        """Get the index of the protein molecule in the system.

        Returns
        -------
        int
            The index of the protein molecule in the system.
        """
        return self.protein_index

    def _setLigandBoundIndex(self, ligand_bound_index):
        """Set the index of the bound ligand molecule in the system.

        Parameters
        ----------

        ligand_bound_index : int
            The index of the bound ligand molecule in the system.
        """
        if not isinstance(ligand_bound_index, int):
            raise ValueError("ligand_bound_index must be an integer.")
        else:
            if ligand_bound_index < 0:
                raise ValueError("ligand_bound_index must be a positive integer")
            if self._system[ligand_bound_index].isWater():
                _warnings.warn(
                    f"The molecule at index {ligand_bound_index} is a water molecule, check your ligand_bound_index."
                )
            self.ligand_bound_index = ligand_bound_index

    def _getLigandBoundIndex(self):
        """Get the index of the bound ligand molecule in the system.

        Returns
        -------
        int
            The index of the bound ligand molecule in the system.
        """
        return self.ligand_bound_index

    def _setLigandFreeIndex(self, ligand_free_index):
        """Set the index of the free ligand molecule in the system.

        Parameters
        ----------

        ligand_free_index : int
            The index of the free ligand molecule in the system.
        """
        if not isinstance(ligand_free_index, int):
            raise ValueError("ligand_free_index must be an integer.")
        else:
            if ligand_free_index < 0:
                raise ValueError("ligand_free_index must be a positive integer")
            if self._system[ligand_free_index].isWater():
                _warnings.warn(
                    f"The molecule at index {ligand_free_index} is a water molecule, check your ligand_free_index."
                )
            self.ligand_free_index = ligand_free_index

    def _getLigandFreeIndex(self):
        """Get the index of the free ligand molecule in the system.

        Returns
        -------
        int
            The index of the free ligand molecule in the system.
        """
        return self.ligand_free_index

    def prepare(
        self,
        ligand_bound_rigid_core,
        ligand_free_rigid_core,
        displacement="20A",
        protein_com_atoms=None,
        ligand_bound_com_atoms=None,
        ligand_free_com_atoms=None,
    ):
        """Prepare the system for the AToM simulation.

        Parameters
        ----------

        ligand_bound_rigid_core : [int]
            A list of three atom indices that define the rigid core of the bound ligand.
            Indices are set relative to the ligand, not the system and are 0-indexed.

        ligand_free_rigid_core : [int]
            A list of three atom indices that define the rigid core of the free ligand.
            Indices are set relative to the ligand, not the system and are 0-indexed.

        displacement : float, string, [float, float, float]
            The diplacement between the bound and free ligands.
            If a float or string is given, BioSimSpace will attempt to find the ideal
            vector along which to displace the ligand by the given magnitude. If a list
            is given, the vector will be used directly.

        protein_com_atoms : [int]
            A list of atom indices that define the center of mass of the protein.
            If None, the center of mass of the protein will be found automatically.

        ligand_bound_com_atoms : [int]
            A list of atom indices that define the center of mass of the bound ligand.
            If None, the center of mass of the bound ligand will be found automatically.

        ligand_free_com_atoms : [int]
            A list of atom indices that define the center of mass of the free ligand.
            If None, the center of mass of the free ligand will be found automatically.

        Returns
        -------

        System : :class:`System <BioSimSpace._SireWrappers.System>`
            The prepared system, including protein and ligands in their correct positions.

        Data : dict
            A dictionary containing the data needed for the AToM simulation.
        """
        if self._is_prepared:
            self._systemInfo()
            self._setLigandBoundRigidCore(ligand_bound_rigid_core)
            self._setLigandFreeRigidCore(ligand_free_rigid_core)
            self._setDisplacement(displacement)
            self._setProtComAtoms(protein_com_atoms)
            self._setLig1ComAtoms(ligand_bound_com_atoms)
            self._setLig2ComAtoms(ligand_free_com_atoms)

            self._findAtomIndices()
            self._makeData()
            serialisable_disp = [
                self.displacement.x(),
                self.displacement.y(),
                self.displacement.z(),
            ]
            temp_data = self.data.copy()
            temp_data["displacement"] = serialisable_disp
            self._system._sire_object.setProperty("atom_data", _json.dumps(temp_data))
            return self._system, self.data

        else:
            # A bit clunky, but setDisplacement needs to be called twice - before and after _makeSystemFromThree
            # the final value will be set after the system is made, but the initial value is needed to make the system
            self._setDisplacement(displacement)
            system, prot_ind, lig1_ind, lig2_ind, dis_vec = self._makeSystemFromThree(
                self._protein, self._ligand_bound, self._ligand_free, self.displacement
            )
            self._setSystem(system, is_prepped=False)
            self._setDisplacement(dis_vec)
            self._setProteinIndex(prot_ind)
            self._setLigandBoundIndex(lig1_ind)
            self._setLigandFreeIndex(lig2_ind)
            self._systemInfo()
            self._setLigandBoundRigidCore(ligand_bound_rigid_core)
            self._setLigandFreeRigidCore(ligand_free_rigid_core)
            self._setProtComAtoms(protein_com_atoms)
            self._setLig1ComAtoms(ligand_bound_com_atoms)
            self._setLig2ComAtoms(ligand_free_com_atoms)
            self._findAtomIndices()
            self._makeData()
            serialisable_disp = [
                self.displacement.x(),
                self.displacement.y(),
                self.displacement.z(),
            ]
            temp_data = self.data.copy()
            temp_data["displacement"] = serialisable_disp
            # encode data in system for consistency
            self._system._sire_object.setProperty("atom_data", _json.dumps(temp_data))
            return self._system, self.data

    @staticmethod
    def _makeSystemFromThree(protein, ligand_bound, ligand_free, displacement):
        """Create a system for AToM simulations.

        Parameters
        ----------

        protein : BioSimSpace._SireWrappers.Molecule
            The protein for the AToM simulation.

        ligand_bound : BioSimSpace._SireWrappers.Molecule
            The bound ligand for the AToM simulation.

        ligand_free : BioSimSpace._SireWrappers.Molecule
            The free ligand for the AToM simulation.

        displacement : BioSimSpace.Types.Length
            The displacement of the ligand along the normal vector.

        Returns
        -------

        BioSimSpace._SireWrappers.System
            The system for the AToM simulation.
        """

        def _findTranslationVector(system, displacement, protein, ligand):

            from sire.legacy.Maths import Vector

            if not isinstance(system, _System):
                raise TypeError("system must be a BioSimSpace system")
            if not isinstance(protein, (_Molecule, type(None))):
                raise TypeError("protein must be a BioSimSpace molecule")
            if not isinstance(ligand, (_Molecule, type(None))):
                raise TypeError("ligand must be a BioSimSpace molecule")

            # Assume that binding sire is the center of mass of the ligand
            binding = _Coordinate(*ligand._getCenterOfMass())

            # Create grid around the binding site
            # This will act as the search region
            grid_length = _Length(20.0, "angstroms")

            num_edges = 5
            search_radius = (grid_length / num_edges) / 2
            grid_min = binding - 0.5 * grid_length
            grid_max = binding + 0.5 * grid_length

            non_protein_coords = Vector()
            # Count grid squares that contain no protein atoms
            num_non_prot = 0

            import numpy as np

            # Loop over the grid
            for x in np.linspace(grid_min.x().value(), grid_max.x().value(), num_edges):
                for y in np.linspace(
                    grid_min.y().value(), grid_max.y().value(), num_edges
                ):
                    for z in np.linspace(
                        grid_min.z().value(), grid_max.z().value(), num_edges
                    ):
                        search = (
                            f"atoms within {search_radius.value()} of ({x}, {y}, {z})"
                        )

                        try:
                            protein.search(search)
                        except:
                            non_protein_coords += Vector(x, y, z)
                            num_non_prot += 1

            non_protein_coords /= num_non_prot
            non_protein_coords = _Coordinate._from_sire_vector(non_protein_coords)

            # Now search out alpha carbons in system
            x = binding.x().angstroms().value()
            y = binding.y().angstroms().value()
            z = binding.z().angstroms().value()
            string = f"(atoms within 10 of {x},{y},{z}) and atomname CA"

            try:
                search = system.search(string)
            except:
                _warnings.warn(
                    "No alpha carbons found in system, falling back on any carbon atoms."
                )
                try:
                    string = f"(atoms within 10 of {x},{y},{z}) and element C"
                    search = system.search(string)
                except:
                    raise ValueError("No carbon atoms found in system")

            com = _Coordinate(_Length(0, "A"), _Length(0, "A"), _Length(0, "A"))
            atoms1 = []
            for atom in search:
                com += atom.coordinates()
                atoms1.append(system.getIndex(atom))
            com /= search.nResults()

            initial_normal_vector = (non_protein_coords - com).toVector().normalise()

            out_of_protein = displacement.value() * initial_normal_vector
            return out_of_protein

        mapping = _matchAtoms(ligand_free, ligand_bound)
        ligand_free_aligned = _rmsdAlign(ligand_free, ligand_bound, mapping)
        prot_lig1 = (protein + ligand_bound).toSystem()

        if isinstance(displacement, _Vector):
            ligand_free_aligned.translate(
                [displacement.x(), displacement.y(), displacement.z()]
            )
            vec = displacement
        else:
            vec = _findTranslationVector(prot_lig1, displacement, protein, ligand_bound)
            ligand_free_aligned.translate([vec.x(), vec.y(), vec.z()])

        sys = (protein + ligand_bound + ligand_free_aligned).toSystem()
        prot_ind = sys.getIndex(protein)
        lig1_ind = sys.getIndex(ligand_bound)
        lig2_ind = sys.getIndex(ligand_free_aligned)
        return sys, prot_ind, lig1_ind, lig2_ind, vec

    def _systemInfo(self):
        """
        If the user gives a pre-prepared AToM system, extract the needed information.
        """
        for p in self.protein_index:
            if self._system[p].isWater():
                _warnings.warn(
                    f"The molecule at index {self.protein_index} appears to be a water molecule."
                    " This should be a protein."
                )
        if self._system[self.ligand_bound_index].isWater():
            _warnings.warn(
                f"The molecule at index {self.ligand_bound_index} appears to be a water molecule."
                " This should be the bound ligand."
            )
        if self._system[self.ligand_free_index].isWater():
            _warnings.warn(
                f"The molecule at index {self.ligand_free_index} appears to be a water molecule."
                " This should be the free ligand."
            )
        self._protein_atomcount = sum(
            self._system[i].nAtoms() for i in self.protein_index
        )
        self.ligand_bound_atomcount = self._system[self.ligand_bound_index].nAtoms()
        self.ligand_free_atomcount = self._system[self.ligand_free_index].nAtoms()

    def _findAtomIndices(self):
        """
        Find the indices of the protein and ligand atoms in the system

        Returns
        -------
        dict
            A dictionary containing the indices of the protein and ligand atoms in the system
        """
        protein_atom_start = self._system[self.protein_index[0]].getAtoms()[0]
        protein_atom_end = self._system[self.protein_index[-1]].getAtoms()[-1]
        self.first_protein_atom_index = self._system.getIndex(protein_atom_start)
        self.last_protein_atom_index = self._system.getIndex(protein_atom_end)

        ligand_bound_atom_start = self._system[self.ligand_bound_index].getAtoms()[0]
        ligand_bound_atom_end = self._system[self.ligand_bound_index].getAtoms()[-1]
        self.first_ligand_bound_atom_index = self._system.getIndex(
            ligand_bound_atom_start
        )
        self.last_ligand_bound_atom_index = self._system.getIndex(ligand_bound_atom_end)

        ligand_free_atom_start = self._system[self.ligand_free_index].getAtoms()[0]
        ligand_free_atom_end = self._system[self.ligand_free_index].getAtoms()[-1]
        self.first_ligand_free_atom_index = self._system.getIndex(
            ligand_free_atom_start
        )
        self.last_ligand_free_atom_index = self._system.getIndex(ligand_free_atom_end)

    def _getProtComAtoms(self):
        """
        Get the atoms that define the center of mass of the protein as a list of ints

        Returns
        -------
        list
            A list of atom indices that define the center of mass of the protein.
        """
        return self._mol1_com_atoms

    def _setProtComAtoms(self, prot_com_atoms):
        """
        Set the atoms that define the center of mass of the protein
        If a list is given, simply set them according to the list.
        If None, find them based on the center of mass of the protein.
        """
        if prot_com_atoms is not None:
            # Make sure its a list of ints
            if not isinstance(prot_com_atoms, list):
                raise TypeError("mol1_com_atoms must be a list")
            if not all(isinstance(x, int) for x in prot_com_atoms):
                raise TypeError("mol1_com_atoms must be a list of ints")
            self._mol1_com_atoms = prot_com_atoms
        else:
            # Find com of the protein
            if self._is_prepared:
                temp_system = self._system._sire_object
                protein = temp_system[self.protein_index[0]]
                for i in self.protein_index[1:]:
                    protein += temp_system[i]
                com = protein.coordinates()
                self._mol1_com_atoms = [
                    a.index().value()
                    for a in protein[f"atoms within 11 angstrom of {com}"]
                ]
                del temp_system
                del protein
            else:
                protein = self._protein
                com = protein._sire_object.coordinates()
                self._mol1_com_atoms = [
                    a.index().value()
                    for a in protein._sire_object[f"atoms within 11 angstrom of {com}"]
                ]

    def _getLig1ComAtoms(self):
        """
        Get the atoms that define the center of mass of the bound ligand as a list of ints

        Returns
        -------
        list
            A list of atom indices that define the center of mass of the bound ligand.
        """
        return self._lig1_com_atoms

    def _setLig1ComAtoms(self, lig1_com_atoms):
        """
        Set the atoms that define the center of mass of the bound ligand
        If a list is given, simply set them according to the list.
        If None, find them based on the center of mass of the bound ligand.
        In most cases this will be all atoms within the ligand
        """
        if lig1_com_atoms is not None:
            # Make sure its a list of ints
            if not isinstance(lig1_com_atoms, list):
                raise TypeError("lig1_com_atoms must be a list")
            if not all(isinstance(x, int) for x in lig1_com_atoms):
                raise TypeError("lig1_com_atoms must be a list of ints")
            self._lig1_com_atoms = lig1_com_atoms
        else:
            # Find com of the ligand
            if self._is_prepared:
                ligand_bound = self._system[self.ligand_bound_index]
            else:
                ligand_bound = self._ligand_bound
            com = ligand_bound._sire_object.coordinates()
            self._lig1_com_atoms = [
                a.index().value()
                for a in ligand_bound._sire_object[f"atoms within 11 angstrom of {com}"]
            ]

    def _getLig2ComAtoms(self):
        """
        Get the atoms that define the center of mass of the free ligand as a list of ints

        Returns
        -------
        list
            A list of atom indices that define the center of mass of the free ligand.
        """
        return self._lig2_com_atoms

    def _setLig2ComAtoms(self, lig2_com_atoms):
        """
        Set the atoms that define the center of mass of the free ligand
        If a list is given, simply set them according to the list.
        If None, find them based on the center of mass of the free ligand.
        In most cases this will be all atoms within the ligand
        """
        if lig2_com_atoms is not None:
            # Make sure its a list of ints
            if not isinstance(lig2_com_atoms, list):
                raise TypeError("lig2_com_atoms must be a list")
            if not all(isinstance(x, int) for x in lig2_com_atoms):
                raise TypeError("lig2_com_atoms must be a list of ints")
            self._lig2_com_atoms = lig2_com_atoms
        else:
            # Find com of the ligand
            if self._is_prepared:
                ligand_free = self._system[self.ligand_free_index]
            else:
                ligand_free = self._ligand_free
            com = ligand_free._sire_object.coordinates()
            self._lig2_com_atoms = [
                a.index().value()
                for a in ligand_free._sire_object[f"atoms within 11 angstrom of {com}"]
            ]

    def _makeData(self):
        """
        Make the data dictionary for the AToM system
        """
        self.data = {}
        self.data["displacement"] = self._getDisplacement()
        self.data["protein_index"] = self._getProteinIndex()
        self.data["ligand_bound_index"] = self._getLigandBoundIndex()
        self.data["ligand_free_index"] = self._getLigandFreeIndex()
        self.data["ligand_bound_rigid_core"] = self._getLigandBoundRigidCore()
        self.data["ligand_free_rigid_core"] = self._getLigandFreeRigidCore()
        self.data["mol1_atomcount"] = self._protein_atomcount
        self.data["ligand_bound_atomcount"] = self.ligand_bound_atomcount
        self.data["ligand_free_atomcount"] = self.ligand_free_atomcount
        self.data["first_protein_atom_index"] = self.first_protein_atom_index
        self.data["last_protein_atom_index"] = self.last_protein_atom_index
        self.data["first_ligand_bound_atom_index"] = self.first_ligand_bound_atom_index
        self.data["last_ligand_bound_atom_index"] = self.last_ligand_bound_atom_index
        self.data["first_ligand_free_atom_index"] = self.first_ligand_free_atom_index
        self.data["last_ligand_free_atom_index"] = self.last_ligand_free_atom_index
        self.data["protein_com_atoms"] = self._mol1_com_atoms
        self.data["ligand_bound_com_atoms"] = self._lig1_com_atoms
        self.data["ligand_free_com_atoms"] = self._lig2_com_atoms

    @staticmethod
    def viewRigidCores(
        system=None,
        ligand_bound=None,
        ligand_free=None,
        ligand_bound_rigid_core=None,
        ligand_free_rigid_core=None,
    ):
        """View the rigid cores of the ligands.
        Carbon atoms of the bound ligand are coloured green, while those of the free ligand are coloured orange.
        Rigid core atoms are colour coded and labelled.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The system for the AToM simulation that has been prepared AToM.prepare().
            All other arguments are ignored if this is provided.

        ligand_bound : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            The bound ligand.

        ligand_free : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            The free ligand.

        ligand_bound_rigid_core : list
            The indices for the rigid core atoms of the bound ligand.

        ligand_free_rigid_core : list
            The indices for the rigid core atoms of the free ligand.
        """
        import math as _math

        def move_to_origin(lig):
            com = _Coordinate(*lig._getCenterOfMass())
            lig.translate([-com.x().value(), -com.y().value(), -com.z().value()])

        def euclidean_distance(point1, point2):
            return _math.sqrt(
                (point1[0] - point2[0]) ** 2
                + (point1[1] - point2[1]) ** 2
                + (point1[2] - point2[2]) ** 2
            )

        def furthest_points(points):
            max_distance = 0
            furthest_pair = None
            n = len(points)

            if n < 2:
                return None, None, 0  # Not enough points to compare

            for i in range(n):
                for j in range(i + 1, n):
                    distance = euclidean_distance(points[i], points[j])
                    if distance > max_distance:
                        max_distance = distance
                        furthest_pair = (points[i], points[j])

            return furthest_pair[0], furthest_pair[1], max_distance

        def vector_from_points(point1, point2):
            dx = point2[0] - point1[0]
            dy = point2[1] - point1[1]
            dz = point2[2] - point1[2]

            magnitude = _math.sqrt(dx**2 + dy**2 + dz**2)
            if magnitude == 0:
                return (0, 0, 0)

            return (dx / magnitude, dy / magnitude, dz / magnitude)

        def find_carbons(ligand):
            # Adding logic to find carbon atosm in sire becuase its easier than trying
            # to find them in nglview
            return [
                i.index().value()
                for i in ligand._sire_object[f"(atom mass > 12) and (atom mass < 13)"]
            ]

        def make_amber_list(l, offset=0):
            # take a list of ints and return a string that is "@<int1>,<int2>,<int3>"
            return "@" + ",".join([str(i + offset) for i in l])

        def _count_num_atoms(ligand):
            return ligand._sire_object.nAtoms()

        # if a system is provided, check that it has the "atom_data" property
        if system is not None:
            sdata = _json.loads(system._sire_object.property("atom_data").value())
            local_s = system.copy()
            ligand_bound = local_s[sdata["ligand_bound_index"]]
            move_to_origin(ligand_bound)
            ligand_free = local_s[sdata["ligand_free_index"]]
            move_to_origin(ligand_free)
            ligand_bound_rigid_core = sdata["ligand_bound_rigid_core"]
            ligand_free_rigid_core = sdata["ligand_free_rigid_core"]

        # if not system provided, ALL other parameters must be provided
        else:
            if ligand_bound is None:
                raise ValueError("ligand_bound must be provided")
            if ligand_free is None:
                raise ValueError("ligand_free must be provided")
            if ligand_bound_rigid_core is None:
                raise ValueError("ligand_bound_rigid_core must be provided")
            if ligand_free_rigid_core is None:
                raise ValueError("ligand_free_rigid_core must be provided")

            if not isinstance(ligand_bound, _Molecule):
                raise TypeError("ligand_bound must be a BioSimSpace molecule")
            if not isinstance(ligand_free, _Molecule):
                raise TypeError("ligand_free must be a BioSimSpace molecule")
            if not isinstance(ligand_bound_rigid_core, list):
                raise TypeError("ligand_bound_rigid_core must be a list")
            elif not len(ligand_bound_rigid_core) == 3:
                raise ValueError("ligand_bound_rigid_core must have length 3")
            if not isinstance(ligand_free_rigid_core, list):
                raise TypeError("ligand_free_rigid_core must be a list")
            elif not len(ligand_free_rigid_core) == 3:
                raise ValueError("ligand_free_rigid_core must have length 3")

            # copy the ligands
            ligand_bound = ligand_bound.copy()
            move_to_origin(ligand_bound)
            ligand_free = ligand_free.copy()
            move_to_origin(ligand_free)

        pre_translation_lig1_core_coords = []

        for i in ligand_bound_rigid_core:
            x = ligand_bound.getAtoms()[i].coordinates().x().value()
            y = ligand_bound.getAtoms()[i].coordinates().y().value()
            z = ligand_bound.getAtoms()[i].coordinates().z().value()
            pre_translation_lig1_core_coords.append((x, y, z))

        point1, point2, distance = furthest_points(pre_translation_lig1_core_coords)
        vector = vector_from_points(point1, point2)

        # need to know the size of ligand_bound
        lig1_coords = []
        for i in ligand_bound.getAtoms():
            x = i.coordinates().x().value()
            y = i.coordinates().y().value()
            z = i.coordinates().z().value()
            lig1_coords.append((x, y, z))

        lig1_point1, lig1_point2, lig1_distance = furthest_points(lig1_coords)

        # Translate ligand_free so they don't overlap
        ligand_free.translate(
            [
                -1.0 * lig1_distance * 2 * vector[0],
                -1.0 * lig1_distance * 2 * vector[1],
                -1.0 * lig1_distance * 2 * vector[2],
            ]
        )
        # Get coords of rigid core atoms
        ligand_bound_core_coords = []
        ligand_free_core_coords = []
        for i in ligand_bound_rigid_core:
            ligand_bound_core_coords.append(ligand_bound.getAtoms()[i].coordinates())
        for i in ligand_free_rigid_core:
            ligand_free_core_coords.append(ligand_free.getAtoms()[i].coordinates())

        # Create molecule containing both ligands
        mol = ligand_bound + ligand_free

        # Create view
        view = _ViewAtoM(mol)

        # Create nglview object
        ngl = view.system(mol)

        ngl.add_ball_and_stick("all")
        ngl.add_ball_and_stick(
            make_amber_list(find_carbons(ligand_bound)), color="green"
        )
        ngl.add_ball_and_stick(
            make_amber_list(
                find_carbons(ligand_free), offset=_count_num_atoms(ligand_bound)
            ),
            color="orange",
        )
        #
        colours = [[1, 1, 0], [1, 0, 1], [0, 1, 1]]
        # Add spheres to rigid core locations
        for coord1, coord2, colour, core_atom_1, core_atom_2 in zip(
            ligand_bound_core_coords,
            ligand_free_core_coords,
            colours,
            ligand_bound_rigid_core,
            ligand_free_rigid_core,
        ):
            ngl.shape.add_sphere(
                [coord1.x().value(), coord1.y().value(), coord1.z().value()],
                colour,
                0.45,
            )
            ngl.shape.add(
                "text",
                [coord1.x().value(), coord1.y().value(), coord1.z().value() - 0.5],
                [0, 0, 0],
                2.5,
                f"{core_atom_1}",
            )
            ngl.shape.add_sphere(
                [coord2.x().value(), coord2.y().value(), coord2.z().value()],
                colour,
                0.45,
            )
            ngl.shape.add(
                "text",
                [coord2.x().value(), coord2.y().value(), coord2.z().value() - 0.5],
                [0, 0, 0],
                2.5,
                f"{core_atom_2}",
            )

        if system is not None:
            del local_s
        return ngl

    @staticmethod
    def run(
        system,
        protocol,
        platform="CPU",
        work_dir=None,
        setup_only=False,
        property_map={},
    ):
        """Run the AToM production simulation(s).

        Parameters
        ----------
        system : :class:`System <BioSimSpace._SireWrappers.System>`
            A prepared AToM system.

        protocol : :class:`Protocol <BioSimSpace.Protocol.AToMProduction>`
            A protocol object that defines the AToM protocol.

        platform : str
            The platform for the simulation.

        work_dir : str
            The working directory for the simulation.

        setup_only : bool
            Whether to only support simulation setup. If True, then no
            simulation processes objects will be created, only the directory
            hierarchy and input files to run a simulation externally. This
            can be useful when you don't intend to use BioSimSpace to run
            the simulation. Note that a 'work_dir' must also be specified.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        """
        runner = _relativeATM(
            system,
            protocol,
            platform,
            property_map=property_map,
            work_dir=work_dir,
            setup_only=setup_only,
        )
        runner.run()
        if setup_only:
            return None
        else:
            return runner._runner

    @staticmethod
    def analyse(
        work_dir,
        method="UWHAM",
        ignore_lower=0,
        ignore_upper=None,
        inflex_indices=None,
    ):
        """Analyse the AToM simulation.

        Parameters
        ----------

        work_dir : str
            The working directory where the AToM simulation is located.

        method : str
            The method to use for the analysis. Currently only UWHAM is supported.

        ignore_lower : int
            Ignore the first N samples when analysing.

        inflex_indices : [int]
            The indices at which the direction changes. For example, if direction=[1,1,-1,-1],
            then inflex_indices=[1,2].
            If None, the inflexion point will be found automatically.

        Returns
        -------

        ddg : :class:`BioSimSpace.Types.Energy`
            The free energy difference between the two ligands.

        ddg_err :class:`BioSimSpace.Types.Energy`
            The error in the free energy difference.
        """
        if not isinstance(ignore_lower, int):
            raise TypeError("'ignore_lower' must be an integer.")
        if ignore_lower < 0:
            raise ValueError("'ignore_lower' must be a positive integer.")
        if ignore_upper is not None:
            if not isinstance(ignore_upper, int):
                raise TypeError("'ignore_upper' must be an integer.")
            if ignore_upper < 0:
                raise ValueError("'ignore_upper' must be a positive integer.")
            if ignore_upper < ignore_lower:
                raise ValueError(
                    "'ignore_upper' must be greater than or equal to 'ignore_lower'."
                )
        if inflex_indices is not None:
            if not isinstance(inflex_indices, list):
                raise TypeError("'inflex_indices' must be a list.")
            if not all(isinstance(x, int) for x in inflex_indices):
                raise TypeError("'inflex_indices' must be a list of integers.")
            if not len(inflex_indices) == 2:
                raise ValueError("'inflex_indices' must have length 2.")
        if method == "UWHAM":
            total_ddg, total_ddg_err = AToM._analyse_UWHAM(
                work_dir, ignore_lower, ignore_upper, inflex_indices
            )
            return total_ddg, total_ddg_err
        if method == "MBAR":
            from ._relative import Relative as _Relative

            # temporary version to check that things are working
            ddg_forward, ddg_reverse = AToM._analyse_MBAR(work_dir)
            ddg_forward = _Relative.difference(ddg_forward)
            ddg_reverse = _Relative.difference(ddg_reverse)
            return ddg_forward, ddg_reverse
        else:
            raise ValueError(f"Method {method} is not supported for analysis.")

    @staticmethod
    def _analyse_UWHAM(work_dir, ignore_lower, ignore_upper, inflex_indices=None):
        """
        Analyse the UWHAM results from the AToM simulation.
        """
        from ._ddg import analyse_UWHAM as _UWHAM

        total_ddg, total_ddg_err = _UWHAM(
            work_dir, ignore_lower, ignore_upper, inflection_indices=inflex_indices
        )
        return total_ddg, total_ddg_err

    @staticmethod
    def _analyse_MBAR(work_dir):
        """
        Analyse the MBAR results from the AToM simulation.
        """
        from ._ddg import analyse_MBAR as _MBAR

        ddg_forward, ddg_reverse = _MBAR(work_dir)
        return ddg_forward, ddg_reverse

    @staticmethod
    def _analyse_test(work_dir):
        """
        Analyse the test results from the AToM simulation.
        """
        from ._ddg import new_MBAR as _test

        ddg_forward, ddg_reverse = _test(work_dir)
        return ddg_forward, ddg_reverse

    @staticmethod
    def _analyse_femto(work_dir):
        from ._ddg import MBAR_hijack_femto

        est, o = MBAR_hijack_femto(work_dir)
        return est, o


class _relativeATM:
    """
    A class for setting up and performing RBFE calculations using AToM
    """

    def __init__(
        self,
        system,
        protocol,
        platform="CPU",
        work_dir=None,
        setup_only=False,
        property_map={},
    ):
        """
        Constructor.

        Parameters
        ----------

        system : BioSimSpace._SireWrappers.System
            A prepared AToM system containing a protein and two ligands, one bound and one free.
            Assumed to already be equilibrated.

        protocol : BioSimSpace.Protocol.AToM
            A protocol object that defines the RBFE protocol.

        platform : str
            The platform for the simulation: “CPU”, “CUDA”, or “OPENCL”.
            For CUDA use the CUDA_VISIBLE_DEVICES environment variable to set the GPUs on which to run,
            e.g. to run on two GPUs indexed 0 and 1 use: CUDA_VISIBLE_DEVICES=0,1.
            For OPENCL, instead use OPENCL_VISIBLE_DEVICES.

        work_dir : str
            The working directory for the simulation.

        setup_only : bool
            Whether to only support simulation setup. If True, then no
            simulation processes objects will be created, only the directory
            hierarchy and input files to run a simulation externally. This
            can be useful when you don't intend to use BioSimSpace to run
            the simulation. Note that a 'work_dir' must also be specified.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        """

        self._system = system.copy()

        # Validate the protocol.
        if protocol is not None:
            from ..Protocol._atm import AToMProduction as _Production

            if not isinstance(protocol, _Production):
                raise TypeError(
                    "'protocol' must be of type 'BioSimSpace.Protocol.AToMProduction'"
                )
            else:
                self._protocol = protocol
        else:
            # No default protocol due to the need for well-defined rigid cores
            raise ValueError("A protocol must be specified")

        # Check the platform.
        if not isinstance(platform, str):
            raise TypeError("'platform' must be of type 'str'.")
        else:
            self._platform = platform

        if not isinstance(setup_only, bool):
            raise TypeError("'setup_only' must be of type 'bool'.")
        else:
            self._setup_only = setup_only

        if work_dir is None and setup_only:
            raise ValueError(
                "A 'work_dir' must be specified when 'setup_only' is True!"
            )

        # Create the working directory.
        self._work_dir = _Utils.WorkDir(work_dir)

        # Check that the map is valid.
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")
        self._property_map = property_map

        self._inititalise_runner(system=self._system)

    def run(self, serial=True):
        """
        Run the simulations.

        Returns
        -------
        Processes : [:class:`Process <BioSimSpace.Process>`]
            A list of process objects.
        """
        # Initialise the runner.
        if not isinstance(serial, bool):
            raise TypeError("'serial' must be of type 'bool'.")

        if self._setup_only:
            _warnings.warn("No processes exist! Object created in 'setup_only' mode.")

        else:
            self._runner.startAll(serial=serial)

    def _inititalise_runner(self, system):
        """
        Internal helper function to initialise the process runner.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.
        """

        # This protocol will have to be minimal - cannot guess rigid core atoms
        if self._protocol is None:
            raise RuntimeError("No protocol has been set - cannot run simulations.")
        # Initialise list to store the processe
        processes = []
        # Get the list of lambda1 values so that the total number of simulations can
        # be asserted
        lambda_list = self._protocol._get_lambda_values()
        # Set index of current simulation to 0
        self._protocol.set_current_index(0)
        lam = lambda_list[0]

        first_dir = "%s/lambda_%5.4f" % (self._work_dir, lam)

        # Create the first simulation, which will be copied and used for future simulations.
        first_process = _OpenMM(
            system=system,
            protocol=self._protocol,
            platform=self._platform,
            work_dir=first_dir,
            property_map=self._property_map,
        )

        if self._setup_only:
            del first_process
        else:
            processes.append(first_process)

        # Remove first index as its already been used
        lambda_list = lambda_list[1:]
        # Enumerate starting at 1 to account for the removal of the first lambda value
        for index, lam in enumerate(lambda_list, 1):
            # Files are named according to index, rather than lambda value
            # This is to avoid confusion arising from the fact that there are multiple lambdas
            # and that the values of lambda1 and lambda2 wont necessarily be go from 0 to 1
            # and may contain duplicates
            new_dir = "%s/lambda_%5.4f" % (self._work_dir, lam)
            # Use absolute path.
            if not _os.path.isabs(new_dir):
                new_dir = _os.path.abspath(new_dir)

            # Delete any existing directories.
            if _os.path.isdir(new_dir):
                _shutil.rmtree(new_dir, ignore_errors=True)

            # Copy the first directory to that of the current lambda value.
            _shutil.copytree(first_dir, new_dir)
            # For speed reasons, additional processes need to be created by copying the first process.
            # this is more difficult than usual due to the number of window-dependent variables
            new_config = []
            # All variables that need to change
            new_lam_1 = self._protocol.getLambda1()[index]
            new_lam_2 = self._protocol.getLambda2()[index]
            new_alpha = self._protocol.getAlpha()[index].value()
            new_uh = self._protocol.getUh()[index].value()
            new_w0 = self._protocol.getW0()[index].value()
            new_direction = self._protocol.getDirection()[index]
            with open(new_dir + "/openmm_script.py", "r") as f:
                for line in f:
                    if line.startswith("lambda1"):
                        new_config.append(f"lambda1 = {new_lam_1}\n")
                    elif line.startswith("lambda2"):
                        new_config.append(f"lambda2 = {new_lam_2}\n")
                    elif line.startswith("alpha"):
                        new_config.append(
                            f"alpha = {new_alpha} * kilocalories_per_mole\n"
                        )
                    elif line.startswith("uh"):
                        new_config.append(f"uh = {new_uh} * kilocalories_per_mole\n")
                    elif line.startswith("w0"):
                        new_config.append(f"w0 = {new_w0} * kilocalories_per_mole\n")
                    elif line.startswith("direction"):
                        new_config.append(f"direction = {new_direction}\n")
                    elif line.startswith("window_index"):
                        new_config.append(f"window_index = {index}\n")
                    else:
                        new_config.append(line)
            with open(new_dir + "/openmm_script.py", "w") as f:
                for line in new_config:
                    f.write(line)

            # biosimspace runner functionality
            if not self._setup_only:
                process = _copy.copy(first_process)
                process._system = first_process._system.copy()
                process._protocol = self._protocol
                process._work_dir = new_dir
                process._stdout_file = new_dir + "/AToM.out"
                process._stderr_file = new_dir + "/AToM.err"
                process._rst_file = new_dir + "/openmm.rst7"
                process._top_file = new_dir + "/openmm.prm7"
                process._traj_file = new_dir + "/openmm.dcd"
                process._config_file = new_dir + "/openmm_script.py"
                process._input_files = [
                    process._config_file,
                    process._rst_file,
                    process._top_file,
                ]
                processes.append(process)

        if not self._setup_only:
            # Initialise process runner.
            self._runner = _ProcessRunner(processes)


class _ViewAtoM(_View):
    """
    Overloads regular view class needed to pass default_representation=False
    into show_file.
    """

    # Initialise super
    def __init__(self, handle, property_map={}, is_lambda1=False):
        super().__init__(handle, property_map, is_lambda1)

    def _create_view(self, system=None, view=None, gui=True, **kwargs):

        if system is None and view is None:
            raise ValueError("Both 'system' and 'view' cannot be 'None'.")

        elif system is not None and view is not None:
            raise ValueError("One of 'system' or 'view' must be 'None'.")

        # Make sure gui flag is valid.
        if gui not in [True, False]:
            gui = True

        # Default to the most recent view.
        if view is None:
            index = self._num_views
        else:
            index = view

        # Create the file name.
        filename = "%s/view_%04d.pdb" % (self._work_dir, index)

        # Increment the number of views.
        if view is None:
            self._num_views += 1

        # Create a PDB object and write to file.
        if system is not None:
            try:
                pdb = _SireIO.PDB2(system, self._property_map)
                pdb.writeToFile(filename)
            except Exception as e:
                msg = "Failed to write system to 'PDB' format."
                if _isVerbose():
                    print(msg)
                    raise IOError(e) from None
                else:
                    raise IOError(msg) from None

        # Import NGLView when it is used for the first time.
        import nglview as _nglview

        # Create the NGLview object.
        view = _nglview.show_file(filename, default_representation=False)

        # Return the view and display it.
        return view.display(gui=gui)
