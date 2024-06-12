######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2023
#
# Authors: Lester Hedges <lester.hedges@gmail.com>, Matthew Burman <matthew@openbiosim.org>
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


from .. import _is_notebook
from .._SireWrappers import Molecule as _Molecule
from .._SireWrappers import System as _System
from .. import _Utils
from ..Types import Length as _Length
from ..Types import Vector as _Vector
from ..Types import Coordinate as _Coordinate
from ..Align import matchAtoms as _matchAtoms
from ..Align import rmsdAlign as _rmsdAlign
from ..Notebook import View as _View
from ..Process._atom_utils import _AToMUtils as _AtomUtils
import warnings as _warnings
import json as _json
import copy as _copy


class AToM:
    """A class for creating, setting up, running and analysing AToM
    simulations."""

    def __init__(
        self,
        system=None,
        protein=None,
        ligand1=None,
        ligand2=None,
        protein_index=0,
        ligand1_index=1,
        ligand2_index=2,
    ):
        """Constructor for the AToM class.

        Parameters
        ----------
        system : BioSimSpace._SireWrappers.System
            A pre-prepared AToM system containing protein and ligands placed in their correct positions.
            If provided takes precedence over protein, ligand1 and ligand2.
        protein : BioSimSpace._SireWrappers.Molecule
            A protein molecule. Will be used along with ligand1 and ligand2 to create a system.
        ligand1 : BioSimSpace._SireWrappers.Molecule
            The bound ligand. Will be used along with protein and ligand2 to create a system.
        ligand2 : BioSimSpace._SireWrappers.Molecule
            The free ligand. Will be used along with protein and ligand1 to create a system.
        protein_index: int, list
            If passing a pre-prepared system, the index (or indices) of the protein molecule in the system (Default 0).
        ligand1_index: int
            If passing a pre-prepared system, the index of the bound ligand molecule in the system (Default 1).
        ligand2_index: int
            If passing a pre-prepared system, the index of the free ligand molecule in the system (Default 2).
        """
        # make sure that either system or protein, ligand1 and ligand2 are given
        if system is None and not all(
            x is not None for x in [protein, ligand1, ligand2]
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
            for x in [protein, ligand1, ligand2]
            if x is not None
        ):
            raise ValueError(
                "The protein, bound ligand and free ligand must be BioSimSpace Molecule objects."
            )
        self._is_prepared = False
        self._setSystem(system)
        if not self._is_prepared:
            self._setProtein(protein)
            self._setLigand1(ligand1)
            self._setLigand2(ligand2)
        else:
            self._setProteinIndex(protein_index)
            self._setLigand1Index(ligand1_index)
            self._setLigand2Index(ligand2_index)

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

    def _setLigand1(self, ligand1):
        """Set the bound ligand for the AToM simulation.

        Parameters
        ----------
        ligand1 : BioSimSpace._SireWrappers.Molecule
            The bound ligand for the AToM simulation.
        """
        if ligand1 is not None:
            if not isinstance(ligand1, _Molecule):
                raise ValueError(
                    "The bound ligand must be a BioSimSpace Molecule object."
                )
            else:
                self._ligand1 = ligand1
        else:
            self._ligand1 = None

    def _getLigand1(self):
        """Get the bound ligand for the AToM simulation.

        Returns
        -------
        BioSimSpace._SireWrappers.Molecule
            The bound ligand for the AToM simulation.
        """
        return self._ligand1

    def _setLigand2(self, ligand2):
        """Set the free ligand for the AToM simulation.

        Parameters
        ----------
        ligand2 : BioSimSpace._SireWrappers.Molecule
            The free ligand for the AToM simulation.
        """
        if ligand2 is not None:
            if not isinstance(ligand2, _Molecule):
                raise ValueError(
                    "The free ligand must be a BioSimSpace Molecule object."
                )
            else:
                self._ligand2 = ligand2
        else:
            self._ligand2 = None

    def _getLigand2(self):
        """Get the free ligand for the AToM simulation.

        Returns
        -------
        BioSimSpace._SireWrappers.Molecule
            The free ligand for the AToM simulation.
        """
        return self._ligand2

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

    def _setLigand1RigidCore(self, ligand1_rigid_core):
        """Set the indices for the rigid core atoms of ligand 1.

        Parameters
        ----------
        ligand1_rigid_core : BioSimSpace._SireWrappers.Molecule
            The rigid core of the bound ligand for the AToM simulation.
        """
        if ligand1_rigid_core is None:
            self.ligand1_rigid_core = None
        else:
            if not isinstance(ligand1_rigid_core, list):
                raise TypeError("ligand1_rigid_core must be a list")
            if len(ligand1_rigid_core) != 3:
                raise ValueError("ligand1_rigid_core must have length 3")
            # make sure all indices are ints
            if not all(isinstance(x, int) for x in ligand1_rigid_core):
                raise TypeError("ligand1_rigid_core must contain only integers")
            if any(x >= self.ligand1_atomcount for x in ligand1_rigid_core):
                raise ValueError(
                    "ligand1_rigid_core contains an index that is greater than the number of atoms in the ligand"
                )
            self.ligand1_rigid_core = ligand1_rigid_core

    def _getLigand1RigidCore(self):
        """Get the indices for the rigid core atoms of ligand 1.

        Returns
        -------
        list
            The indices for the rigid core atoms of ligand 1.
        """
        return self.ligand1_rigid_core

    def _setLigand2RigidCore(self, ligand2_rigid_core):
        """Set the indices for the rigid core atoms of ligand 2.

        Parameters
        ----------
        ligand2_rigid_core : BioSimSpace._SireWrappers.Molecule
            The rigid core of the free ligand for the AToM simulation.
        """
        if ligand2_rigid_core is None:
            self.ligand2_rigid_core = None
        else:
            if not isinstance(ligand2_rigid_core, list):
                raise TypeError("ligand2_rigid_core must be a list")
            if len(ligand2_rigid_core) != 3:
                raise ValueError("ligand2_rigid_core must have length 3")
            # make sure all indices are ints
            if not all(isinstance(x, int) for x in ligand2_rigid_core):
                raise TypeError("ligand2_rigid_core must contain only integers")
            if any(x >= self.ligand2_atomcount for x in ligand2_rigid_core):
                raise ValueError(
                    "ligand2_rigid_core contains an index that is greater than the number of atoms in the ligand"
                )
            self.ligand2_rigid_core = ligand2_rigid_core

    def _getLigand2RigidCore(self):
        """Get the indices for the rigid core atoms of ligand 2.

        Returns
        -------
        list
            The indices for the rigid core atoms of ligand 2.
        """
        return self.ligand2_rigid_core

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

    def _setLigand1Index(self, ligand1_index):
        """Set the index of the bound ligand molecule in the system.

        Parameters
        ----------
        ligand1_index : int
            The index of the bound ligand molecule in the system.
        """
        if not isinstance(ligand1_index, int):
            raise ValueError("ligand1_index must be an integer.")
        else:
            if ligand1_index < 0:
                raise ValueError("ligand1_index must be a positive integer")
            if self._system[ligand1_index].isWater():
                _warnings.warn(
                    f"The molecule at index {ligand1_index} is a water molecule, check your ligand1_index."
                )
            self.ligand1_index = ligand1_index

    def _getLigand1Index(self):
        """Get the index of the bound ligand molecule in the system.

        Returns
        -------
        int
            The index of the bound ligand molecule in the system.
        """
        return self.ligand1_index

    def _setLigand2Index(self, ligand2_index):
        """Set the index of the free ligand molecule in the system.

        Parameters
        ----------
        ligand2_index : int
            The index of the free ligand molecule in the system.
        """
        if not isinstance(ligand2_index, int):
            raise ValueError("ligand2_index must be an integer.")
        else:
            if ligand2_index < 0:
                raise ValueError("ligand2_index must be a positive integer")
            if self._system[ligand2_index].isWater():
                _warnings.warn(
                    f"The molecule at index {ligand2_index} is a water molecule, check your ligand2_index."
                )
            self.ligand2_index = ligand2_index

    def _getLigand2Index(self):
        """Get the index of the free ligand molecule in the system.

        Returns
        -------
        int
            The index of the free ligand molecule in the system.
        """
        return self.ligand2_index

    def prepare(
        self,
        ligand1_rigid_core,
        ligand2_rigid_core,
        displacement="20A",
        protein_com_atoms=None,
        ligand1_com_atoms=None,
        ligand2_com_atoms=None,
    ):
        """Prepare the system for the AToM simulation.

        Parameters
        ----------
        ligand1_rigid_core : list
            A list of three atom indices that define the rigid core of the bound ligand.
            Indices are set relative to the ligand, not the system and are 0-indexed.
        ligand2_rigid_core : list
            A list of three atom indices that define the rigid core of the free ligand.
            Indices are set relative to the ligand, not the system and are 0-indexed.
        displacement : float, string, list
            The diplacement between the bound and free ligands.
            If a float or string is given, BioSimSpace will attempt to find the ideal vector along which to displace the ligand by the given magnitude.
            If a list is given, the vector will be used directly.
            Lengths should always be given in angstroms.
            Default is 20A.
        protein_com_atoms : list
            A list of atom indices that define the center of mass of the protein.
            If None, the center of mass of the protein will be found automatically.
        ligand1_com_atoms : list
            A list of atom indices that define the center of mass of the bound ligand.
            If None, the center of mass of the bound ligand will be found automatically.
        ligand2_com_atoms : list
            A list of atom indices that define the center of mass of the free ligand.
            If None, the center of mass of the free ligand will be found automatically.
        """
        if self._is_prepared:
            self._systemInfo()
            self._setLigand1RigidCore(ligand1_rigid_core)
            self._setLigand2RigidCore(ligand2_rigid_core)
            self._setDisplacement(displacement)
            self._setProtComAtoms(protein_com_atoms)
            self._setLig1ComAtoms(ligand1_com_atoms)
            self._setLig2ComAtoms(ligand2_com_atoms)

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
                self._protein, self._ligand1, self._ligand2, self.displacement
            )
            self._setSystem(system, is_prepped=False)
            self._setDisplacement(dis_vec)
            self._setProteinIndex(prot_ind)
            self._setLigand1Index(lig1_ind)
            self._setLigand2Index(lig2_ind)
            self._systemInfo()
            self._setLigand1RigidCore(ligand1_rigid_core)
            self._setLigand2RigidCore(ligand2_rigid_core)
            self._setProtComAtoms(protein_com_atoms)
            self._setLig1ComAtoms(ligand1_com_atoms)
            self._setLig2ComAtoms(ligand2_com_atoms)
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
    def _makeSystemFromThree(protein, ligand1, ligand2, displacement):
        """Create a system for AToM simulations.

        Parameters
        ----------
        protein : BioSimSpace._SireWrappers.Molecule
            The protein for the AToM simulation.
        ligand1 : BioSimSpace._SireWrappers.Molecule
            The bound ligand for the AToM simulation.
        ligand2 : BioSimSpace._SireWrappers.Molecule
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

        mapping = _matchAtoms(ligand2, ligand1)
        ligand2_aligned = _rmsdAlign(ligand2, ligand1, mapping)
        prot_lig1 = (protein + ligand1).toSystem()

        if isinstance(displacement, _Vector):
            ligand2_aligned.translate(
                [displacement.x(), displacement.y(), displacement.z()]
            )
            vec = displacement
        else:
            vec = _findTranslationVector(prot_lig1, displacement, protein, ligand1)
            ligand2_aligned.translate([vec.x(), vec.y(), vec.z()])

        system = (protein + ligand1 + ligand2_aligned).toSystem()
        prot_ind = system.getIndex(protein)
        lig1_ind = system.getIndex(ligand1)
        lig2_ind = system.getIndex(ligand2_aligned)
        return system, prot_ind, lig1_ind, lig2_ind, vec

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
        if self._system[self.ligand1_index].isWater():
            _warnings.warn(
                f"The molecule at index {self.ligand1_index} appears to be a water molecule."
                " This should be the bound ligand."
            )
        if self._system[self.ligand2_index].isWater():
            _warnings.warn(
                f"The molecule at index {self.ligand2_index} appears to be a water molecule."
                " This should be the free ligand."
            )
        self._protein_atomcount = sum(
            self._system[i].nAtoms() for i in self.protein_index
        )
        self.ligand1_atomcount = self._system[self.ligand1_index].nAtoms()
        self.ligand2_atomcount = self._system[self.ligand2_index].nAtoms()

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

        ligand1_atom_start = self._system[self.ligand1_index].getAtoms()[0]
        ligand1_atom_end = self._system[self.ligand1_index].getAtoms()[-1]
        self.first_ligand1_atom_index = self._system.getIndex(ligand1_atom_start)
        self.last_ligand1_atom_index = self._system.getIndex(ligand1_atom_end)

        ligand2_atom_start = self._system[self.ligand2_index].getAtoms()[0]
        ligand2_atom_end = self._system[self.ligand2_index].getAtoms()[-1]
        self.first_ligand2_atom_index = self._system.getIndex(ligand2_atom_start)
        self.last_ligand2_atom_index = self._system.getIndex(ligand2_atom_end)

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
                ligand1 = self._system[self._ligand1_index]
            else:
                ligand1 = self._ligand1
            com = ligand1._sire_object.coordinates()
            self._lig1_com_atoms = [
                a.index().value()
                for a in ligand1._sire_object[f"atoms within 11 angstrom of {com}"]
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
                ligand2 = self._system[self._ligand2_index]
            else:
                ligand2 = self._ligand2
            com = ligand2._sire_object.coordinates()
            self._lig2_com_atoms = [
                a.index().value()
                for a in ligand2._sire_object[f"atoms within 11 angstrom of {com}"]
            ]

    def _makeData(self):
        """
        Make the data dictionary for the AToM system
        """
        self.data = {}
        self.data["displacement"] = self._getDisplacement()
        self.data["protein_index"] = self._getProteinIndex()
        self.data["ligand1_index"] = self._getLigand1Index()
        self.data["ligand2_index"] = self._getLigand2Index()
        self.data["ligand1_rigid_core"] = self._getLigand1RigidCore()
        self.data["ligand2_rigid_core"] = self._getLigand2RigidCore()
        self.data["mol1_atomcount"] = self._protein_atomcount
        self.data["ligand1_atomcount"] = self.ligand1_atomcount
        self.data["ligand2_atomcount"] = self.ligand2_atomcount
        self.data["first_protein_atom_index"] = self.first_protein_atom_index
        self.data["last_protein_atom_index"] = self.last_protein_atom_index
        self.data["first_ligand1_atom_index"] = self.first_ligand1_atom_index
        self.data["last_ligand1_atom_index"] = self.last_ligand1_atom_index
        self.data["first_ligand2_atom_index"] = self.first_ligand2_atom_index
        self.data["last_ligand2_atom_index"] = self.last_ligand2_atom_index
        self.data["protein_com_atoms"] = self._mol1_com_atoms
        self.data["ligand1_com_atoms"] = self._lig1_com_atoms
        self.data["ligand2_com_atoms"] = self._lig2_com_atoms

    @staticmethod
    def viewRigidCores(
        system=None,
        ligand1=None,
        ligand2=None,
        ligand1_rigid_core=None,
        ligand2_rigid_core=None,
    ):
        """View the rigid cores of the ligands.

        Parameters
        ----------
        system : BioSimSpace._SireWrappers.System
            The system for the AToM simulation that has been prepared AToM.prepare().
            All other parameters are ignored if this is provided.
        ligand1 : BioSimSpace._SireWrappers.Molecule
            The bound ligand.
        ligand2 : BioSimSpace._SireWrappers.Molecule
            The free ligand.
        ligand1_rigid_core : list
            The indices for the rigid core atoms of the bound ligand.
        ligand2_rigid_core : list
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
            return (point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2])

        # if a system is provided, check that it has the "atom_data" property
        if system is not None:
            sdata = _json.loads(system._sire_object.property("atom_data").value())
            local_s = system.copy()
            ligand1 = local_s[sdata["ligand1_index"]]
            move_to_origin(ligand1)
            ligand2 = local_s[sdata["ligand2_index"]]
            move_to_origin(ligand2)
            ligand1_rigid_core = sdata["ligand1_rigid_core"]
            ligand2_rigid_core = sdata["ligand2_rigid_core"]

        # if not system provided, ALL other parameters must be provided
        else:
            if ligand1 is None:
                raise ValueError("ligand1 must be provided")
            if ligand2 is None:
                raise ValueError("ligand2 must be provided")
            if ligand1_rigid_core is None:
                raise ValueError("ligand1_rigid_core must be provided")
            if ligand2_rigid_core is None:
                raise ValueError("ligand2_rigid_core must be provided")

            if not isinstance(ligand1, _Molecule):
                raise TypeError("ligand1 must be a BioSimSpace molecule")
            if not isinstance(ligand2, _Molecule):
                raise TypeError("ligand2 must be a BioSimSpace molecule")
            if not isinstance(ligand1_rigid_core, list):
                raise TypeError("ligand1_rigid_core must be a list")
            elif not len(ligand1_rigid_core) == 3:
                raise ValueError("ligand1_rigid_core must have length 3")
            if not isinstance(ligand2_rigid_core, list):
                raise TypeError("ligand2_rigid_core must be a list")
            elif not len(ligand2_rigid_core) == 3:
                raise ValueError("ligand2_rigid_core must have length 3")

            # copy the ligands
            ligand1 = ligand1.copy()
            move_to_origin(ligand1)
            ligand2 = ligand2.copy()
            move_to_origin(ligand2)

        pre_translation_lig1_core_coords = []

        for i in ligand1_rigid_core:
            x = ligand1.getAtoms()[i].coordinates().x().value()
            y = ligand1.getAtoms()[i].coordinates().y().value()
            z = ligand1.getAtoms()[i].coordinates().z().value()
            pre_translation_lig1_core_coords.append((x, y, z))

        point1, point2, distance = furthest_points(pre_translation_lig1_core_coords)
        vector = vector_from_points(point1, point2)

        # need to know the size of ligand1
        lig1_coords = []
        for i in ligand1.getAtoms():
            x = i.coordinates().x().value()
            y = i.coordinates().y().value()
            z = i.coordinates().z().value()
            lig1_coords.append((x, y, z))

        lig1_point1, lig1_point2, lig1_distance = furthest_points(lig1_coords)

        # Translate ligand2 so they don't overlap
        ligand2.translate(
            [
                lig1_distance * vector[0],
                lig1_distance * vector[1],
                lig1_distance * vector[2],
            ]
        )
        # Get coords of rigid core atoms
        ligand1_core_coords = []
        ligand2_core_coords = []
        for i in ligand1_rigid_core:
            ligand1_core_coords.append(ligand1.getAtoms()[i].coordinates())
        for i in ligand2_rigid_core:
            ligand2_core_coords.append(ligand2.getAtoms()[i].coordinates())

        # Create molecule containing both ligands
        mol = ligand1 + ligand2

        # Create view
        view = _View(mol)

        # Create nglview object
        ngl = view.system(mol)
        colours = [[1, 1, 0], [1, 0, 1], [0, 1, 1]]
        # Add spheres to rigid core locations
        for coord1, coord2, colour in zip(
            ligand1_core_coords, ligand2_core_coords, colours
        ):
            ngl.shape.add_sphere(
                [coord1.x().value(), coord1.y().value(), coord1.z().value()],
                colour,
                0.7,
            )
            ngl.shape.add_sphere(
                [coord2.x().value(), coord2.y().value(), coord2.z().value()],
                colour,
                0.7,
            )
        if system is not None:
            del local_s
        return ngl
