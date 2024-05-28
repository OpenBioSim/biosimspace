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

__all__ = [
    "makeSystem",
    "viewRigidCores",
    "relativeATM",
    "Minimise",
    "Equilibrate",
    "Anneal",
]

from ... import _is_notebook
from ..._SireWrappers import Molecule as _Molecule
from ..._SireWrappers import System as _System
from ... import _Utils
from ...Types import Length as _Length
from ...Types import Vector as _Vector
from ...Types import Coordinate as _Coordinate
from ...Align import matchAtoms as _matchAtoms
from ...Align import rmsdAlign as _rmsdAlign
from ... import Process as _Process
from ...Notebook import View as _View
from ...Process._atom_utils import _AToMUtils as _AtomUtils
import os as _os
import shutil as _shutil
import copy as _copy
import warnings as _warnings
import math as _math


class makeSystem:
    """
    A class for creating AToM systems
    """

    def __init__(
        self,
        mol1,
        ligand1=None,
        ligand2=None,
        displacement="20A",
        protein_index=0,
        ligand1_index=1,
        ligand2_index=2,
        ligand1_rigid_core=None,
        ligand2_rigid_core=None,
        prot_com_atoms=None,
        lig1_com_atoms=None,
        lig2_com_atoms=None,
    ):
        """
        Constructor for AToMSystem class

        Parameters
        ----------
        mol1 : BioSimSpace._SireWrappers.Molecule, BioSimSpace._SireWrappers.Molecules, BioSimSpace._SireWrappers.System
            The first molecule/system - this can be one of the following:
            - A single protein (BioSimSpace._SireWrappers.Molecule)
            - A pre-prepared protein-ligand-ligand system ready for use in AToM (BioSimSpace._SireWrappers.System).
        ligand1 : BioSimSpace._SireWrappers.Molecule
            The bound ligand.
        ligand1 : BioSimSpace._SireWrappers.Molecule
            The free ligand.
        displacement : float, list
            The diplacement between the bound and free ligands.
            If a float is given, BioSimSpace will attempt to find the ideal vector along which to displace the ligand by the given magnitude.
            If a list is given, the vector will be used directly.
            Lengths should always be given in angstroms.
        protein_index : int, list of int
            The index of the protein in the system (only needed if passing in a pre-prepared system).
        ligand1_index : int
            The index of the bound ligand in the system (only needed if passing in a pre-prepared system).
        ligand2_index : int
            The index of the free ligand in the system (only needed if passing in a pre-prepared system).
        ligand1_rigid_core : list
            A list of three atom indices that define the rigid core of the bound ligand.
            Indices are set relative to the ligand, not the system and are 0-indexed.
        ligand2_rigid_core : list
            A list of three atom indices that define the rigid core of the free ligand.
            Indices are set relative to the ligand, not the system and are 0-indexed.
        protein_com_atoms: list
            A list of atom indices that define the center of mass of the protein.
            If None, the center of mass will be calculated.
        lig1_com_atoms: list
            A list of atom indices that define the center of mass of the bound ligand.
            If None, the center of mass will be calculated.
        lig2_com_atoms: list
            A list of atom indices that define the center of mass of the free ligand.
            If None, the center of mass will be calculated.
        """
        if isinstance(mol1, _Molecule) and (
            not isinstance(ligand1, _Molecule) or not isinstance(ligand2, _Molecule)
        ):
            raise ValueError(
                "ligand1 and ligand2 must be specified if mol1 is a Molecule"
            )
        self.data = {}
        self._is_pre_prepared = False
        self._is_made = False
        self._setProteinIndex(protein_index)
        self._setLigand1Index(ligand1_index)
        self._setLigand2Index(ligand2_index)
        self._setMol1(mol1)
        self._setLigand1(ligand1)
        self._setLigand2(ligand2)
        self.setLigand1RigidCore(ligand1_rigid_core)
        self.setLigand2RigidCore(ligand2_rigid_core)
        self._setDisplacement(displacement)
        if isinstance(mol1, _Molecule):
            print("Making from three")
            self._makeSystemFromThree()
        else:
            # These will be updated if/when needed
            self._setProteinIndex(protein_index)
            self._setLigand1Index(ligand1_index)
            self._setLigand2Index(ligand2_index)
        self._set_mol1_com_atoms(prot_com_atoms)
        self._set_lig1_com_atoms(lig1_com_atoms)
        self._set_lig2_com_atoms(lig2_com_atoms)

    def getProteinIndex(self):
        """
        Get the index of the protein in the system

        Returns
        -------
        int
            The index of the protein in the system.
        """
        return self.protein_index

    def _setProteinIndex(self, protein_index):
        """
        Set the index of the protein in the system

        Parameters
        ----------
        protein_index : int
            The index of the protein in the system.
        """
        if isinstance(protein_index, list):
            # check that all elements are ints
            if not all(isinstance(x, int) for x in protein_index):
                raise TypeError("protein_index must be a list of ints or a single int")
            self.protein_index = protein_index
        elif isinstance(protein_index, int):
            self.protein_index = [protein_index]
        else:
            raise TypeError("protein_index must be an int or a list of ints")

    def getLigand1Index(self):
        """
        Get the index of the bound ligand in the system

        Returns
        -------
        int
            The index of the bound ligand in the system.
        """
        return self.ligand1_index

    def _setLigand1Index(self, ligand1_index):
        """
        Set the index of the bound ligand in the system

        Parameters
        ----------
        ligand1_index : int
            The index of the bound ligand in the system.
        """
        if not isinstance(ligand1_index, int):
            raise TypeError("ligand1_index must be an int")
        else:
            self.ligand1_index = ligand1_index

    def getLigand2Index(self):
        """
        Get the index of the free ligand in the system

        Returns
        -------
        int
            The index of the free ligand in the system.
        """
        return self.ligand2_index

    def _setLigand2Index(self, ligand2_index):
        """
        Set the index of the free ligand in the system

        Parameters
        ----------
        ligand2_index : int
            The index of the free ligand in the system.
        """
        if not isinstance(ligand2_index, int):
            raise TypeError("ligand2_index must be an int")
        else:
            self.ligand2_index = ligand2_index

    def getLigand1RigidCore(self):
        """
        Get the user-defined rigid core atom indices for the bound ligand

        Returns
        -------
        list
            A list of three atom indices that define the rigid core of the bound ligand.
        """
        return self.ligand1_rigid_core

    def setLigand1RigidCore(self, ligand1_rigid_core):
        """
        Set the user-defined rigid core atom indices for the bound ligand
        """
        if not ligand1_rigid_core:
            self.ligand1_rigid_core = None
        else:
            if not isinstance(ligand1_rigid_core, list):
                raise TypeError("ligand1_rigid_core must be a list")
            if len(ligand1_rigid_core) != 3:
                raise ValueError("ligand1_rigid_core must have length 3")
            if any(x >= self.ligand1_atomcount for x in ligand1_rigid_core):
                raise ValueError(
                    "ligand1_rigid_core contains an index that is greater than the number of atoms in the ligand"
                )
            self.ligand1_rigid_core = ligand1_rigid_core

    def getLigand2RigidCore(self):
        """
        Get the user-defined rigid core atom indices for the free ligand

        Returns
        -------
        list
            A list of three atom indices that define the rigid core of the free ligand.
        """
        return self.ligand2_rigid_core

    def setLigand2RigidCore(self, ligand2_rigid_core):
        """
        Set the user-defined rigid core atom indices for the free ligand
        """
        if not ligand2_rigid_core:
            self.ligand2_rigid_core = None
        else:
            if not isinstance(ligand2_rigid_core, list):
                raise TypeError("ligand2_rigid_core must be a list")
            if len(ligand2_rigid_core) != 3:
                raise ValueError("ligand2_rigid_core must have length 3")
            if any(x >= self.ligand2_atomcount for x in ligand2_rigid_core):
                raise ValueError(
                    "ligand2_rigid_core contains an index that is greater than the number of atoms in the ligand"
                )
            self.ligand2_rigid_core = ligand2_rigid_core

    def _setMol1(self, mol1):
        """
        Set the first molecule/system

        Parameters
        ----------
        mol1 : BioSimSpace._SireWrappers.Molecule, BioSimSpace._SireWrappers.System
            The first molecule/system - this can be one of the following:
            - A single protein (BioSimSpace._SireWrappers.Molecule)
            - A pre-prepared protein-ligand-ligand system ready for use in AToM (BioSimSpace._SireWrappers.System).
        """
        if not isinstance(mol1, (_Molecule, _System)):
            raise TypeError("mol1 must be a Molecule or System")

        if isinstance(mol1, _Molecule):
            if mol1.isWater():
                print("Looks like mol1 is a water molecule - it should be a protein")
            self.mol1_atomcount = mol1.nAtoms()
            self.mol1 = mol1
        elif isinstance(mol1, _System):
            print("Assuming that mol1 is a pre-prepared AToM system")
            self._is_pre_prepared = True
            self.system = mol1
            self._systemInfo()

    def _systemInfo(self):
        """
        If the user gives a pre-prepared AToM system, extract the needed information

        Returns
        -------
        dict
            A dictionary containing information on the AToM system
        """
        for p in self.protein_index:
            if self.system[p].isWater():
                print(
                    f"The molecule at index {self.protein_index} appears to be a water molecule."
                    " This should be a protein."
                )
        if self.system[self.ligand1_index].isWater():
            print(
                f"The molecule at index {self.ligand1_index} appears to be a water molecule."
                " This should be the bound ligand."
            )
        if self.system[self.ligand2_index].isWater():
            print(
                f"The molecule at index {self.ligand2_index} appears to be a water molecule."
                " This should be the free ligand."
            )
        self.mol1_atomcount = sum(self.system[i].nAtoms() for i in self.protein_index)
        self.ligand1_atomcount = self.system[self.ligand1_index].nAtoms()
        self.ligand2_atomcount = self.system[self.ligand2_index].nAtoms()

    def _setLigand1(self, ligand1):
        if self._is_pre_prepared:
            print("Pre-prepared system given...Ignoring ligand1")
        else:
            if not isinstance(ligand1, _Molecule):
                raise TypeError("ligand1 must be a Molecule")
            if ligand1.isWater():
                print(
                    "Looks like ligand1 is a water molecule - it should the bound ligand"
                )
            self.ligand1 = ligand1
            self.ligand1_atomcount = ligand1.nAtoms()

    def _setLigand2(self, ligand2):
        if self._is_pre_prepared:
            print("Pre-prepared system given...Ignoring ligand2")
        else:
            if not isinstance(ligand2, _Molecule):
                raise TypeError("ligand2 must be a Molecule")
            if ligand2.isWater():
                print(
                    "Looks like ligand2 is a water molecule - it should the free ligand"
                )
            self.ligand2 = ligand2
            self.ligand2_atomcount = ligand2.nAtoms()

    def _setDisplacement(self, displacement):
        if isinstance(displacement, str):
            self.displacement = _Length(displacement)
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
        else:
            raise TypeError("displacement must be a string, BSS length or list")

    def _makeSystemFromThree(self):
        """
        Make the AToM system from a protein and two ligands
        """
        # First step is to align the ligands
        # align ligand2 to ligand1
        mapping = _matchAtoms(self.ligand2, self.ligand1)
        self.ligand2 = _rmsdAlign(self.ligand2, self.ligand1, mapping)

        # Concatenate the molecules
        system_ligand1 = (self.mol1 + self.ligand1).toSystem()

        if isinstance(self.displacement, _Vector):
            self.ligand2.translate(
                [self.displacement.x(), self.displacement.y(), self.displacement.z()]
            )
            self.data["displacement"] = self.displacement
        else:
            # Now we need to translate ligand2 so that it is separated from the protein/ligand1
            # by the desired distance
            vec = self._findTranslationVector(
                system_ligand1,
                self.displacement,
                protein=self.mol1,
                ligand=self.ligand1,
            )
            self.ligand2.translate([vec.x(), vec.y(), vec.z()])
            self.data["displacement"] = [vec.x(), vec.y(), vec.z()]

        self.system = (self.mol1 + self.ligand1 + self.ligand2).toSystem()
        self._setProteinIndex(self.system.getIndex(self.mol1))
        self._setLigand1Index(self.system.getIndex(self.ligand1))
        self._setLigand2Index(self.system.getIndex(self.ligand2))

    def _findAtomIndices(self):
        """
        Find the indices of the protein and ligand atoms in the system

        Returns
        -------
        dict
            A dictionary containing the indices of the protein and ligand atoms in the system
        """
        print(f"CURRENT INDEX: {self.protein_index}")
        protein_atom_start = self.system[self.protein_index[0]].getAtoms()[0]
        protein_atom_end = self.system[self.protein_index[-1]].getAtoms()[-1]
        self.first_protein_atom_index = self.system.getIndex(protein_atom_start)
        self.last_protein_atom_index = self.system.getIndex(protein_atom_end)

        ligand1_atom_start = self.system[self.ligand1_index].getAtoms()[0]
        ligand1_atom_end = self.system[self.ligand1_index].getAtoms()[-1]
        self.first_ligand1_atom_index = self.system.getIndex(ligand1_atom_start)
        self.last_ligand1_atom_index = self.system.getIndex(ligand1_atom_end)

        ligand2_atom_start = self.system[self.ligand2_index].getAtoms()[0]
        ligand2_atom_end = self.system[self.ligand2_index].getAtoms()[-1]
        self.first_ligand2_atom_index = self.system.getIndex(ligand2_atom_start)
        self.last_ligand2_atom_index = self.system.getIndex(ligand2_atom_end)

    @staticmethod
    def _findTranslationVector(
        system,
        displacement,
        protein=None,
        ligand=None,
        protein_index=None,
        ligand1_index=None,
    ):
        """
        Finds the vector along which the free ligand is to be translated.
        Based on the funnel making logic of biosimspace

        Parameters
        ----------
        system : BioSimSpace._SireWrappers.System
            The system containing the protein and ligand.

        protein : BioSimSpace._SireWrappers.Molecule
            The protein molecule.
        ligand : BioSimSpace._SireWrappers.Molecule
            The ligand molecule.
        displacement : float
            The desired displacement between the ligands in angstroms.

        Returns
        -------
        BioSimSpace.Types.Vector
            The vector along which the ligand is to be translated.
        """
        from sire.legacy.Maths import Vector

        if not isinstance(system, _System):
            raise TypeError("system must be a BioSimSpace system")
        if not isinstance(protein, (_Molecule, type(None))):
            raise TypeError("protein must be a BioSimSpace molecule")
        if not isinstance(ligand, (_Molecule, type(None))):
            raise TypeError("ligand must be a BioSimSpace molecule")

        if protein is None:
            protein = system[protein_index]
        if ligand is None:
            ligand = system[ligand1_index]
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
            for y in np.linspace(grid_min.y().value(), grid_max.y().value(), num_edges):
                for z in np.linspace(
                    grid_min.z().value(), grid_max.z().value(), num_edges
                ):
                    search = f"atoms within {search_radius.value()} of ({x}, {y}, {z})"

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

    def get_prot_com_atoms(self):
        """
        Get the atoms that define the center of mass of the protein as a list of ints

        Returns
        -------
        list
            A list of atom indices that define the center of mass of the protein.
        """
        return self._mol1_com_atoms

    def _set_mol1_com_atoms(self, mol1_com_atoms):
        """
        Set the atoms that define the center of mass of the protein
        If a list is given, simply set them according to the list.
        If None, find them based on the center of mass of the protein.
        """
        if mol1_com_atoms is not None:
            # Make sure its a list of ints
            if not isinstance(mol1_com_atoms, list):
                raise TypeError("mol1_com_atoms must be a list")
            if not all(isinstance(x, int) for x in mol1_com_atoms):
                raise TypeError("mol1_com_atoms must be a list of ints")
            self._mol1_com_atoms = mol1_com_atoms
        else:
            # Find com of the protein
            if self._is_pre_prepared:
                temp_system = self.system._sire_object
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
                protein = self.mol1
                com = protein._sire_object.coordinates()
                self._mol1_com_atoms = [
                    a.index().value()
                    for a in protein._sire_object[f"atoms within 11 angstrom of {com}"]
                ]

    def get_lig1_com_atoms(self):
        """
        Get the atoms that define the center of mass of the bound ligand as a list of ints

        Returns
        -------
        list
            A list of atom indices that define the center of mass of the bound ligand.
        """
        return self._lig1_com_atoms

    def _set_lig1_com_atoms(self, lig1_com_atoms):
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
            if self._is_pre_prepared:
                ligand1 = self.system[self.ligand1_index]
            else:
                ligand1 = self.ligand1
            com = ligand1._sire_object.coordinates()
            self._lig1_com_atoms = [
                a.index().value()
                for a in ligand1._sire_object[f"atoms within 11 angstrom of {com}"]
            ]

    def get_lig2_com_atoms(self):
        """
        Get the atoms that define the center of mass of the free ligand as a list of ints

        Returns
        -------
        list
            A list of atom indices that define the center of mass of the free ligand.
        """
        return self._lig2_com_atoms

    def _set_lig2_com_atoms(self, lig2_com_atoms):
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
            if self._is_pre_prepared:
                ligand2 = self.system[self.ligand2_index]
            else:
                ligand2 = self.ligand2
            com = ligand2._sire_object.coordinates()
            self._lig2_com_atoms = [
                a.index().value()
                for a in ligand2._sire_object[f"atoms within 11 angstrom of {com}"]
            ]

    def getSystem(self):
        """
        Get the AToM system

        Returns
        -------
        BioSimSpace._SireWrappers.System
            The AToM system
        """
        self._makeData()
        return self.system, self.data

    def _makeData(self):
        """
        Make the data dictionary for the AToM system
        """
        self._findAtomIndices()
        if self._is_pre_prepared:
            if not isinstance(self.displacement, _Vector):
                raise ValueError(
                    "Displacement must be a vector or list if a pre-prepared system is given"
                )
            else:
                self.data["displacement"] = self.displacement
        self.data["protein_index"] = self.protein_index
        self.data["ligand1_index"] = self.ligand1_index
        self.data["ligand2_index"] = self.ligand2_index
        self.data["ligand1_rigid_core"] = self.ligand1_rigid_core
        self.data["ligand2_rigid_core"] = self.ligand2_rigid_core
        self.data["mol1_atomcount"] = self.mol1_atomcount
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


class relativeATM:
    """
    A class for setting up and performing RBFE calculations using AToM
    """

    def __init__(
        self,
        system,
        protocol=None,
        platform="CPU",
        work_dir=None,
        setup_only=False,
        ignore_warnings=False,
        show_errors=True,
        extra_options={},
        extra_lines=[],
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

        ignore_warnings : bool
            Whether to ignore warnings when generating the binary run file.
            This option is specific to GROMACS and will be ignored when a
            different molecular dynamics engine is chosen.

        show_errors : bool
            Whether to show warning/error messages when generating the binary
            run file. This option is specific to GROMACS and will be ignored
            when a different molecular dynamics engine is chosen.

        extra_options : dict
            A dictionary containing extra options. Overrides the defaults generated
            by the protocol.

        extra_lines : [str]
            A list of extra lines to put at the end of the configuration file.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        """
        if not isinstance(system, _System):
            raise TypeError(
                "'system' must be of type 'BioSimSpace._SireWrappers.System'"
            )
        else:
            # Store a copy of solvated system.
            self._system = system.copy()

        # Validate the protocol.
        if protocol is not None:
            from ...Protocol._AToM import AToMProduction as _Production

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

        if not isinstance(ignore_warnings, bool):
            raise ValueError("'ignore_warnings' must be of type 'bool.")
        self._ignore_warnings = ignore_warnings

        if not isinstance(show_errors, bool):
            raise ValueError("'show_errors' must be of type 'bool.")
        self._show_errors = show_errors

        # Check the extra options.
        if not isinstance(extra_options, dict):
            raise TypeError("'extra_options' must be of type 'dict'.")
        else:
            keys = extra_options.keys()
            if not all(isinstance(k, str) for k in keys):
                raise TypeError("Keys of 'extra_options' must be of type 'str'.")
        self._extra_options = extra_options

        # Check the extra lines.
        if not isinstance(extra_lines, list):
            raise TypeError("'extra_lines' must be of type 'list'.")
        else:
            if not all(isinstance(line, str) for line in extra_lines):
                raise TypeError("Lines in 'extra_lines' must be of type 'str'.")
        self._extra_lines = extra_lines

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
        list of :class:`Process <BioSimSpace.Process>`
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
        # TODO: generate generic protocol if None is passed.
        # This protocol will have to be minimal - cannot guess rigid core atoms
        if self._protocol is None:
            raise RuntimeError("No protocol has been set - cannot run simulations.")
        # Initialise list to store the processe
        processes = []
        # Get the list of lambda1 values so that the total number of simulations can
        # be asserted
        lambda_list = self._protocol._get_lambda_values()
        # Set index of current simulation to 0
        self._protocol._set_current_index(0)
        lam = lambda_list[0]

        first_dir = "%s/lambda_%5.4f" % (self._work_dir, lam)

        # Create the first simulation, which will be copied and used for future simulations.
        first_process = _Production(
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
            # TODO: Support for simulations restarting from a checkpoint.
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
            new_direction = self._protocol.getDirections()[index]
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
            self._runner = _Process.ProcessRunner(processes)


def viewRigidCores(system, data):
    """
    View the aligned ligands with rigid core atoms defined by the user during system creation.
    """

    if (
        "ligand1_rigid_core" not in data.keys()
        or "ligand2_rigid_core" not in data.keys()
    ):
        raise ValueError(
            "ligand1_rigid_core and ligand2_rigid_core must be defined in data"
        )

    def move_to_origin(lig):
        com = _Coordinate(*lig._getCenterOfMass())
        lig.translate([-com.x().value(), -com.y().value(), -com.z().value()])

    if not _is_notebook:
        raise RuntimeError("This function can only be used in a Jupyter notebook")

    if not isinstance(system, _System):
        raise TypeError("system must be a BioSimSpace system")

    if not isinstance(data, dict):
        raise TypeError("data must be a dictionary")

    # copy the ligands
    ligand1 = system[data["ligand1_index"]].copy()
    move_to_origin(ligand1)
    ligand2 = system[data["ligand2_index"]].copy()
    move_to_origin(ligand2)

    # Translate ligand2 so they don't overlap
    ligand2.translate([10.0, 0, 0])

    # Get coords of rigid core atoms
    ligand1_core_coords = []
    ligand2_core_coords = []

    for i in data["ligand1_rigid_core"]:
        ligand1_core_coords.append(ligand1.getAtoms()[i].coordinates())

    for i in data["ligand2_rigid_core"]:
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
            [coord1.x().value(), coord1.y().value(), coord1.z().value()], colour, 0.7
        )
        ngl.shape.add_sphere(
            [coord2.x().value(), coord2.y().value(), coord2.z().value()], colour, 0.7
        )
    return ngl


class Minimise(_Process.OpenMM):
    """A class for minimising AToM systems.
    Includes restraints."""

    def __init__(
        self,
        system,
        protocol,
        exe=None,
        platform="CPU",
        seed=None,
        work_dir=None,
        property_map={},
    ):
        # Validate the protocol.
        if protocol is not None:
            from ...Protocol._AToM import AToMMinimisation as _Minimisation

            if not isinstance(protocol, _Minimisation):
                raise TypeError(
                    "'protocol' must be of type 'BioSimSpace.Protocol.AToM'"
                )
            else:
                self._protocol = protocol
        else:
            # No default protocol due to the need for well-defined rigid cores
            raise ValueError("A protocol must be specified")
        super().__init__(
            system=system,
            protocol=protocol,
            exe=exe,
            name="minimise",
            platform=platform,
            work_dir=work_dir,
            seed=seed,
            property_map=property_map,
        )

    def _generate_config(self):
        """Custom config generation, mixing in restraints
        with minimisation protocol. AToM force is not applied here."""

        util = _AtomUtils(self._protocol)
        # Clear the existing configuration list.
        self._config = []

        # Get the "space" property from the user mapping.
        prop = self._property_map.get("space", "space")

        # Check whether the system contains periodic box information.
        if prop in self._system._sire_object.propertyKeys():
            try:
                # Make sure that we have a periodic box. The system will now have
                # a default cartesian space.
                box = self._system._sire_object.property(prop)
                has_box = box.isPeriodic()
            except:
                has_box = False
        else:
            _warnings.warn("No simulation box found. Assuming gas phase simulation.")
            has_box = False

        # Write the OpenMM import statements.
        self._add_config_imports()
        self._add_config_monkey_patches()
        # Load the input files.
        self.addToConfig("\n# Load the topology and coordinate files.")
        self.addToConfig(
            "\n# We use ParmEd due to issues with the built in AmberPrmtopFile for certain triclinic spaces."
        )
        self.addToConfig(
            f"prm = parmed.load_file('{self._name}.prm7', '{self._name}.rst7')"
        )

        # Don't use a cut-off if this is a vacuum simulation or if box information
        # is missing.
        self.addToConfig("\n# Initialise the molecular system.")
        if not has_box or not self._has_water:
            self.addToConfig("system = prm.createSystem(nonbondedMethod=NoCutoff,")
        else:
            self.addToConfig("system = prm.createSystem(nonbondedMethod=PME,")
        self.addToConfig("                          nonbondedCutoff=1*nanometer,")
        self.addToConfig("                          constraints=HBonds)")

        # Set the integrator. (Use zero-temperature as this is just a dummy step.)
        self.addToConfig("\n# Define the integrator.")
        self.addToConfig("integrator = LangevinMiddleIntegrator(0*kelvin,")
        self.addToConfig("                                1/picosecond,")
        self.addToConfig("                                0.002*picoseconds)")

        # Add the platform information.
        self._add_config_platform()

        # Add any position restraints.
        if self._protocol.getRestraint() is not None:
            restraint = self._protocol.getRestraint()
            # Search for the atoms to restrain by keyword.
            if isinstance(restraint, str):
                restrained_atoms = self._system.getRestraintAtoms(restraint)
            # Use the user-defined list of indices.
            else:
                restrained_atoms = restraint
            self.addToConfig("\n# Add position restraints.")
            frc = util.create_flat_bottom_restraint(restrained_atoms)
            self.addToConfig(frc)

        # Add the atom-specific restraints.
        disp = util.createDisplacement()
        self.addToConfig(disp)
        if self._protocol.getCoreAlignment():
            alignment = util.createAlignmentForce()
            self.addToConfig("\n# Add alignment force.")
            self.addToConfig(alignment)
        if self._protocol.getCMCMRestraint():
            CMCM = util.createCOMRestraint()
            self.addToConfig("\n# Add COM restraint.")
            self.addToConfig(CMCM)

        # Set up the simulation object.
        self.addToConfig("\n# Initialise and configure the simulation object.")
        self.addToConfig("simulation = Simulation(prm.topology,")
        self.addToConfig("                        system,")
        self.addToConfig("                        integrator,")
        self.addToConfig("                        platform,")
        self.addToConfig("                        properties)")
        if self._protocol.getRestraint() is not None:
            self.addToConfig("simulation.context.setPositions(positions)")
        else:
            self.addToConfig("simulation.context.setPositions(prm.positions)")
        self.addToConfig("if prm.box_vectors is not None:")
        self.addToConfig("    box_vectors = reducePeriodicBoxVectors(prm.box_vectors)")
        self.addToConfig("    simulation.context.setPeriodicBoxVectors(*box_vectors)")
        self.addToConfig(
            f"simulation.minimizeEnergy(maxIterations={self._protocol.getSteps()})"
        )

        # Add the reporters.
        self.addToConfig("\n# Add reporters.")
        self._add_config_reporters(state_interval=1, traj_interval=1)

        # Now run the simulation.
        self.addToConfig(
            "\n# Run a single simulation step to allow us to get the system and energy."
        )
        self.addToConfig(f"simulation.step(1)")

        # Flag that this isn't a custom protocol.
        self._protocol._setCustomised(False)


class Equilibrate(_Process.OpenMM):
    """A class for running a pre-production equilibration step.
    Required for most systems to preven atom overlapping issues."""

    def __init__(
        self,
        system,
        protocol,
        exe=None,
        platform="CUDA",
        seed=None,
        work_dir=None,
        property_map={},
    ):
        # Validate the protocol.
        if protocol is not None:
            from ...Protocol._AToM import AToMEquilibration as _Equilibration

            if not isinstance(protocol, _Equilibration):
                raise TypeError(
                    "'protocol' must be of type 'BioSimSpace.Protocol.AToMEquilibration'"
                )
            else:
                self._protocol = protocol
        else:
            # No default protocol due to the need for well-defined rigid cores
            raise ValueError("A protocol must be specified")
        super().__init__(
            system=system,
            protocol=protocol,
            exe=exe,
            name="equilibrate",
            platform=platform,
            work_dir=work_dir,
            seed=seed,
            property_map=property_map,
        )

    def _generate_config(self):

        util = _AtomUtils(self._protocol)
        # Clear the existing configuration list.
        self._config = []

        # Get the "space" property from the user mapping.
        prop = self._property_map.get("space", "space")

        # Check whether the system contains periodic box information.
        if prop in self._system._sire_object.propertyKeys():
            try:
                # Make sure that we have a periodic box. The system will now have
                # a default cartesian space.
                box = self._system._sire_object.property(prop)
                has_box = box.isPeriodic()
            except:
                has_box = False
        else:
            _warnings.warn("No simulation box found. Assuming gas phase simulation.")
            has_box = False

        # Write the OpenMM import statements and monkey-patches.
        self._add_config_imports()
        self._add_config_monkey_patches()

        # Load the input files.
        self.addToConfig("\n# Load the topology and coordinate files.")
        self.addToConfig(
            "\n# We use ParmEd due to issues with the built in AmberPrmtopFile for certain triclinic spaces."
        )
        self.addToConfig(
            f"prm = parmed.load_file('{self._name}.prm7', '{self._name}.rst7')"
        )

        # Don't use a cut-off if this is a vacuum simulation or if box information
        # is missing.
        is_periodic = True
        self.addToConfig("\n# Initialise the molecular system.")
        if not has_box or not self._has_water:
            is_periodic = False
            self.addToConfig("system = prm.createSystem(nonbondedMethod=NoCutoff,")
        else:
            self.addToConfig("system = prm.createSystem(nonbondedMethod=PME,")
        self.addToConfig("                          nonbondedCutoff=1*nanometer,")
        self.addToConfig("                          constraints=HBonds)")

        # Get the starting temperature and system pressure.
        temperature = self._protocol.getStartTemperature().kelvin().value()
        pressure = self._protocol.getPressure()

        # Add a Monte Carlo barostat if the simulation is at constant pressure.
        is_const_pressure = False
        if pressure is not None:
            # Cannot use a barostat with a non-periodic system.
            if not is_periodic:
                _warnings.warn(
                    "Cannot use a barostat for a vacuum or non-periodic simulation"
                )
            else:
                is_const_pressure = True

                # Convert to bar and get the value.
                pressure = pressure.bar().value()

                # Create the barostat and add its force to the system.
                self.addToConfig("\n# Add a barostat to run at constant pressure.")
                self.addToConfig(
                    f"barostat = MonteCarloBarostat({pressure}*bar, {temperature}*kelvin)"
                )
                if self._is_seeded:
                    self.addToConfig(f"barostat.setRandomNumberSeed({self._seed})")
                self.addToConfig("system.addForce(barostat)")

        # Add any position restraints.
        if self._protocol.getRestraint() is not None:
            restraint = self._protocol.getRestraint()
            # Search for the atoms to restrain by keyword.
            if isinstance(restraint, str):
                restrained_atoms = self._system.getRestraintAtoms(restraint)
            # Use the user-defined list of indices.
            else:
                restrained_atoms = restraint
            self.addToConfig("\n# Add position restraints.")
            frc = util.create_flat_bottom_restraint(restrained_atoms)
            self.addToConfig(frc)

            # Add the atom-specific restraints.
        disp = util.createDisplacement()
        self.addToConfig(disp)
        if self._protocol.getUseATMForce():
            atm = util.createATMForce(index=None)
            self.addToConfig(atm)

        if self._protocol.getCoreAlignment():
            alignment = util.createAlignmentForce()
            self.addToConfig("\n# Add alignment force.")
            self.addToConfig(alignment)
        if self._protocol.getCMCMRestraint():
            CMCM = util.createCOMRestraint()
            self.addToConfig("\n# Add COM restraint.")
            self.addToConfig(CMCM)

        # Get the integration time step from the protocol.
        timestep = self._protocol.getTimeStep().picoseconds().value()

        # Set the integrator.
        self.addToConfig("\n# Define the integrator.")
        self.addToConfig(f"integrator = LangevinMiddleIntegrator({temperature}*kelvin,")
        friction = 1 / self._protocol.getThermostatTimeConstant().picoseconds().value()
        self.addToConfig(f"                                {friction:.5f}/picosecond,")
        self.addToConfig(f"                                {timestep}*picoseconds)")
        if self._is_seeded:
            self.addToConfig(f"integrator.setRandomNumberSeed({self._seed})")

        # Add the platform information.
        self._add_config_platform()

        # Set up the simulation object.
        self.addToConfig("\n# Initialise and configure the simulation object.")
        self.addToConfig("simulation = Simulation(prm.topology,")
        self.addToConfig("                        system,")
        self.addToConfig("                        integrator,")
        self.addToConfig("                        platform,")
        self.addToConfig("                        properties)")
        if self._protocol.getRestraint() is not None:
            self.addToConfig("simulation.context.setPositions(positions)")
        else:
            self.addToConfig("simulation.context.setPositions(prm.positions)")
        self.addToConfig("if prm.box_vectors is not None:")
        self.addToConfig("    box_vectors = reducePeriodicBoxVectors(prm.box_vectors)")
        self.addToConfig("    simulation.context.setPeriodicBoxVectors(*box_vectors)")

        # Set initial velocities from temperature distribution.
        self.addToConfig("\n# Setting initial system velocities.")
        self.addToConfig(
            f"simulation.context.setVelocitiesToTemperature({temperature})"
        )

        # Work out the number of integration steps.
        steps = _math.ceil(self._protocol.getRunTime() / self._protocol.getTimeStep())

        # Get the report and restart intervals.
        report_interval = self._protocol.getReportInterval()
        restart_interval = self._protocol.getRestartInterval()

        # Cap the intervals at the total number of steps.
        if report_interval > steps:
            report_interval = steps
        if restart_interval > steps:
            restart_interval = steps

        # Add the reporters.
        self.addToConfig("\n# Add reporters.")
        self._add_config_reporters(
            state_interval=report_interval,
            traj_interval=restart_interval,
            is_restart=False,
        )

        # Now run the simulation.
        self.addToConfig("\n# Run the simulation.")

        # Constant temperature equilibration.
        if self._protocol.isConstantTemp():
            self.addToConfig(f"simulation.step({steps})")

        # Heating / cooling cycle.
        else:
            # Adjust temperature every 100 cycles, assuming that there at
            # least that many cycles.
            if steps > 100:
                # Work out the number of temperature cycles.
                temp_cycles = _math.ceil(steps / 100)

                # Work out the temperature change per cycle.
                delta_temp = (
                    self._protocol.getEndTemperature().kelvin().value()
                    - self._protocol.getStartTemperature().kelvin().value()
                ) / temp_cycles

                self.addToConfig(f"start_temperature = {temperature}")
                self.addToConfig(f"for x in range(0, {temp_cycles}):")
                self.addToConfig(f"    temperature = {temperature} + x*{delta_temp}")
                self.addToConfig(f"    integrator.setTemperature(temperature*kelvin)")
                if is_const_pressure:
                    self.addToConfig(
                        f"    barostat.setDefaultTemperature(temperature*kelvin)"
                    )
                self.addToConfig("    simulation.step(100)")
            else:
                # Work out the temperature change per step.
                delta_temp = (
                    self._protocol.getEndTemperature().kelvin().value()
                    - self._protocol.getStartTemperature().kelvin().value()
                ) / steps

                self.addToConfig(f"start_temperature = {temperature}")
                self.addToConfig(f"for x in range(0, {steps}):")
                self.addToConfig(f"    temperature = {temperature} + x*{delta_temp}")
                self.addToConfig(f"    integrator.setTemperature(temperature*kelvin)")
                if is_const_pressure:
                    self.addToConfig(
                        f"    barostat.setDefaultTemperature(temperature*kelvin)"
                    )
                self.addToConfig("    simulation.step(1)")


class Anneal(_Process.OpenMM):
    """A class for pre-production annealing of AToM systems."""

    def __init__(
        self,
        system,
        protocol,
        exe=None,
        platform="CUDA",
        seed=None,
        work_dir=None,
        property_map={},
    ):
        # Validate the protocol.
        if protocol is not None:
            from ...Protocol._AToM import AToMAnnealing as _Annealing

            if not isinstance(protocol, _Annealing):
                raise TypeError(
                    "'protocol' must be of type 'BioSimSpace.Protocol.AToMEquilibration'"
                )
            else:
                self._protocol = protocol
        else:
            # No default protocol due to the need for well-defined rigid cores
            raise ValueError("A protocol must be specified")

        if protocol.getAnnealValues() is None or protocol.getAnnealNumCycles is None:
            raise ValueError("Anneal values and cycles must be specified")

        # Use this to set default values for the protocol
        protocol._set_current_index(0)
        super().__init__(
            system=system,
            protocol=protocol,
            exe=exe,
            name="anneal",
            platform=platform,
            work_dir=work_dir,
            seed=seed,
            property_map=property_map,
        )

    def _generate_config(self):
        util = _AtomUtils(self._protocol)
        # Clear the existing configuration list.
        self._config = []

        # Get the "space" property from the user mapping.
        prop = self._property_map.get("space", "space")

        # Check whether the system contains periodic box information.
        if prop in self._system._sire_object.propertyKeys():
            try:
                # Make sure that we have a periodic box. The system will now have
                # a default cartesian space.
                box = self._system._sire_object.property(prop)
                has_box = box.isPeriodic()
            except:
                has_box = False
        else:
            _warnings.warn("No simulation box found. Assuming gas phase simulation.")
            has_box = False

        # Write the OpenMM import statements.
        self._add_config_imports()
        self._add_config_monkey_patches()

        # Add standard openMM config
        self.addToConfig("from glob import glob")
        self.addToConfig("import math")
        self.addToConfig("import os")
        self.addToConfig("import shutil")

        # Load the input files.
        self.addToConfig("\n# Load the topology and coordinate files.")
        self.addToConfig(
            "\n# We use ParmEd due to issues with the built in AmberPrmtopFile for certain triclinic spaces."
        )
        self.addToConfig(
            f"prm = parmed.load_file('{self._name}.prm7', '{self._name}.rst7')"
        )

        # Don't use a cut-off if this is a vacuum simulation or if box information
        # is missing.
        is_periodic = True
        self.addToConfig("\n# Initialise the molecular system.")
        if not has_box or not self._has_water:
            is_periodic = False
            self.addToConfig("system = prm.createSystem(nonbondedMethod=NoCutoff,")
        else:
            self.addToConfig("system = prm.createSystem(nonbondedMethod=PME,")
        self.addToConfig("                          nonbondedCutoff=1*nanometer,")
        self.addToConfig("                          constraints=HBonds)")

        # Get the starting temperature and system pressure.
        temperature = self._protocol.getTemperature().kelvin().value()
        pressure = self._protocol.getPressure()

        # Add a Monte Carlo barostat if the simulation is at constant pressure.
        is_const_pressure = False
        if pressure is not None:
            # Cannot use a barostat with a non-periodic system.
            if not is_periodic:
                _warnings.warn(
                    "Cannot use a barostat for a vacuum or non-periodic simulation"
                )
            else:
                is_const_pressure = True

                # Convert to bar and get the value.
                pressure = pressure.bar().value()

                # Create the barostat and add its force to the system.
                self.addToConfig("\n# Add a barostat to run at constant pressure.")
                self.addToConfig(
                    f"barostat = MonteCarloBarostat({pressure}*bar, {temperature}*kelvin)"
                )
                if self._is_seeded:
                    self.addToConfig(f"barostat.setRandomNumberSeed({self._seed})")
                self.addToConfig("system.addForce(barostat)")

        # Add any position restraints.
        if self._protocol.getRestraint() is not None:
            restraint = self._protocol.getRestraint()
            # Search for the atoms to restrain by keyword.
            if isinstance(restraint, str):
                restrained_atoms = self._system.getRestraintAtoms(restraint)
            # Use the user-defined list of indices.
            else:
                restrained_atoms = restraint
            self.addToConfig("\n# Add position restraints.")
            frc = util.create_flat_bottom_restraint(restrained_atoms)
            self.addToConfig(frc)

        # Get the integration time step from the protocol.
        timestep = self._protocol.getTimeStep().picoseconds().value()

        # Set the integrator.
        self.addToConfig("\n# Define the integrator.")
        self.addToConfig(f"integrator = LangevinMiddleIntegrator({temperature}*kelvin,")
        friction = 1 / self._protocol.getThermostatTimeConstant().picoseconds().value()
        self.addToConfig(f"                                {friction:.5f}/picosecond,")
        self.addToConfig(f"                                {timestep}*picoseconds)")
        if self._is_seeded:
            self.addToConfig(f"integrator.setRandomNumberSeed({self._seed})")

        # Add the platform information.
        self._add_config_platform()

        # Use utils to create AToM-specific forces
        # Atom force is the only window-dependent force
        disp = util.createDisplacement()
        self.addToConfig(disp)
        self.addToConfig("\n# Add AToM Force.")
        self.addToConfig(util.createATMForce(self._protocol._get_window_index()))
        if self._protocol.getCoreAlignment():
            alignment = util.createAlignmentForce()
            self.addToConfig("\n# Add alignment force.")
            self.addToConfig(alignment)

        if self._protocol.getCMCMRestraint():
            CMCM = util.createCOMRestraint()
            self.addToConfig("\n# Add COM restraint.")
            self.addToConfig(CMCM)

        self.addToConfig("\n# Initialise and configure the simulation object.")
        self.addToConfig("simulation = Simulation(prm.topology,")
        self.addToConfig("                        system,")
        self.addToConfig("                        integrator,")
        self.addToConfig("                        platform,")
        self.addToConfig("                        properties)")
        if self._protocol.getRestraint() is not None:
            self.addToConfig("simulation.context.setPositions(positions)")
        else:
            self.addToConfig("simulation.context.setPositions(prm.positions)")
        self.addToConfig("if prm.box_vectors is not None:")
        self.addToConfig("    box_vectors = reducePeriodicBoxVectors(prm.box_vectors)")
        self.addToConfig("    simulation.context.setPeriodicBoxVectors(*box_vectors)")

        # Set initial velocities from temperature distribution.
        self.addToConfig("\n# Setting initial system velocities.")
        self.addToConfig(
            f"simulation.context.setVelocitiesToTemperature({temperature})"
        )

        # Check for a restart file and load the simulation state.
        is_restart, step = self._add_config_restart()

        # Work out the number of integration steps.
        total_steps = _math.ceil(
            self._protocol.getRunTime() / self._protocol.getTimeStep()
        )

        # Subtract the current number of steps.
        steps = total_steps - step

        # Exit if the simulation has already finished.
        if steps <= 0:
            print("The simulation has already finished!")
            return

        # Inform user that a restart was loaded.
        self.addToConfig("\n# Print restart information.")
        self.addToConfig("if is_restart:")
        self.addToConfig(f"    steps = {total_steps}")
        self.addToConfig("    percent_complete = 100 * (step / steps)")
        self.addToConfig("    print('Loaded state from an existing simulation.')")
        self.addToConfig("    print(f'Simulation is {percent_complete}% complete.')")

        # Get the report and restart intervals.
        report_interval = self._protocol.getReportInterval()
        restart_interval = self._protocol.getRestartInterval()

        # Cap the intervals at the total number of steps.
        if report_interval > steps:
            report_interval = steps
        if restart_interval > steps:
            restart_interval = steps

        # Add the reporters.
        self.addToConfig("\n# Add reporters.")
        self._add_config_reporters(
            state_interval=report_interval,
            traj_interval=restart_interval,
            is_restart=is_restart,
        )

        # Work out the total simulation time in picoseconds.
        run_time = steps * timestep

        # Work out the number of cycles in 100 picosecond intervals.
        cycles = _math.ceil(run_time / 100)

        # Work out the number of steps per cycle.
        steps_per_cycle = int(steps / cycles)

        # get annealing protocol from atom utils
        annealing_protocol = util.createAnnealingProtocol()
        self.addToConfig(annealing_protocol)


class _Production(_Process.OpenMM):
    """An internal class for running a single production AToM simulation."""

    def __init__(
        self,
        system,
        protocol,
        exe=None,
        platform="CUDA",
        seed=None,
        work_dir=None,
        property_map={},
    ):
        # Validate the protocol.
        if protocol is not None:
            from ...Protocol._AToM import AToMProduction as _Production

            if not isinstance(protocol, _Production):
                raise TypeError(
                    "'protocol' must be of type 'BioSimSpace.Protocol.AToMEquilibration'"
                )
            else:
                self._protocol = protocol
        else:
            # No default protocol due to the need for well-defined rigid cores
            raise ValueError("A protocol must be specified")

        # Use this to set default values for the protocol
        protocol._set_current_index(0)
        super().__init__(
            system=system,
            protocol=protocol,
            exe=exe,
            name="openmm",
            platform=platform,
            work_dir=work_dir,
            seed=seed,
            property_map=property_map,
        )

    def _generate_config(self):
        # Clear the existing configuration list.
        self._config = []

        # Get the "space" property from the user mapping.
        prop = self._property_map.get("space", "space")

        # Check whether the system contains periodic box information.
        if prop in self._system._sire_object.propertyKeys():
            try:
                # Make sure that we have a periodic box. The system will now have
                # a default cartesian space.
                box = self._system._sire_object.property(prop)
                has_box = box.isPeriodic()
            except:
                has_box = False
        else:
            _warnings.warn("No simulation box found. Assuming gas phase simulation.")
            has_box = False
        # TODO: check extra_options, extra_lines and property_map
        util = _AtomUtils(self._protocol)
        if self._protocol._get_window_index() is None:
            raise _IncompatibleError(
                "AToM protocol requires the current window index to be set."
            )

        # Write the OpenMM import statements.
        self._add_config_imports()
        self.addToConfig("import pandas as pd")
        self.addToConfig("import numpy as np")
        self._add_config_monkey_patches()
        self.addToConfig("\n")
        self.addToConfig(util.createSoftcorePertE())
        # Add standard openMM config
        self.addToConfig("from glob import glob")
        self.addToConfig("import math")
        self.addToConfig("import os")
        self.addToConfig("import shutil")

        # Load the input files.
        self.addToConfig("\n# Load the topology and coordinate files.")
        self.addToConfig(
            "\n# We use ParmEd due to issues with the built in AmberPrmtopFile for certain triclinic spaces."
        )
        self.addToConfig(
            f"prm = parmed.load_file('{self._name}.prm7', '{self._name}.rst7')"
        )

        # Don't use a cut-off if this is a vacuum simulation or if box information
        # is missing.
        is_periodic = True
        self.addToConfig("\n# Initialise the molecular system.")
        if not has_box or not self._has_water:
            is_periodic = False
            self.addToConfig("system = prm.createSystem(nonbondedMethod=NoCutoff,")
        else:
            self.addToConfig("system = prm.createSystem(nonbondedMethod=PME,")
        self.addToConfig("                          nonbondedCutoff=1*nanometer,")
        self.addToConfig("                          constraints=HBonds)")

        # Get the starting temperature and system pressure.
        temperature = self._protocol.getTemperature().kelvin().value()
        pressure = self._protocol.getPressure()

        # Add a Monte Carlo barostat if the simulation is at constant pressure.
        is_const_pressure = False
        if pressure is not None:
            # Cannot use a barostat with a non-periodic system.
            if not is_periodic:
                _warnings.warn(
                    "Cannot use a barostat for a vacuum or non-periodic simulation"
                )
            else:
                is_const_pressure = True

                # Convert to bar and get the value.
                pressure = pressure.bar().value()

                # Create the barostat and add its force to the system.
                self.addToConfig("\n# Add a barostat to run at constant pressure.")
                self.addToConfig(
                    f"barostat = MonteCarloBarostat({pressure}*bar, {temperature}*kelvin)"
                )
                if self._is_seeded:
                    self.addToConfig(f"barostat.setRandomNumberSeed({self._seed})")
                self.addToConfig("system.addForce(barostat)")

        # Add any position restraints.
        if self._protocol.getRestraint() is not None:
            restraint = self._protocol.getRestraint()
            # Search for the atoms to restrain by keyword.
            if isinstance(restraint, str):
                restrained_atoms = self._system.getRestraintAtoms(restraint)
            # Use the user-defined list of indices.
            else:
                restrained_atoms = restraint
            self.addToConfig("\n# Add position restraints.")
            frc = util.create_flat_bottom_restraint(restrained_atoms)
            self.addToConfig(frc)

        # Get the integration time step from the protocol.
        timestep = self._protocol.getTimeStep().picoseconds().value()

        # Set the integrator.
        self.addToConfig("\n# Define the integrator.")
        self.addToConfig(f"integrator = LangevinMiddleIntegrator({temperature}*kelvin,")
        friction = 1 / self._protocol.getThermostatTimeConstant().picoseconds().value()
        self.addToConfig(f"                                {friction:.5f}/picosecond,")
        self.addToConfig(f"                                {timestep}*picoseconds)")
        if self._is_seeded:
            self.addToConfig(f"integrator.setRandomNumberSeed({self._seed})")

        # Add the platform information.
        self._add_config_platform()

        # Use utils to create AToM-specific forces
        # Atom force is the only window-dependent force
        disp = util.createDisplacement()
        self.addToConfig(disp)
        self.addToConfig("\n# Add AToM Force.")
        self.addToConfig(util.createATMForce(self._protocol._get_window_index()))
        if self._protocol.getCoreAlignment():
            alignment = util.createAlignmentForce()
            self.addToConfig("\n# Add alignment force.")
            self.addToConfig(alignment)

        if self._protocol.getCMCMRestraint():
            CMCM = util.createCOMRestraint()
            self.addToConfig("\n# Add COM restraint.")
            self.addToConfig(CMCM)

        self.addToConfig("\n# Initialise and configure the simulation object.")
        self.addToConfig("simulation = Simulation(prm.topology,")
        self.addToConfig("                        system,")
        self.addToConfig("                        integrator,")
        self.addToConfig("                        platform,")
        self.addToConfig("                        properties)")
        if self._protocol.getRestraint() is not None:
            self.addToConfig("simulation.context.setPositions(positions)")
        else:
            self.addToConfig("simulation.context.setPositions(prm.positions)")
        self.addToConfig("if prm.box_vectors is not None:")
        self.addToConfig("    box_vectors = reducePeriodicBoxVectors(prm.box_vectors)")
        self.addToConfig("    simulation.context.setPeriodicBoxVectors(*box_vectors)")

        # Set initial velocities from temperature distribution.
        self.addToConfig("\n# Setting initial system velocities.")
        self.addToConfig(
            f"simulation.context.setVelocitiesToTemperature({temperature})"
        )

        # Check for a restart file and load the simulation state.
        is_restart, step = self._add_config_restart()

        # Work out the number of integration steps.
        total_steps = _math.ceil(
            self._protocol.getRunTime() / self._protocol.getTimeStep()
        )

        # Subtract the current number of steps.
        steps = total_steps - step

        # Exit if the simulation has already finished.
        if steps <= 0:
            print("The simulation has already finished!")
            return

        # Inform user that a restart was loaded.
        self.addToConfig("\n# Print restart information.")
        self.addToConfig("if is_restart:")
        self.addToConfig(f"    steps = {total_steps}")
        self.addToConfig("    percent_complete = 100 * (step / steps)")
        self.addToConfig("    print('Loaded state from an existing simulation.')")
        self.addToConfig("    print(f'Simulation is {percent_complete}% complete.')")

        # Get the report and restart intervals.
        report_interval = self._protocol.getReportInterval()
        restart_interval = self._protocol.getRestartInterval()

        # Cap the intervals at the total number of steps.
        if report_interval > steps:
            report_interval = steps
        if restart_interval > steps:
            restart_interval = steps

        # Add the reporters.
        self.addToConfig("\n# Add reporters.")
        self._add_config_reporters(
            state_interval=report_interval,
            traj_interval=restart_interval,
            is_restart=is_restart,
        )

        # Work out the total simulation time in picoseconds.
        run_time = steps * timestep

        # Work out the number of cycles in 100 picosecond intervals.
        cycles = _math.ceil(run_time / (report_interval * timestep))

        # Work out the number of steps per cycle.
        steps_per_cycle = int(steps / cycles)

        self.addToConfig(f"\ntemperature = {temperature}")
        # Now run the simulation.
        self.addToConfig(
            util.createSoftcorePertELoop(
                self._name, cycles, steps_per_cycle, report_interval, timestep, step
            )
        )
        # Now run the simulation.
        # self.addToConfig("\n# Run the simulation in 100 picosecond cycles.")
        # self.addToConfig(f"for x in range(0, {cycles}):")
        # self.addToConfig(f"    simulation.step({steps_per_cycle})")
        # self.addToConfig(f"    simulation.saveState('{self._name}.xml')")
