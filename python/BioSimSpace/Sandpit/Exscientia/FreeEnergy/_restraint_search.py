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

"""
Functionality for selecting receptor-ligand restraints from unrestrained
simulations.
"""

__author__ = "Finlay Clark"
__email__ = "finlay.clark@ed.ac.uk"

__all__ = ["RestraintSearch"]

import matplotlib.pyplot as _plt
import numpy as _np
import os as _os
from scipy.stats import circmean as _circmean
from sklearn.cluster import KMeans as _KMeans
import warnings as _warnings

from sire.legacy.Base import getBinDir as _getBinDir
from sire.legacy.Base import getShareDir as _getShareDir
from sire.legacy.Units import k_boltz as _k_boltz  # kcal / (mol K)

from .._Exceptions import AnalysisError as _AnalysisError
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from ..MD._md import _find_md_engines
from ._restraint import Restraint as _Restraint
from .._SireWrappers import System as _System
from ..Trajectory._trajectory import Trajectory as _Trajectory
from ..Types import Length as _Length
from ..Types import Temperature as _Temperature
from ..Types import Length as _Length
from .. import Units as _Units
from ..Units.Length import angstrom as _angstrom
from ..Units.Angle import radian as _radian
from ..Units.Angle import degree as _degree
from ..Units.Angle import radian as _radian
from ..Units.Energy import kcal_per_mol as _kcal_per_mol
from ..Units.Length import angstrom as _angstrom
from .. import _gmx_exe
from .. import _is_notebook
from .. import Process as _Process
from .. import Protocol as _Protocol
from .. import Units as _Units

if _is_notebook:
    from tqdm.notebook import tqdm as _tqdm
else:
    from tqdm import tqdm as _tqdm

from .._Utils import _try_import, _have_imported, WorkDir as _WorkDir
from .... import _isVerbose

from ..MD._md import _find_md_engines

if _is_notebook:
    from IPython.display import FileLink as _FileLink
    from tqdm.notebook import tqdm as _tqdm
else:
    from tqdm import tqdm as _tqdm


_mda = _try_import("MDAnalysis")

if _have_imported(_mda):
    from MDAnalysis.analysis.distances import dist as _dist
    from MDAnalysis.lib.distances import calc_angles as _calc_angles
    from MDAnalysis.lib.distances import calc_dihedrals as _calc_dihedrals


_MDRestraintsGenerator = _try_import(
    "MDRestraintsGenerator",
    install_command="pip install MDRestraintsGenerator",
)
if _have_imported(_MDRestraintsGenerator):
    from Bio import BiopythonDeprecationWarning

    _warnings.simplefilter("ignore", BiopythonDeprecationWarning)
    from MDRestraintsGenerator import search as _search
    from MDRestraintsGenerator.restraints import (
        FindBoreschRestraint as _FindBoreschRestraint,
    )

is_MDRestraintsGenerator = _have_imported(_MDRestraintsGenerator)


class RestraintSearch:
    """
    Class for running unrestrained simulations from which receptor-ligand
    restraints are selected.
    """

    # Create a list of supported molecular dynamics engines.
    _engines = ["GROMACS", "SOMD"]

    # Create a list of supported restraint types.
    _restraint_types = ["boresch", "multiple_distance"]

    def __init__(
        self,
        system,
        protocol=None,
        work_dir=None,
        engine=None,
        gpu_support=False,
        ignore_warnings=False,
        show_errors=True,
        extra_options=None,
        extra_lines=None,
        property_map={},
        **kwargs,
    ):
        """
        Constructor.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system for the ABFE simulation. This must contain
            a single decoupled molecule and is assumed to have already
            been equilibrated.

        protocol : :class:`Protocol.Production <BioSimSpace.Protocol.Production>`, \
            The simulation production protocol.

        work_dir : str
            The working directory for the ABFE restraint generation
            simulation.

        engine : str
            The molecular dynamics engine used to run the simulation. Available
            options are "AMBER", "GROMACS", or "SOMD". If this argument is omitted
            then BioSimSpace will choose an appropriate engine for you.

        gpu_support : bool
            Whether the engine must have GPU support.

        ignore_warnings : bool
            Whether to ignore warnings when generating the binary run file.
            This option is specific to GROMACS and will be ignored when a
            different molecular dynamics engine is chosen.

        show_errors : bool
            Whether to show warning/error messages when generating the binary
            run file. This option is specific to GROMACS and will be ignored
            when a different molecular dynamics engine is chosen.

        extra_options : dict
            A dictionary containing extra options. Overrides the ones generated from the protocol.

        extra_lines : list
            A list of extra lines to be put at the end of the script.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        kwargs :
            Keyword arguments to be passed to the BSS.Process.
        """

        # Validate the input.
        if not _have_imported(_mda):
            raise _MissingSoftwareError(
                "Cannot perform a RestraintSearch because MDAnalysis is "
                "not installed!"
            )

        if not isinstance(system, _System):
            raise TypeError(
                "'system' must be of type 'BioSimSpace._SireWrappers.System'"
            )
        else:
            # Store a copy of solvated system.
            self._system = system.copy()

        if protocol is not None:
            if isinstance(protocol, _Protocol.Production):
                self._protocol = protocol
            else:
                raise TypeError(
                    "'protocol' must be of type 'BioSimSpace.Protocol.Production'"
                )
        else:
            self._protocol = _Protocol.Production()

        if extra_options is None:
            self._extra_options = {}
        elif not isinstance(extra_options, dict):
            raise ValueError("'extra_options' should be a dict.")
        else:
            self._extra_options = extra_options

        if extra_lines is None:
            self._extra_lines = []
        elif not isinstance(extra_lines, list):
            raise ValueError("'extra_lines' should be a list.")
        else:
            self._extra_lines = extra_lines

        # Create the working directory.
        self._work_dir = _WorkDir(work_dir)

        # There must be a single molecule to be decoupled (or annihilated).
        if system.nDecoupledMolecules() != 1:
            raise ValueError(
                "The system must contain a single molecule to be decoupled! "
                "Use the 'BioSimSpace.Align.Decouple' function to mark a molecule"
                " to be decoupled."
            )

        # Validate the user specified molecular dynamics engine.
        self._exe = None
        if engine is not None:
            if not isinstance(engine, str):
                raise TypeError("'engine' must be of type 'str'.")

            # Strip whitespace from engine and convert to upper case.
            engine = engine.replace(" ", "").upper()

            # Check that the engine is supported.
            if engine not in self._engines:
                raise ValueError(
                    "Unsupported molecular dynamics engine '%s'. "
                    "Supported engines are: %r." % ", ".join(self._engines)
                )

            # Make sure GROMACS is installed if GROMACS engine is selected.
            if engine == "GROMACS":
                if _gmx_exe is None:
                    raise _MissingSoftwareError(
                        "Cannot use GROMACS engine as GROMACS is not installed!"
                    )
                self._exe = _gmx_exe

            elif engine == "AMBER":
                # Find a molecular dynamics engine and executable.
                engines, exes = _find_md_engines(system, protocol, engine, gpu_support)
                if not exes:
                    raise _MissingSoftwareError(
                        "Cannot use AMBER engine as AMBER is not installed!"
                    )
                elif len(exes) > 1:
                    _warnings.warn(
                        f"Multiple AMBER engines were found. Proceeding with {exes[0]}..."
                    )
                self._exe = exes[0]

        else:
            # Use SOMD as a default.
            engine = "SOMD"

        # Set the engine.
        self._engine = engine

        if not isinstance(ignore_warnings, bool):
            raise ValueError("'ignore_warnings' must be of type 'bool'.")
        self._ignore_warnings = ignore_warnings

        if not isinstance(show_errors, bool):
            raise ValueError("'show_errors' must be of type 'bool'.")
        self._show_errors = show_errors

        # Check that the map is valid.
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'.")
        self._property_map = property_map

        # Create fake instance methods for 'analyse'. These
        # pass instance data through to the staticmethod versions.
        self.analyse = self._analyse

        # Initialise the process.
        self._initialise_process(self._system, gpu_support, **kwargs)

    def start(self):
        """Start the simulation."""
        self._process.start()

    def wait(self):
        """Wait for the simulation to finish."""
        self._process.wait()

    def kill(self):
        """Kill the process."""
        self._process.kill()

    def isRunning(self):
        """Check if the process is running."""
        return self._process.isRunning()

    def workDir(self):
        """
        Return the working directory.

        Returns
        -------

        work_dir : str
            The path of the working directory.
        """
        return self._work_dir

    def _analyse(
        self,
        restraint_type="Boresch",
        method="MDRestraintsGenerator",
        append_to_ligand_selection="",
        receptor_selection_str="protein and name CA C N",
        force_constant=None,
        cutoff=10 * _angstrom,
        restraint_idx=0,
        block="AUTO",
    ):
        """
        Analyse trajectory and select restraints which best mimic strongest
        receptor-ligand interactions.

         Parameters
         ----------

         restraint_type: str
             The type of restraints to select, from "Boresch" or "multiple_distance".
             Default is 'Boresch'.

         method: str
             The method to use to derive the restraints. 'BSS' or 'MDRestraintsGenerator'.
             BSS uses the native BioSimSpace derivation.

         append_to_ligand_selection: str
             Appends the supplied string to the default atom selection which chooses
             the atoms in the ligand to consider as potential anchor points. The default
             atom selection is f'resname {ligand_resname} and not name H*'. Uses the
             mdanalysis atom selection language. For example, 'not name O*' will result
             in an atom selection of f'resname {ligand_resname} and not name H* and not
             name O*'.

         receptor_selection_str: str
             The selection string for the atoms in the receptor to consider
             as potential anchor points. The default atom selection is
             'protein and name CA C N'. Uses the mdanalysis atom selection
             language.

         force_constant: BioSimSpace.Types.Energy / BioSimSpace.Types.Area
             The force constant to use for all restraints. For angles, the units of
             area will be converted to A-2 and exchanged for rad-2. If None,
             the default force constants are used, which are 10 kcal mol-1 A-2 [rad-2] when
             method == "MDRestraintsGenerator", or fit to fluctuations observed during
             the simulation is method == "BSS".

        cutoff: BioSimSpace.Types.Length
            The greatest distance between ligand and receptor anchor atoms.
            Only affects behaviour when method == "BSS" Receptor anchors
            further than cutoff Angstroms from the closest ligand anchors will not
            be included in the search for potential anchor points.

         restraint_idx: int
             The index of the restraint from a list of candidate restraints ordered by
             suitability. restraint_idx != 0 is only valid if method == 'BSS'.

         block : bool
             Whether to block until the process has finished running.

         Returns
         -------

         restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
             The restraints of `restraint_type` which best mimic the strongest receptor-ligand
             interactions.
        """
        # Wait for the process to finish.
        if block is True or block == "AUTO":
            self.wait()

        # Return the result of calling the staticmethod, passing in the working
        # directory of this object.
        return RestraintSearch.analyse(
            str(self._work_dir),
            self._system,
            self._process.getTrajectory(),
            self._protocol.getTemperature(),
            restraint_type=restraint_type,
            method=method,
            append_to_ligand_selection=append_to_ligand_selection,
            receptor_selection_str=receptor_selection_str,
            cutoff=cutoff,
            force_constant=force_constant,
            restraint_idx=restraint_idx,
        )

    def _initialise_process(self, system, gpu_support, **kwargs):
        """
        Internal helper function to initialise the process.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.

        gpu_support : bool
            Whether the engine must have GPU support.

        kwargs :
            Keyword arguments to be passed to the BSS.Process.
        """

        # Convert to an appropriate AMBER topology. (Required by SOMD for its
        # FEP setup.)
        if self._engine == "SOMD":
            system._set_water_topology("AMBER", property_map=self._property_map)

        # Setup all of the simulation processes.

        # SOMD.
        if self._engine == "SOMD":
            # Check for GPU support.
            if "CUDA_VISIBLE_DEVICES" in _os.environ:
                platform = "CUDA"
            else:
                if gpu_support:
                    raise ValueError(
                        "gpu_support cannot be True if CUDA_VISIBLE_DEVICES is not set."
                    )
                platform = "CPU"

            self._process = _Process.Somd(
                system,
                self._protocol,
                platform=platform,
                work_dir=self._work_dir,
                property_map=self._property_map,
                extra_options=self._extra_options,
                extra_lines=self._extra_lines,
                **kwargs,
            )

        # GROMACS.
        elif self._engine == "GROMACS":
            self._process = _Process.Gromacs(
                system,
                self._protocol,
                work_dir=self._work_dir,
                ignore_warnings=self._ignore_warnings,
                show_errors=self._show_errors,
                extra_options=self._extra_options,
                extra_lines=self._extra_lines,
                **kwargs,
            )
            if gpu_support:
                self._process.setArg("-update", "gpu")

        # AMBER.
        elif self._engine == "AMBER":
            self._process = _Process.Amber(
                system,
                self._protocol,
                exe=self._exe,
                work_dir=self._work_dir,
                extra_options=self._extra_options,
                extra_lines=self._extra_lines,
                **kwargs,
            )

    @staticmethod
    def analyse(
        work_dir,
        system,
        traj,
        temperature,
        restraint_type="Boresch",
        method="MDRestraintsGenerator",
        append_to_ligand_selection="",
        receptor_selection_str="protein and name CA C N",
        force_constant=None,
        cutoff=10 * _angstrom,
        restraint_idx=0,
    ):
        """
        Analyse existing trajectory from a simulation working directory and
        select restraints which best mimic the strongest receptor-ligand
        interactions.

           Parameters
           ----------

           work_dir : str
               The working directory for the simulation.

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system for the ABFE restraint calculation. This
               must contain a single decoupled molecule and is assumed to have
               already been equilibrated.

           traj : :class:`System <BioSimSpace.Trajectory._trajectory.Trajectory>`
               The trajectory for analysis.

           temperature : :class:`System <BioSimSpace.Types.Temperature>`
               The temperature of the system

           restraint_type : str
               The type of restraints to select, from 'Boresch' or 'multiple_distance'.
               Default is "Boresch".

           method : str
                The method to use to derive the restraints. 'BSS' or 'MDRestraintsGenerator'.
                BSS uses the native BioSimSpace derivation.

           append_to_ligand_selection : str
               Appends the supplied string to the default atom selection which chooses
               the atoms in the ligand to consider as potential anchor points. The default
               atom selection is f'resname {ligand_resname} and not name H*'. Uses the
               mdanalysis atom selection language. For example, 'not name O*' will result
               in an atom selection of f'resname {ligand_resname} and not name H* and not
               name O*'.

           receptor_selection_str : str
               The selection string for the atoms in the receptor to consider
               as potential anchor points. The default atom selection is
               'protein and name CA C N'. Uses the mdanalysis atom selection
               language.

           force_constant: BioSimSpace.Types.Energy / BioSimSpace.Types.Area
               The force constant to use for all restraints. For angles, the units of
               area will be converted to A-2 and exchanged for rad-2. If None,
               the default force constants are used, which are 10 kcal mol-1 A-2 [rad-2] when
               method == "MDRestraintsGenerator", or fit to fluctuations observed during
               the simulation is method == "BSS".

           cutoff: BioSimSpace.Types.Length,
               The greatest distance between ligand and receptor anchor atoms.
               Only affects behaviour when method == "BSS" Receptor anchors
               further than cutoff Angstroms from the closest ligand anchors will not
               be included in the search for potential anchor points.

            restraint_idx: int
                The index of the restraint from a list of candidate restraints ordered by
                suitability. restraint_idx != 0 is only valid if method == 'BSS'.

           Returns
           -------

           restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
               The restraints of `restraint_type` which best mimic the strongest receptor-ligand
               interactions.
        """
        _supported_methods = {
            "boresch": ["BSS", "MDRestraintsGenerator"],
            "multiple_distance": ["BSS"],
        }

        # Check all inputs

        if not isinstance(work_dir, str):
            raise TypeError(f"work_dir: {work_dir} must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError(f"work_dir: {work_dir} doesn't exist!")

        if not isinstance(system, _System):
            raise TypeError(
                f"system {type(system)} must be of type 'BioSimSpace._SireWrappers.System'"
            )
        else:
            # Store a copy of solvated system.
            _system = system.copy()

        if not isinstance(traj, _Trajectory):
            raise TypeError(
                f"traj {type(traj)} must be of type 'BioSimSpace.Trajectory._trajectory.Trajectory'"
            )

        if not isinstance(temperature, _Temperature):
            raise ValueError(
                f"temperature {type(temperature)} must be of type 'BioSimSpace.Types.Temperature'"
            )

        if not isinstance(restraint_type, str):
            raise TypeError(
                f"restraint_type {type(restraint_type)} must be of type 'str'."
            )

        if not restraint_type.lower() in RestraintSearch._restraint_types:
            raise NotImplementedError(
                f"Restraint type {restraint_type} is not implemented. "
                f"Please choose from {RestraintSearch._restraint_types}"
            )

        if not isinstance(method, str):
            raise TypeError(f"method {type(method)} must be of type 'str'.")
        if not method.lower() in ["mdrestraintsgenerator", "bss"]:
            raise NotImplementedError(
                "Deriving restraints using 'MDRestraintsGenerator'"
                "or 'BSS' are the only options implemented."
            )
        # Check that the selected method is supported for the selected restraint type.
        if not method in _supported_methods[restraint_type.lower()]:
            raise NotImplementedError(
                f"Method {method} is not supported for restraint type {restraint_type}. "
                f"Please choose from {_supported_methods[restraint_type.lower()]}"
            )

        if method.lower() == "mdrestraintsgenerator":
            if not is_MDRestraintsGenerator:
                raise ValueError(
                    "Please install MDRestraintsGenerator to search for restraints with it. "
                    "Alternatively, use the 'BSS' method."
                )

        if not isinstance(append_to_ligand_selection, str):
            raise TypeError(
                f"append_to_lig_selection {type(append_to_ligand_selection)} must be of type 'str'."
            )

        if not isinstance(receptor_selection_str, str):
            raise TypeError(
                f"append_to_recept_selection {type(receptor_selection_str)} must be of type 'str'."
            )

        if not isinstance(cutoff, _Length):
            raise TypeError(
                f"cutoff {type(cutoff)} must be of type 'BioSimSpace.Types.Length.'"
            )

        if force_constant:
            dim = force_constant.dimensions()
            if dim != (1, 0, -2, 0, 0, -1, 0):
                raise ValueError(
                    "force_constant must be of type "
                    "'BioSimSpace.Types.Energy'/'BioSimSpace.Types.Length^2'"
                    " or NoneType"
                )

        # There must be a single molecule to be decoupled (or annihilated).
        if system.nDecoupledMolecules() != 1:
            raise ValueError(
                "The system must contain a single molecule to be decoupled! "
                "Use the 'BioSimSpace.Align.Decouple' function to mark a molecule"
                " to be decoupled."
            )

        # Extract restraints from simulation
        # Get mdanalysis universe object
        u = traj.getTrajectory(format="mdanalysis")

        # Find decoupled molecule and use it to create ligand selection
        decoupled_mol = _system.getDecoupledMolecules()[0]
        decoupled_resname = decoupled_mol.getResidues()[0].name()

        ligand_selection_str = f"((resname {decoupled_resname}) and (not name H*))"
        if append_to_ligand_selection:
            ligand_selection_str += " and "
            ligand_selection_str += append_to_ligand_selection

        if restraint_type.lower() == "boresch":
            return RestraintSearch._boresch_restraint(
                u,
                system,
                temperature,
                ligand_selection_str,
                receptor_selection_str,
                method,
                work_dir,
                force_constant,
                cutoff,
                restraint_idx=restraint_idx,
            )

        if restraint_type.lower() == "multiple_distance":
            return RestraintSearch._multiple_distance_restraint(
                u,
                system,
                temperature,
                ligand_selection_str,
                receptor_selection_str,
                method,
                work_dir,
                cutoff,
            )

    @staticmethod
    def _boresch_restraint(
        u,
        system,
        temperature,
        ligand_selection_str,
        receptor_selection_str,
        method,
        work_dir,
        force_constant,
        cutoff,
        restraint_idx=0,
    ):
        """
        Generate the Boresch Restraint.

        Parameters
        ----------

        u : MDAnalysis.Universe
            The trajectory for the ABFE restraint calculation as a
            MDAnalysis.Universe object.

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system for the ABFE restraint calculation. This
            must contain a single decoupled molecule and is assumed to have
            already been equilibrated.

        temperature : :class:`System <BioSimSpace.Types.Temperature>`
            The temperature of the system

        ligand_selection_str : str
            The selection string for the atoms in the ligand to consider
            as potential anchor points.

        receptor_selection_str : str
            The selection string for the atoms in the receptor to consider
            as potential anchor points. Uses the mdanalysis atom selection
            language.

        method : str
            The method to use to derive the restraints. 'BSS' or 'MDRestraintsGenerator'.
            BSS uses the native BioSimSpace derivation.

        work_dir : str
            The working directory for the simulation.

        force_constant: BioSimSpace.Types.Energy / BioSimSpace.Types.Area
            The force constant to use for all restraints. For angles, the units of
            area will be converted to A-2 and exchanged for rad-2. If None,
            the default force constants are used, which are 10 kcal mol-1 A-2 [rad-2] when
            method == "MDRestraintsGenerator", or fit to fluctuations observed during
            the simulation is method == "BSS".

        cutoff: BioSimSpace.Types.Length
            The greatest distance between ligand and receptor anchor atoms.
            Only affects behaviour when method == "BSS" Receptor anchors
            further than cutoff Angstroms from the closest ligand anchors will not
            be included in the search for potential anchor points.

        restraint_idx: int
            The index of the restraint from a list of candidate restraints ordered by
            suitability. restraint_idx != 0 is only valid if method == 'BSS'.

        cutoff : BioSimSpace.Types.Length
            The greatest distance between ligand and receptor anchor atoms, in
            Angstrom. Receptor anchors further than cutoff Angstroms from the closest
            ligand anchors will not be included in the search for potential anchor points.

        Returns
        -------

        restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
            The restraints of `restraint_type` which best mimic the strongest receptor-ligand
            interactions.
        """
        if method == "MDRestraintsGenerator":
            if not is_MDRestraintsGenerator:
                raise ValueError(
                    "Please install MDRestraintsGenerator to search for restraints with it. "
                    "Alternatively, use the 'BSS' method."
                )
            if restraint_idx != 0:
                raise ValueError("restraint_idx must be 0 for MDRestraintsGenerator.")
            else:
                return RestraintSearch._boresch_restraint_MDRestraintsGenerator(
                    u,
                    system,
                    temperature,
                    ligand_selection_str,
                    receptor_selection_str,
                    force_constant,
                    work_dir,
                )

        elif method == "BSS":
            return RestraintSearch._boresch_restraint_BSS(
                u,
                system,
                temperature,
                ligand_selection_str,
                receptor_selection_str,
                work_dir,
                force_constant,
                cutoff,
                restraint_idx=restraint_idx,
            )

    @staticmethod
    def _boresch_restraint_MDRestraintsGenerator(
        u,
        system,
        temperature,
        ligand_selection_str,
        receptor_selection_str,
        force_constant,
        work_dir,
    ):
        """
        Generate the Boresch Restraint using MDRestraintsGenerator.

        Parameters
        ----------

        u : MDAnalysis.Universe
            The trajectory for the ABFE restraint calculation as a
            MDAnalysis.Universe object.

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system for the ABFE restraint calculation. This
            must contain a single decoupled molecule and is assumed to have
            already been equilibrated.

        temperature : :class:`System <BioSimSpace.Types.Temperature>`
            The temperature of the system

        ligand_selection_str: str
            The selection string for the atoms in the ligand to consider
            as potential anchor points.

        receptor_selection_str: str
            The selection string for the protein in the ligand to consider
            as potential anchor points.

        force_constant: BioSimSpace.Types.Energy / BioSimSpace.Types.Area
            The force constant to use for all restraints. For angles, the units of
            area will be converted to A-2 and exchanged for rad-2. If None,
            the default force constants are used, which are 10 kcal mol-1 A-2 [rad-2] when
            method == "MDRestraintsGenerator", or fit to fluctuations observed during
            the simulation is method == "BSS".

        work_dir : str
            The working directory for the simulation.

        Returns
        -------

        restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
            The restraints of `restraint_type` which best mimic the strongest receptor-ligand
            interactions.
        """
        print(
            "Using MDRestraintsGenerator to generate Boresch restraints. If you publish "
            "any results using this method, please cite: 10.5281/zenodo.4570556 "
            "and https://doi.org/10.1038/s42004-022-00721-4"
        )

        ligand_atoms = _search.find_ligand_atoms(
            u, l_selection=ligand_selection_str, p_align=receptor_selection_str
        )
        # find protein atoms
        atom_set = []
        for l_atoms in ligand_atoms:
            psearch = _search.FindHostAtoms(
                u, l_atoms[0], p_selection=receptor_selection_str
            )
            psearch.run()
            atom_set.extend([(l_atoms, p) for p in psearch.host_atoms])
        # Create the boresch finder analysis object
        boresch = _FindBoreschRestraint(u, atom_set)
        # Run the restraint analysis
        boresch.run()

        # Save the analysis results
        boresch.restraint.plot(path=work_dir)
        # Write out the intermolecular section to a topology
        boresch.restraint.write(path=work_dir)
        dG_off = boresch.restraint.standard_state()
        with open(f"{work_dir}/dG_off.dat", "w") as writer:
            writer.write(f"{dG_off}")

        # This just shows how is the data being recorded, so the
        # index gets reassigned multiple times.
        # https://github.com/IAlibay/MDRestraintsGenerator/blob/master/MDRestraintsGenerator/datatypes.py#L810-L825
        l1_idx, r1_idx = boresch.restraint.bond.atomgroup.atoms.ix
        l2_idx, l1_idx, r1_idx = boresch.restraint.angles[0].atomgroup.atoms.ix
        l1_idx, r1_idx, r2_idx = boresch.restraint.angles[1].atomgroup.atoms.ix
        l3_idx, l2_idx, l1_idx, r1_idx = boresch.restraint.dihedrals[
            0
        ].atomgroup.atoms.ix
        l2_idx, l1_idx, r1_idx, r2_idx = boresch.restraint.dihedrals[
            1
        ].atomgroup.atoms.ix
        l1_idx, r1_idx, r2_idx, r3_idx = boresch.restraint.dihedrals[
            2
        ].atomgroup.atoms.ix

        # Select force constants
        if force_constant:
            k_dist = force_constant
            k_ang = (
                (force_constant / (_kcal_per_mol / (_angstrom**2)))
                * _kcal_per_mol
                / (_radian**2)
            )
        else:
            k_dist = 10 * _kcal_per_mol / (_angstrom**2)
            k_ang = 10 * _kcal_per_mol / (_radian**2)

        # The index of the best frame
        index = boresch.restraint.min_frame
        # r1-l1 (r0, kr)
        r0 = boresch.restraint.bond.values[index] * _angstrom
        kr = k_dist
        # r2-r1-l1 (thetaA0, kthetaA)
        thetaA0 = boresch.restraint.angles[1].values[index] * _degree
        kthetaA = k_ang
        # r1-l1-l2 (thetaB0, kthetaB)
        thetaB0 = boresch.restraint.angles[0].values[index] * _degree
        kthetaB = k_ang
        # r3-r2-r1-l1 (phiA0, kphiA)
        phiA0 = boresch.restraint.dihedrals[2].values[index] * _degree
        kphiA = k_ang
        # r2-r1-l1-l2 (phiB0, kphiB)
        phiB0 = boresch.restraint.dihedrals[1].values[index] * _degree
        kphiB = k_ang
        # r1-l1-l2-l3 (phiC0, kphiC)
        phiC0 = boresch.restraint.dihedrals[0].values[index] * _degree
        kphiC = k_ang

        restraint_dict = {
            # The default index is in the format of numpy.int64
            # So need to convert to int
            "anchor_points": {
                "r1": system.getAtom(int(r1_idx)),
                "r2": system.getAtom(int(r2_idx)),
                "r3": system.getAtom(int(r3_idx)),
                "l1": system.getAtom(int(l1_idx)),
                "l2": system.getAtom(int(l2_idx)),
                "l3": system.getAtom(int(l3_idx)),
            },
            "equilibrium_values": {
                "r0": r0,
                "thetaA0": thetaA0,
                "thetaB0": thetaB0,
                "phiA0": phiA0,
                "phiB0": phiB0,
                "phiC0": phiC0,
            },
            "force_constants": {
                "kr": kr,
                "kthetaA": kthetaA,
                "kthetaB": kthetaB,
                "kphiA": kphiA,
                "kphiB": kphiB,
                "kphiC": kphiC,
            },
        }
        # TODO: extract the best frame and feed it into Restraint
        # Waiting for the BSS to fix the getFrames
        # best_frame = traj.getFrames(index)
        best_frame = system
        restraint = _Restraint(
            best_frame, restraint_dict, temperature, restraint_type="Boresch"
        )
        return restraint

    @staticmethod
    def _findOrderedPairs(u, ligand_selection_str, receptor_selection_str, cutoff):
        """
        Return a list of receptor-ligand anchor atoms pairs in the form
        (lig atom index, receptor atom index), where the pairs are ordered
        from low to high variance of distance over the trajectory.

        Parameters
        ----------

        u : MDAnalysis.Universe
            The trajectory for the ABFE restraint calculation as a
            MDAnalysis.Universe object.

        ligand_selection_str: str
            The selection string for the atoms in the ligand to consider
            as potential anchor points.

        receptor_selection_str: str
            The selection string for the protein in the ligand to consider
            as potential anchor points.

        cutoff: BioSimSpace.Types.Length
            The greatest distance between ligand and receptor anchor atoms.
            Only affects behaviour when method == "BSS" Receptor anchors
            further than cutoff Angstroms from the closest ligand anchors will not
            be included in the search for potential anchor points.

        Returns
        -------

        pairs_ordered_sd : list of tuples
            List of receptor-ligand atom pairs ordered by increasing variance of distance over
            the trajectory.
        """

        lig_selection = u.select_atoms(ligand_selection_str)
        pair_variance_dict = {}

        # Get all receptor atoms within specified distance of cutoff
        for lig_atom in lig_selection:
            for prot_atom in u.select_atoms(
                f"{receptor_selection_str} and (around {cutoff / _angstrom} index {lig_atom.index})"
            ):
                pair_variance_dict[(lig_atom.index, prot_atom.index)] = {}
                pair_variance_dict[(lig_atom.index, prot_atom.index)]["dists"] = []

        # Compute Average Distance and SD
        for frame in _tqdm(
            u.trajectory, desc="Searching for low variance pairs. Frame no: "
        ):
            for lig_atom_index, prot_atom_index in pair_variance_dict.keys():
                distance = _dist(
                    _mda.AtomGroup([u.atoms[lig_atom_index]]),
                    _mda.AtomGroup([u.atoms[prot_atom_index]]),
                    box=frame.dimensions,
                )[2][0]
                pair_variance_dict[(lig_atom_index, prot_atom_index)]["dists"].append(
                    distance
                )

        # change lists to numpy arrays
        for pair in pair_variance_dict.keys():
            pair_variance_dict[pair]["dists"] = _np.array(
                pair_variance_dict[pair]["dists"]
            )

        # calculate SD
        for pair in pair_variance_dict.keys():
            pair_variance_dict[pair]["sd"] = pair_variance_dict[pair]["dists"].std()

        # get n pairs with lowest SD
        pairs_ordered_sd = []
        for item in sorted(pair_variance_dict.items(), key=lambda item: item[1]["sd"]):
            pairs_ordered_sd.append(item[0])

        if len(pairs_ordered_sd) == 0:
            raise _AnalysisError(
                "No receptor-ligand atom pairs found within specified cutoff."
            )

        return pairs_ordered_sd

    @staticmethod
    def _getDistance(idx1, idx2, u):
        """
        Get the distance between two atoms in a universe.

        Parameters
        ----------
        idx1 : int
            Index of the first atom
        idx2 : int
            Index of the second atom
        u : MDAnalysis.Universe
            The MDA universe containing the atoms and
            trajectory.

        Returns
        -------
        distance : float
            The distance between the two atoms in Angstroms.
        """
        distance = _dist(
            _mda.AtomGroup([u.atoms[idx1]]),
            _mda.AtomGroup([u.atoms[idx2]]),
            box=u.dimensions,
        )[2][0]
        return distance

    @staticmethod
    def _boresch_restraint_BSS(
        u,
        system,
        temperature,
        ligand_selection_str,
        receptor_selection_str,
        work_dir,
        force_constant,
        cutoff,
        restraint_idx=0,
    ):
        """
        Generate the Boresch Restraint. This method was inspired by Irfan Alibay's
        MDRestraintsGenerator. Please see:

        https://doi.org/10.5281/zenodo.4570556
        https://doi.org/10.1038/s42004-022-007

        Some of the main differences compared to MDRestraintsGenerator are:

            - Scoring by configurational volume, rather than total variance
            - Setting the force constants based on variances, rather than uniformly
            - Checking the restraints for instabilities based on the energy penalty
              for approaching points of instability, rather than absolute values of
              e.g. angles.

        Parameters
        ----------

        u : MDAnalysis.Universe
            The trajectory for the ABFE restraint calculation as a
            MDAnalysis.Universe object.

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system for the ABFE restraint calculation. This
            must contain a single decoupled molecule and is assumed to have
            already been equilibrated.

        temperature : :class:`System <BioSimSpace.Types.Temperature>`
            The temperature of the system

        ligand_selection_str: str
            The selection string for the atoms in the ligand to consider
            as potential anchor points.

        receptor_selection_str: str
            The selection string for the protein in the ligand to consider
            as potential anchor points.

        work_dir : str
            The working directory for the simulation.

        force_constant: BioSimSpace.Types.Energy / BioSimSpace.Types.Area
            The force constant to use for all restraints. For angles, the units of
            area will be converted to A-2 and exchanged for rad-2. If None,
            the default force constants are used, which are 10 kcal mol-1 A-2 [rad-2] when
            method == "MDRestraintsGenerator", or fit to fluctuations observed during
            the simulation is method == "BSS".

        cutoff: BioSimSpace.Types.Length
            The greatest distance between ligand and receptor anchor atoms.
            Only affects behaviour when method == "BSS" Receptor anchors
            further than cutoff Angstroms from the closest ligand anchors will not
            be included in the search for potential anchor points.

        restraint_idx: int
            The index of the restraint from a list of candidate restraints ordered by
            increasing configurational volume - a lower index gives a stronger restraint.

        Returns
        -------

        restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
            The restraints of `restraint_type` which best mimic the strongest receptor-ligand
            interactions.
        """

        def _getAnchorAts(a1_idx, selection_str, u):
            """
            Takes in index of anchor atom 1 (in either the receptor or ligand)
            and universe and returns list of all three anchor atoms, which are chosen
            to be contiguous and to satisfy the selection string. Only one set of anchor
            points is chosen per index so as to search a wider variery of anchor points
            without dramatically increasing the set of candidate anchor points searched.

            Parameters
            ----------

            a1_idx : int
                Index of the first anchor atom

            selection_str : str
                The selection of atoms from which anchor atoms
                may be selected. Uses MDAnalysis atom selection
                language.

            u : MDAnalysis.Universe
                The trajectory for the ABFE restraint calculation as a
                MDAnalysis.Universe object.

            Returns
            -------

            inds : List of ints
                The indices of all three selected anchor atoms to be used in
                Boresch restraints.
            """
            # Get the atoms bonded to the first anchor point which satisfy
            # selection string
            a1_at = u.atoms[a1_idx]
            bonded_at_a1 = a1_at.bonded_atoms.select_atoms(selection_str)
            if len(bonded_at_a1) == 0:
                raise _AnalysisError(
                    "Could not find anchor points matching search critera"
                )

            # Take the first bonded atom to be the second anchor point
            a2_idx = bonded_at_a1[0].index

            # Try to take the second atom for the list of those bonded to a2
            if len(bonded_at_a1) > 1:
                a3_idx = bonded_at_a1[1].index
            # Otherwise take from atoms bonded to a2
            else:
                bonded_at_a2 = bonded_at_a1[0].bonded_atoms.select_atoms(selection_str)
                bonded_at_a2 -= a1_at  # Ensure we do not select a1 again
                if len(bonded_at_a2) == 0:
                    raise _AnalysisError(
                        "Could not find anchor points matching search critera"
                    )
                else:
                    a3_idx = bonded_at_a2[0].index

            return a1_idx, a2_idx, a3_idx

        def _getAngle(idx1, idx2, idx3, u):
            """
            Get the angle between three atoms in a universe.

            Parameters
            ----------
            idx1 : int
                Index of the first atom
            idx2 : int
                Index of the second atom
            idx3 : int
                Index of the third atom
            u : MDAnalysis.Universe
                The MDA universe containing the atoms and
                trajectory.

            Returns
            -------
            angle : float
                The angle between the three atoms in radians.
            """
            angle = sum(
                u.atoms[idx] for idx in [idx1, idx2, idx3]
            ).angle.value()  # Degrees
            angle = _np.deg2rad(angle)  # Radians
            return angle

        def _getDihedral(idx1, idx2, idx3, idx4, u):
            """
            Get the dihedral angle between four atoms in a universe.

            Parameters
            ----------
            idx1 : int
                Index of the first atom
            idx2 : int
                Index of the second atom
            idx3 : int
                Index of the third atom
            idx4 : int
                Index of the fourth atom
            u : MDAnalysis.Universe
                The MDA universe containing the atoms and
                trajectory.

            Returns
            -------
            dihedral : float
                The dihedral angle between the four atoms in radians.
            """
            dihedral = sum(
                u.atoms[idx] for idx in [idx1, idx2, idx3, idx4]
            ).dihedral.value()  # Degrees
            dihedral = _np.deg2rad(dihedral)  # Radians
            return dihedral

        def _getBoreschDOF(l1, l2, l3, r1, r2, r3, u):
            """Calculate Boresch degrees of freedom from indices of anchor atoms"""
            # Ordering of connection of anchors is r3,r2,r1,l1,l2,l3
            r = RestraintSearch._getDistance(r1, l1, u)
            thetaA = _getAngle(r2, r1, l1, u)
            thetaB = _getAngle(r1, l1, l2, u)
            phiA = _getDihedral(r3, r2, r1, l1, u)
            phiB = _getDihedral(r2, r1, l1, l2, u)
            phiC = _getDihedral(r1, l1, l2, l3, u)
            # Not restrained but distance from collinearity must be checked
            thetaR = _getAngle(r3, r2, r1, u)  # Receptor internal angle
            thetaL = _getAngle(l1, l2, l3, u)  # Ligand internal angle
            return r, thetaA, thetaB, phiA, phiB, phiC, thetaR, thetaL

        def _getConfigVol(equil_vals, force_consts, temp):
            """
            Find the configurational volume accessible to the
            decoupled restrained ligand based on the Boresch restraint.
            Based on part of Eqn 32 of J. Phys. Chem. B 2003, 107, 35, 9535–9551.

            Parameters
            ----------

            equil_vals : dict
                Dictionary of equilibrium values for the Boresch restraint. Must have units
                of Angstrom and radians. Of the form {"r": r0, "thetaA": thetaA0,
                "thetaB": thetaB0, "phiA": phiA0, "phiB": phiB0, "phiC": phiC0}.

            force_consts : dict
                Dictionary of force constants for the Boresch restraint. Force constants
                must have units of kcal /mol. Of the form {"r": kr, "thetaA": kthetaA,
                "thetaB": kthetaB, "phiA": kphiA0, "phiB": kphiB, "phiC": kphiC}.

            temp : float
                The temperature, in K.

            Returns
            -------

            config_vol : float
                The configurational volume accessible to the restrained decoupled ligand,
                in Angstrom^3.
            """
            RT = _k_boltz.value() * temp  # in kcal / mol
            numerator1 = (
                (equil_vals["r"] ** 2)
                * _np.sin(equil_vals["thetaA"])
                * _np.sin(equil_vals["thetaB"])
            )  # Units: A**2
            numerator2 = (2 * _np.pi * RT) ** 3  # Units: (kcal / mol )**3
            denominator = _np.sqrt(
                _np.array([val for val in force_consts.values()]).prod()
            )
            # Units: (kcal / mol)**3 A**-1
            config_vol = numerator1 * numerator2 / denominator  # Units: A**3
            return config_vol

        def _findOrderedBoresch(
            u,
            ligand_selection_str,
            receptor_selection_str,
            pair_list,
            temp,
            force_constant,
            no_pairs=50,
        ):
            """
            Calculate a list of Boresch restraints and associated
            statistics over the trajectory.

            Parameters
            ----------

            u : MDAnalysis.Universe
                The trajectory for the ABFE restraint calculation as a
                MDAnalysis.Universe object.

            ligand_selection_str: str
                The selection string for the atoms in the ligand to consider
                as potential anchor points.

            receptor_selection_str: str
                The selection string for the protein in the ligand to consider
                as potential anchor points.

            pair_list : List of tuples
                List of receptor-ligand atom pairs to be used as the r1 and l1
                anchor points in candidate Boresch restraints.

            temp : float
                The temperature, in K.

            force_constant: BioSimSpace.Types.Energy / BioSimSpace.Types.Area
                The force constant to use for all restraints. For angles, the units of
                area will be converted to A-2 and exchanged for rad-2. If None,
                the default force constants are used, which are 10 kcal mol-1 A-2 [rad-2] when
                method == "MDRestraintsGenerator", or fit to fluctuations observed during
                the simulation is method == "BSS".

            no_pairs : int
                Number of pairs to be used in the calculation. Pairs in pair_list
                after the index supplied will be ignored.

            Returns
            -------

            pairs_ordered_boresch : List of tuples
                List of receptor-ligand atom pairs from pair_list ordered by priority
                the Boresch restraints they generate. pairs_ordered_boresch[0] labels
                the optimal restraint according to the algorithm implemented.

            boresch_dof_data: dict
                Dictionary of statistics for the Boresch restraints obtained over the
                trajectory. Keys are the pair tuples supplied in pair_list.
            """
            boresch_dof_list = [
                "r",
                "thetaA",
                "thetaB",
                "phiA",
                "phiB",
                "phiC",
                "thetaR",
                "thetaL",
            ]  # thetaR and thetaL are the internal
            # angles of the receptor and ligand

            # get values of degrees of freedom for lowest SD pairs across whole trajectory

            boresch_dof_data = {}
            for pair in _tqdm(
                pair_list[:no_pairs],
                desc="Scoring candidate Boresch anchor points. Anchor set no: ",
            ):
                l1_idx, r1_idx = pair
                try:
                    _, l2_idx, l3_idx = _getAnchorAts(l1_idx, ligand_selection_str, u)
                    _, r2_idx, r3_idx = _getAnchorAts(r1_idx, receptor_selection_str, u)
                except (
                    _AnalysisError
                ):  # Failed to find full set of anchor points for this pair
                    continue

                boresch_dof_data[pair] = {}
                boresch_dof_data[pair]["anchor_ats"] = [
                    l1_idx,
                    l2_idx,
                    l3_idx,
                    r1_idx,
                    r2_idx,
                    r3_idx,
                ]

                # Add sub dictionaries for each Boresch degree of freedom
                for dof in boresch_dof_list:
                    boresch_dof_data[pair][dof] = {}
                    boresch_dof_data[pair][dof]["values"] = []

                for i, _ in enumerate(
                    u.trajectory
                ):  # TODO: Use MDA.analysis.base instead?
                    (
                        r,
                        thetaA,
                        thetaB,
                        phiA,
                        phiB,
                        phiC,
                        thetaR,
                        thetaL,
                    ) = _getBoreschDOF(
                        l1_idx, l2_idx, l3_idx, r1_idx, r2_idx, r3_idx, u
                    )
                    boresch_dof_data[pair]["r"]["values"].append(r)
                    boresch_dof_data[pair]["thetaA"]["values"].append(thetaA)
                    boresch_dof_data[pair]["thetaB"]["values"].append(thetaB)
                    boresch_dof_data[pair]["phiA"]["values"].append(phiA)
                    boresch_dof_data[pair]["phiB"]["values"].append(phiB)
                    boresch_dof_data[pair]["phiC"]["values"].append(phiC)
                    boresch_dof_data[pair]["thetaR"]["values"].append(thetaR)
                    boresch_dof_data[pair]["thetaL"]["values"].append(thetaL)

                # Calculate statistics for each Boresch degree of freedom
                for dof in boresch_dof_list:
                    boresch_dof_data[pair][dof]["values"] = _np.array(
                        boresch_dof_data[pair][dof]["values"]
                    )
                    # Check not dihedral
                    if not dof[:3] == "phi":
                        boresch_dof_data[pair][dof]["avg"] = boresch_dof_data[pair][
                            dof
                        ]["values"].mean()
                        boresch_dof_data[pair][dof]["var"] = boresch_dof_data[pair][
                            dof
                        ]["values"].var()
                    # If dihedral, have to calculate circular stats
                    else:
                        circmean = _circmean(
                            boresch_dof_data[pair][dof]["values"],
                            high=_np.pi,
                            low=-_np.pi,
                        )
                        boresch_dof_data[pair][dof]["avg"] = circmean

                        # Cannot use scipy's circvar as later than v 1.8
                        # as this is calculated in the range 0 - 1
                        corrected_values = []
                        for val in boresch_dof_data[pair][dof]["values"]:
                            dtheta = abs(val - circmean)
                            corrected_values.append(min(dtheta, 2 * _np.pi - dtheta))
                        corrected_values = _np.array(corrected_values)
                        boresch_dof_data[pair][dof]["var"] = corrected_values.var()

                    # Assume Gaussian distributions and calculate force constants for harmonic potentials
                    # so as to reproduce these distributions at 298 K
                    boresch_dof_data[pair][dof]["k"] = (
                        _k_boltz.value() * temp / (boresch_dof_data[pair][dof]["var"])
                    )  # Force constants in kcal mol-1 A-2 [rad-2]

                # Calculate the configurational volume accessible based on each restraint
                equil_vals = {
                    dof: boresch_dof_data[pair][dof]["avg"] for dof in boresch_dof_list
                }
                force_consts = {
                    dof: boresch_dof_data[pair][dof]["k"] for dof in boresch_dof_list
                }
                boresch_dof_data[pair]["config_vol"] = _getConfigVol(
                    equil_vals, force_consts, temp
                )

                # Now, after we've used the fluctuation-derived force constants to calculate the
                # configurational volume, set to user-supplied value if specified.
                if force_constant:
                    k = force_constant / (
                        _kcal_per_mol / (_angstrom**2)
                    )  # Convert to kcal mol-1 A-2,
                    # and use same value for angle force
                    # constants in kcal mol-1 rad-2
                    boresch_dof_data[pair]["r"]["k"] = k
                    for restrained_angle in [
                        "thetaA",
                        "thetaB",
                        "phiA",
                        "phiB",
                        "phiC",
                    ]:  # Do not change thetaR, thetaL
                        boresch_dof_data[pair][restrained_angle]["k"] = k

            # Order pairs according to configurational volume - smaller volume indicates stronger
            # restraints mimicking stronger native interactions
            pairs_ordered_boresch_var = []
            for item in sorted(
                boresch_dof_data.items(), key=lambda item: item[1]["config_vol"]
            ):
                pairs_ordered_boresch_var.append(item[0])

            # Filter out force restraints with with 10 kT of collinearity or r = 0
            # Convert 10 kT to angle
            R = _k_boltz.value()  # molar gas constant in kcal mol-1 K-1
            min_stable_dist = lambda k: _np.sqrt(
                (20 * R * temp) / k
            )  # Get the "distance" at which
            # restraint penalty is 10 kT
            pairs_ordered_boresch = []
            for pair in pairs_ordered_boresch_var:
                # Check equil distance
                r0 = boresch_dof_data[pair]["r"]["avg"]
                kr = boresch_dof_data[pair]["r"]["k"]
                stable_distance = r0 > min_stable_dist(kr)
                stable_angle = True
                for angle in ["thetaA", "thetaB", "thetaR", "thetaL"]:
                    # Check equil angle
                    ang0 = boresch_dof_data[pair][f"{angle}"]["avg"]
                    kang = boresch_dof_data[pair][f"{angle}"]["k"]
                    # Check minimum distance to collinearity
                    min_dist = min([abs(ang0 - 0), abs(ang0 - _np.pi)])
                    if min_dist < min_stable_dist(kang):
                        stable_angle = False
                # If no likely instabilities, add pair
                if stable_distance and stable_angle:
                    pairs_ordered_boresch.append(pair)

            if len(pairs_ordered_boresch) == 0:
                raise _AnalysisError(
                    "No candidate sets of Boresch restraints are suitable. Please expand "
                    "search criteria or increase force constants."
                )

            return pairs_ordered_boresch, boresch_dof_data

        def _plotDOF(
            ordered_restraint_labels,
            dof_data,
            restraint_idx=0,
            dof_to_plot=["r", "thetaA", "thetaB", "phiA", "phiB", "phiC"],
        ):
            """
            Plot histograms and variation with time of DOF over a trajectory.

            Parameters
            ----------

            ordered_restraint_labels : list
                List of labels for candidate restraints, which will be used to
                look up the restraints in dof_data. Presumed to be ordered by
                some measure of desirability for use.

            dof_data : dict
                Dictionary of data (obtained from the trajectory) for each restraint
                with a label given in ordered_restraint_labels.

            restraint_idx : int
                Index of the restraint in ordered_restraint_labels to use.

            dof_to_plot : list
                List of DOF to plot.
            """
            # The labels for each DOF to use on the plots
            dof_labels = {
                "r": r"$r$ / $\mathrm{\AA}$",
                "thetaA": r"$\theta_{\mathrm{A}}$ / rad",
                "thetaB": r"$\theta_{\mathrm{B}}$ / rad",
                "phiA": r"$\phi_{\mathrm{A}}$ / rad",
                "phiB": r"$\phi_{\mathrm{B}}$ / rad",
                "phiC": r"$\phi_{\mathrm{C}}$ / rad",
            }

            n_dof = len(dof_to_plot)
            label = ordered_restraint_labels[restraint_idx]

            # Plot histograms
            fig, axs = _plt.subplots(1, n_dof, figsize=(16, 4), dpi=500)
            for i, dof in enumerate(dof_to_plot):
                axs[i].hist(dof_data[label][dof]["values"], bins=10, density=True)
                axs[i].axvline(
                    x=dof_data[label][dof]["avg"],
                    color="r",
                    linestyle="dashed",
                    linewidth=2,
                    label="mean",
                )
                axs[i].set_xlabel(dof_labels[dof])
                axs[i].set_ylabel("Probability")
                if i == n_dof - 1:  # Only add legend to last plot
                    axs[i].legend()
            fig.tight_layout()
            fig.savefig(
                f"{work_dir}/restraint_idx{restraint_idx}_dof_hist.png",
                facecolor="white",
            )
            _plt.close(fig)

            # Plot variation with time to see if there are slow DOF
            fig, axs = _plt.subplots(1, n_dof, figsize=(16, 4), dpi=500)
            for i, dof in enumerate(dof_to_plot):
                axs[i].plot(
                    [x for x in range(len(dof_data[label][dof]["values"]))],
                    dof_data[label][dof]["values"],
                )
                axs[i].axhline(
                    y=dof_data[label][dof]["avg"],
                    color="r",
                    linestyle="dashed",
                    linewidth=2,
                    label="mean",
                )  # No need to add legend as has been
                # added to histograms
                axs[i].set_ylabel(dof_labels[dof])
                axs[i].set_xlabel("Frame No")
            fig.tight_layout()
            fig.savefig(
                f"{work_dir}/restraint_idx{restraint_idx}_dof_time.png",
                facecolor="white",
            )
            _plt.close(fig)

        def _getBoreschRestraint(pair, boresch_dof_data):
            """
            Get the Boresch restraints for a specified pair in a form compatible
            with BSS.

            Parameters
            ----------

            pair : tuple
                The (receptor_idx, ligand_idx) anchor atom pair labelling
                the restraint in boresch_dof_data

            restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
                The restraint defined by boresch_dof_data and labelled by pair.

            """
            anchor_idxs = {
                "l1": boresch_dof_data[pair]["anchor_ats"][0],
                "l2": boresch_dof_data[pair]["anchor_ats"][1],
                "l3": boresch_dof_data[pair]["anchor_ats"][2],
                "r1": boresch_dof_data[pair]["anchor_ats"][3],
                "r2": boresch_dof_data[pair]["anchor_ats"][4],
                "r3": boresch_dof_data[pair]["anchor_ats"][5],
            }

            anchor_ats = {k: system.getAtom(int(v)) for k, v in anchor_idxs.items()}

            # Check we have found all anchors
            if not len(anchor_ats) == 6:
                raise _AnalysisError("Could not find all anchor atoms in system")

            # Get remaining parameters
            r0 = boresch_dof_data[pair]["r"]["avg"]
            thetaA0 = boresch_dof_data[pair]["thetaA"]["avg"]
            thetaB0 = boresch_dof_data[pair]["thetaB"]["avg"]
            phiA0 = boresch_dof_data[pair]["phiA"]["avg"]
            phiB0 = boresch_dof_data[pair]["phiB"]["avg"]
            phiC0 = boresch_dof_data[pair]["phiC"]["avg"]
            kr = boresch_dof_data[pair]["r"]["k"]
            kthetaA = boresch_dof_data[pair]["thetaA"]["k"]
            kthetaB = boresch_dof_data[pair]["thetaB"]["k"]
            kphiA = boresch_dof_data[pair]["phiA"]["k"]
            kphiB = boresch_dof_data[pair]["phiB"]["k"]
            kphiC = boresch_dof_data[pair]["phiC"]["k"]

            restraint_dict = {
                "anchor_points": {
                    "r1": anchor_ats["r1"],
                    "r2": anchor_ats["r2"],
                    "r3": anchor_ats["r3"],
                    "l1": anchor_ats["l1"],
                    "l2": anchor_ats["l2"],
                    "l3": anchor_ats["l3"],
                },
                "equilibrium_values": {
                    "r0": r0 * _Units.Length.angstrom,
                    "thetaA0": thetaA0 * _radian,
                    "thetaB0": thetaB0 * _radian,
                    "phiA0": phiA0 * _radian,
                    "phiB0": phiB0 * _radian,
                    "phiC0": phiC0 * _radian,
                },
                "force_constants": {
                    "kr": kr * _kcal_per_mol / _angstrom**2,
                    "kthetaA": kthetaA * _kcal_per_mol / (_radian * _radian),
                    "kthetaB": kthetaB * _kcal_per_mol / (_radian * _radian),
                    "kphiA": kphiA * _kcal_per_mol / (_radian * _radian),
                    "kphiB": kphiB * _kcal_per_mol / (_radian * _radian),
                    "kphiC": kphiC * _kcal_per_mol / (_radian * _radian),
                },
            }

            restraint = _Restraint(
                system,
                restraint_dict,
                temperature=temperature,
                restraint_type="Boresch",
            )
            return restraint

        # Find pairs with lowest SD
        pairs_ordered_sd = RestraintSearch._findOrderedPairs(
            u, ligand_selection_str, receptor_selection_str, cutoff
        )

        # Convert to Boresch anchors, order by correction, and filter
        pairs_ordered_boresch, boresch_dof_data = _findOrderedBoresch(
            u,
            ligand_selection_str,
            receptor_selection_str,
            pairs_ordered_sd,
            temperature.value(),
            force_constant,
        )

        # Plot
        _plotDOF(pairs_ordered_boresch, boresch_dof_data, restraint_idx=restraint_idx)

        # Convert to BSS compatible dictionary
        restraint = _getBoreschRestraint(
            pairs_ordered_boresch[restraint_idx], boresch_dof_data
        )

        return restraint

    @staticmethod
    def _multiple_distance_restraint(
        u,
        system,
        temperature,
        ligand_selection_str,
        receptor_selection_str,
        method,
        work_dir,
        cutoff,
    ):
        """
        Generate the multiple distance restraints Restraint using the native
        BioSimSpace implementation. This involves:

        1. Finding the pairs of receptor-ligand atoms that have the lowest standard
           deviation of distances
        2. Selecting pairs based on lowest standard deviation of distances and diversity
           of direction. One pair is selected for each selected ligand atom.
        3. Deriving restraint parameters from the selected pairs so that the flat-bottomed
           region encompasses 95 % of the sampled distances. The pair with the lowest
           variance of the distance is selected to be the "permanent" distance restraint.

        Parameters
        ----------

        u : MDAnalysis.Universe
            The trajectory for the ABFE restraint calculation as a
            MDAnalysis.Universe object.

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system for the ABFE restraint calculation. This
            must contain a single decoupled molecule and is assumed to have
            already been equilibrated.

        temperature : :class:`System <BioSimSpace.Types.Temperature>`
            The temperature of the system

        ligand_selection_str : str
            The selection string for the atoms in the ligand to consider
            as potential anchor points.

        receptor_selection_str : str
            The selection string for the atoms in the receptor to consider
            as potential anchor points. Uses the mdanalysis atom selection
            language.

        method : str
            The method to use to derive the restraints. Currently only 'BSS'
            is implemented, which uses the native BioSimSpace derivation.

        work_dir : str
            The working directory for the simulation.

        cutoff: BioSimSpace.Types.Length
            The greatest distance between ligand and receptor anchor atoms.
            Only affects behaviour when method == "BSS" Receptor anchors
            further than cutoff Angstroms from the closest ligand anchors will not
            be included in the search for potential anchor points.

        Returns
        -------

        restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
            The restraints of `restraint_type` which best mimic the strongest receptor-ligand
            interactions.
        """
        if method == "MDRestraintsGenerator":
            raise NotImplementedError(
                "MDRestraintsGenerator is not implemented for multiple distance restraints."
            )

        elif method == "BSS":
            return RestraintSearch._multiple_distance_restraint_BSS(
                u,
                system,
                temperature,
                ligand_selection_str,
                receptor_selection_str,
                work_dir,
                cutoff,
            )

        else:
            raise ValueError("method must be 'MDRestraintsGenerator' or 'BSS'.")

    @staticmethod
    def _multiple_distance_restraint_BSS(
        u,
        system,
        temperature,
        ligand_selection_str,
        receptor_selection_str,
        work_dir,
        cutoff,
    ):
        """
        Generate the multiple distance restraints Restraint using the native
        BioSimSpace implementation. This involves:

        1. Finding the pairs of receptor-ligand atoms that have the lowest standard
           deviation of distances
        2. Selecting pairs based on lowest standard deviation of distances and diversity
           of direction. One pair is selected for each selected ligand atom.
        3. Deriving restraint parameters from the selected pairs so that the flat-bottomed
           region encompasses 95 % of the sampled distances. The pair with the lowest
           variance of the distance is selected to be the "permanent" distance restraint.

        Parameters
        ----------

        u : MDAnalysis.Universe
            The trajectory for the ABFE restraint calculation as a
            MDAnalysis.Universe object.

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system for the ABFE restraint calculation. This
            must contain a single decoupled molecule and is assumed to have
            already been equilibrated.

        temperature : :class:`System <BioSimSpace.Types.Temperature>`
            The temperature of the system

        ligand_selection_str: str
            The selection string for the atoms in the ligand to consider
            as potential anchor points.

        receptor_selection_str: str
            The selection string for the protein in the ligand to consider
            as potential anchor points.

        work_dir : str
            The working directory for the simulation.

        cutoff: BioSimSpace.Types.Length
            The greatest distance between ligand and receptor anchor atoms.
            Only affects behaviour when method == "BSS" Receptor anchors
            further than cutoff Angstroms from the closest ligand anchors will not
            be included in the search for potential anchor points.

        Returns
        -------

        restraint : :class:`Restraint <BioSimSpace.Sandpit.Exscientia.FreeEnergy.Restraint>`
            The restraints of `restraint_type` which best mimic the strongest receptor-ligand
            interactions.
        """

        def _get_norm_vector(frame, pair):
            """
            Get the normalised interatomic vector between two atoms.

            Parameters
            ----------
            frame : MDAnalysis.Timestep
                The frame of the trajectory to use.
            pair : tuple
                The pair of atoms to consider.

            Returns
            -------
            vector : np.ndarray
                The normalised interatomic vector between the two atoms.
            """
            inter_vec = frame.positions[pair[1]] - frame.positions[pair[0]]
            # Check the length of this is not excessive (e.g. due to PBCs)
            norm = _np.linalg.norm(inter_vec)
            if norm > 15:
                raise ValueError(
                    f"The interatomic vector between atoms {pair[0]} and {pair[1]} is "
                    f"excessive ({norm:.2f} A). This is likely due to PBCs. "
                    "Please check your input trajectory."
                )
            else:
                return inter_vec / norm

        def _clusterPairsByDirection(u, pairs_ordered_sd, n_clusters):
            """
            Cluster the pairs of atoms by direction (based on the final structure).

            Parameters
            ----------
            u : MDAnalysis.Universe
                The trajectory for the ABFE restraint calculation as a
                MDAnalysis.Universe object.
            pairs_ordered_sd : list
                The pairs of atoms ordered by SD.
            n_clusters : int
                The number of clusters to use.

            Returns
            -------
            pairs_ordered_clustered : list[tuple]
                The pairs of atoms ordered by SD and clustered by direction.
                The third element of each tuple is the cluster number.
            """
            # Get the normalised interatomic vectors for each pair of atoms,
            # using the final frame of the trajectory.
            vectors_dict = {
                pair: _get_norm_vector(u.trajectory[-1], pair)
                for pair in pairs_ordered_sd
            }
            # Cluster the vectors by direction, using k-means clustering. Predict the cluster
            # indices at the same time.
            cluster_indices = _KMeans(n_clusters=n_clusters, n_init=10).fit_predict(
                list(vectors_dict.values())
            )

            # Create a new list of pairs of atoms, ordered by SD and including the cluster number.
            pairs_ordered_clustered = [
                (pair[0], pair[1], cluster_indices[i])
                for i, pair in enumerate(pairs_ordered_sd)
            ]

            return pairs_ordered_clustered

        def _filterPairsByCluster(pairs_ordered_clustered, selected_ligand_atoms):
            """
            Filter the pairs of atoms by cluster. This is done by checking each pair
            in order of decreasing variance. If the ligand heavy atom has not already
            been used in a restraint and no other restraints from the given directional
            cluster have been selected, the restraint is accepted.

            Parameters
            ----------
            pairs_ordered_clustered : list[tuple]
                The pairs of atoms ordered by SD and clustered by direction.
                The first and second elements of each tuple are the atom indices.
                The third element of each tuple is the cluster number.

            selected_ligand_atoms : list[int]
                The indices of the selected ligand atoms.

            Returns
            -------
            pairs_ordered_clustered_filtered : list[tuple]
                The pairs of atoms ordered by SD and clustered by direction, filtered
                by cluster. The tuples contain the atom indices only.
            """
            filtered_pairs = []
            selected_clusters = []
            available_ligand_atoms = [at.index for at in selected_ligand_atoms]

            # Loop over the pairs of atoms in order of decreasing variance.
            for pair in pairs_ordered_clustered:
                # If the ligand atom has not already been used in a restraint and no
                # other restraints from the given directional cluster have been selected,
                # accept the restraint.
                if (
                    pair[0] in available_ligand_atoms
                    and pair[2] not in selected_clusters
                ):
                    filtered_pairs.append((pair[0], pair[1]))
                    selected_clusters.append(pair[2])
                    available_ligand_atoms.remove(pair[0])

                if len(available_ligand_atoms) == 0:
                    return filtered_pairs

            # If not all ligand atoms have been used in a restraint, return nothing.
            return None

        def _plotDistanceRestraints(distance_dict):
            """
            Plot the variation of the the distances for each pair
            over the trajectory.

            Parameters
            ----------
            distance_dict : dict
                A dictionary containing the distances for each pair, where
                the keys are the pair indices and the values are lists of
                distances in Angstroms.
            """
            n_pairs = len(distance_dict)
            n_columns = 6
            n_rows = int(_np.ceil(n_pairs / n_columns))

            # Plot histograms
            fig, axs = _plt.subplots(
                n_rows, n_columns, figsize=(2.7 * n_columns, 4 * n_rows), dpi=500
            )
            axs = axs.flatten()
            for i, pair in enumerate(distance_dict):
                axs[i].hist(distance_dict[pair], bins=10, density=True)
                axs[i].axvline(
                    x=_np.mean(distance_dict[pair]),
                    color="r",
                    linestyle="dashed",
                    linewidth=2,
                    label="mean",
                )
                axs[i].set_xlabel(f"({pair[0]}, {pair[1]}) Distance / $\mathrm{{\AA}}$")
                axs[i].set_ylabel("Probability")
                if i == n_pairs - 1:  # Only add legend to last plot
                    axs[i].legend()
                # Hide unused axes
                if i >= n_pairs:
                    axs[i].axis("off")

            fig.tight_layout()
            fig.savefig(
                f"{work_dir}/distance_restraints_hist.png",
                facecolor="white",
            )
            _plt.close(fig)

            # Plot variation with time to see if there are slow DOF
            fig, axs = _plt.subplots(
                n_rows, n_columns, figsize=(2.7 * n_columns, 4 * n_rows), dpi=500
            )
            axs = axs.flatten()
            for i, pair in enumerate(distance_dict):
                axs[i].plot(
                    list(range(len(distance_dict[pair]))),
                    distance_dict[pair],
                )
                axs[i].axhline(
                    y=_np.mean(distance_dict[pair]),
                    color="r",
                    linestyle="dashed",
                    linewidth=2,
                    label="mean",
                )  # No need to add legend as has been
                # added to histograms
                axs[i].set_ylabel(f"({pair[0]}, {pair[1]}) Distance / $\mathrm{{\AA}}$")
                axs[i].set_xlabel("Frame Index")
            fig.tight_layout()
            fig.savefig(
                f"{work_dir}/distance_restraints_time.png",
                facecolor="white",
            )
            _plt.close(fig)

        def _getRestraintDict(u, pairs_ordered):
            """
            Generate the multiple distance restraints dictionary given
            a list of pairs of atoms to restrain, ordered by decreasing
            variance of intramolecular distance. The first pair will be
            selected to be the "permanent" distance restraint. The restraints
            are selected so that the flat-bottomed region covers 95 % of
            the distances sampled during the trajectory.

            Parameters
            ----------
            u : MDAnalysis.Universe
                The trajectory for the ABFE restraint calculation as a
                MDAnalysis.Universe object.

            pairs_ordered : list[tuple]
                The pairs of atoms ordered by SD. The first and second
                elements of each tuple are the atom indices.

            Returns
            -------
            restraint_dict : dict
                The multiple distance restraints dictionary.
            """
            # For each pair, get a list of the distances over the trajectory.
            distances = {pair: [] for pair in pairs_ordered}
            for _ in u.trajectory:
                for pair in pairs_ordered:
                    distances[pair].append(
                        RestraintSearch._getDistance(pair[0], pair[1], u)
                    )

            # Plot the distances over the trajectory.
            _plotDistanceRestraints(distances)

            # For each pair, get the restraint parameters and add them to the
            # restraint dictionary.
            restraint_dict = {"distance_restraints": []}
            for i, pair in enumerate(pairs_ordered):
                r1 = system.getAtom(int(pair[1]))
                l1 = system.getAtom(int(pair[0]))
                r0 = _np.mean(distances[pair]) * _angstrom  # Equilibrium distance
                kr = 40 * _kcal_per_mol / (_angstrom**2)  # Default of 40 kcal/mol/A^2
                r_fb = (
                    _np.percentile(distances[pair], 97.5)
                    - _np.percentile(distances[pair], 2.5)
                ) * _angstrom  # Flat-bottomed region
                individual_restraint_dict = {
                    "r1": r1,
                    "l1": l1,
                    "r0": r0,
                    "kr": kr,
                    "r_fb": r_fb,
                }
                # If this is the first pair, add it as the permanent distance restraint.
                if i == 0:
                    restraint_dict["permanent_distance_restraint"] = (
                        individual_restraint_dict
                    )
                else:
                    restraint_dict["distance_restraints"].append(
                        individual_restraint_dict
                    )

            return restraint_dict

        # Find pairs with lowest SD, ordered by SD.
        pairs_ordered_sd = RestraintSearch._findOrderedPairs(
            u, ligand_selection_str, receptor_selection_str, cutoff
        )

        # Cluster pairs according to direction (based on initial structure). Choose
        # the number of clusters to be equal to the number of selected heavy atoms in the ligand.
        selected_ligand_atoms = u.select_atoms(ligand_selection_str).select_atoms(
            "not name H*"
        )
        # Choose initial number of clusters to to the number of selected heavy atoms in the ligand.
        n_clusters_initial = len(selected_ligand_atoms)
        n_clusters_max = 3 * n_clusters_initial
        n_clusters_to_try = [
            round(n) for n in _np.linspace(n_clusters_initial, n_clusters_max, 5)
        ]

        for n_clusters in n_clusters_to_try:
            pairs_ordered_clustered = _clusterPairsByDirection(
                u, pairs_ordered_sd, n_clusters=n_clusters
            )

            # Check each pair in order of reverse variance. If the ligand heavy atom has not
            # already been used in a restraint and no other restraints from the given cluster
            # have been selected, accept the restraint.
            pairs_ordered_clustered_filtered = _filterPairsByCluster(
                pairs_ordered_clustered, selected_ligand_atoms
            )

            if pairs_ordered_clustered_filtered is not None:
                break

            if (
                n_clusters_to_try[-1] == n_clusters
                and pairs_ordered_clustered_filtered is None
            ):
                raise _AnalysisError(
                    "Could not find suitable multiple distance restraints. "
                    "Please try increasing the number of candidate pairs."
                )

        # Get restraint parameters for each pair by selecting the flat-bottomed region
        # to encompass 95 % of sampled configurations. This also plots the data.
        restraint_dict = _getRestraintDict(u, pairs_ordered_clustered_filtered)

        # Create the restraint object and return it.
        restraint = _Restraint(
            system,
            restraint_dict,
            temperature=temperature,
            restraint_type="multiple_distance",
        )

        return restraint
