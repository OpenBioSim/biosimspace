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
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""Functionality for running simulations with GROMACS."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Gromacs"]

import glob as _glob
import os as _os

from .._Utils import _try_import

_pygtail = _try_import("pygtail")
import shutil as _shutil
import subprocess as _subprocess
import timeit as _timeit
import warnings as _warnings
from tempfile import TemporaryDirectory as _TemporaryDirectory
from pathlib import Path as _Path

import numpy as _np

from sire.legacy import Base as _SireBase
from sire.legacy import IO as _SireIO
from sire.legacy import Maths as _SireMaths
from sire.legacy import Units as _SireUnits
from sire.legacy import Vol as _SireVol

from .._Utils import _assert_imported, _have_imported, _try_import

# alchemlyb isn't available on all variants of Python that we support, so we
# need to try_import it.
_alchemlyb = _try_import("alchemlyb")

if _have_imported(_alchemlyb):
    from alchemlyb.parsing.gmx import extract as _extract

from .. import _gmx_exe
from .. import _isVerbose
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from ..Types._type import Type as _Type

from .. import IO as _IO
from .. import Protocol as _Protocol
from .. import Trajectory as _Trajectory
from .. import Types as _Types
from .. import Units as _Units
from .. import _Utils

from . import _process

from ._plumed import Plumed as _Plumed


class Gromacs(_process.Process):
    """A class for running simulations using GROMACS."""

    def __init__(
        self,
        system,
        protocol,
        exe=None,
        name="gromacs",
        work_dir=None,
        seed=None,
        extra_options=None,
        extra_lines=None,
        reference_system=None,
        property_map={},
        restraint=None,
        ignore_warnings=False,
        show_errors=True,
        checkpoint_file=None,
    ):
        """Constructor.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.

        protocol : :class:`Protocol <BioSimSpace.Protocol>`
            The protocol for the GROMACS process.

        exe : str
            The full path to the GROMACS executable.

        name : str
            The name of the process.

        work_dir :
            The working directory for the process.

        seed : int
            A random number seed.

        extra_options : dict
            A dictionary containing extra options. Overrides the ones generated from the protocol.

        extra_lines : list
            A list of extra lines to be put at the end of the script.

        reference_system : :class:`System <BioSimSpace._SireWrappers.System>` or None
            An optional system to use as a source of reference coordinates, if applicable.
            It is assumed that this system has the same topology as "system". If this is
            None, then "system" is used as a reference.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        restraint : :class:`Restraint <BioSimSpace.FreeEnergy.Restraint>`
            The Restraint object that contains information for the ABFE
            calculations.

        ignore_warnings : bool
            Whether to ignore warnings when generating the binary run file
            with 'gmx grompp'. By default, these warnings are elevated to
            errors and will halt the program.

        show_errors : bool
            Whether to show warning/error messages when generating the binary
            run file.

        checkpoint_file : str
           The path to a checkpoint file from a previous run. This can be used
           to continue an existing simulation. Currently we only support the
           use of checkpoint files for Equilibration protocols.
        """

        # Call the base class constructor.
        super().__init__(
            system,
            protocol,
            name,
            work_dir,
            seed,
            extra_options,
            extra_lines,
            property_map,
            restraint,
        )

        # Set the package name.
        self._package_name = "GROMACS"

        # This process can generate trajectory data.
        self._has_trajectory = True

        # Use GROMACS executable from environment.
        if exe is None:
            if _gmx_exe is not None:
                self._exe = _gmx_exe
            else:
                raise _MissingSoftwareError(
                    "'BioSimSpace.Process.Gromacs' is not supported. "
                    "Please install GROMACS (http://www.gromacs.org)."
                )
        # Use user-specified executable.
        else:
            # Make sure executable exists.
            if _os.path.isfile(exe):
                self._exe = exe
            else:
                raise IOError("GROMACS executable doesn't exist: '%s'" % exe)

        if not isinstance(ignore_warnings, bool):
            raise ValueError("'ignore_warnings' must be of type 'bool'.")
        self._ignore_warnings = ignore_warnings

        if not isinstance(show_errors, bool):
            raise ValueError("'show_errors' must be of type 'bool'.")
        self._show_errors = show_errors

        if restraint and not isinstance(protocol, _Protocol._FreeEnergyMixin):
            raise ValueError(
                "'BioSimSpace.Process.Gromacs' requires a "
                "FreeEnergy protocol for running with a restraint!"
            )

        # Initialise the energy dictionary and title header.
        self._energy_dict = (
            dict()
        )  # cannot figure out how to set value for _process._MultiDict()

        # Store the name of the GROMACS log file.
        self._log_file = "%s/%s.log" % (self._work_dir, name)
        self._eng_file = "%s/%s.edr" % (self._work_dir, name)

        # The names of the input files.
        self._gro_file = "%s/%s.gro" % (self._work_dir, name)
        self._top_file = "%s/%s.top" % (self._work_dir, name)

        # The name of the trajectory file.
        self._traj_file = "%s/%s.trr" % (self._work_dir, name)

        # The name of the output coordinate file.
        self._crd_file = "%s/%s_out.gro" % (self._work_dir, name)

        # Set the path for the GROMACS configuration file.
        self._config_file = "%s/%s.mdp" % (self._work_dir, name)

        # Set the reference system
        self._ref_file = f"{self._work_dir}/{name}_ref.gro"
        self._ref_system = reference_system

        # Create the list of input files.
        self._input_files = [self._config_file, self._gro_file, self._top_file]

        # Initialise the PLUMED interface object.
        self._plumed = None

        # Set the path of Gromacs checkpoint file.
        self._checkpoint_file = None
        if checkpoint_file is not None:
            if not isinstance(checkpoint_file, str):
                raise ValueError("'checkpoint_file' must be of type 'str'.")
            else:
                if _os.path.isfile(checkpoint_file):
                    self._checkpoint_file = checkpoint_file
                else:
                    raise IOError(
                        "GROMACS checkpoint file doesn't exist: '%s'" % checkpoint_file
                    )

        # Now set up the working directory for the process.
        self._setup()

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...
        self._write_system(
            self._system, coord_file=self._gro_file, topol_file=self._top_file
        )

        # Create the reference file
        if self._ref_system is not None:
            if (
                self._system.getAlchemicalIon()
                or self._protocol.getRestraint() is not None
            ):
                self._write_system(self._ref_system, ref_file=self._ref_file)
        else:
            _shutil.copy(self._gro_file, self._ref_file)

        # Create the binary input file name.
        self._tpr_file = "%s/%s.tpr" % (self._work_dir, self._name)
        self._input_files.append(self._tpr_file)

        # Generate the GROMACS configuration file.
        # Skip if the user has passed a custom config.
        if isinstance(self._protocol, _Protocol.Custom):
            self.setConfig(self._protocol.getConfig())
        else:
            self._generate_config()
        self.writeConfig(self._config_file)

        # Generate the dictionary of command-line arguments.
        self._generate_args()

        # Return the list of input files.
        return self._input_files

    def _write_system(self, system, coord_file=None, topol_file=None, ref_file=None):
        """Validates an input system and makes some internal modifications to it,
        if needed, before writing it out to a coordinate and/or a topology file.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.

        coord_file : str or None
            The coordinate file to which to write out the system.

        topol_file : str or None
            The topology file to which to write out the system.

        ref_file : str or None
            The file to which to write out the reference system for position restraints.
        """
        # Create a copy of the system.
        system = system.copy()

        if isinstance(self._protocol, _Protocol._FreeEnergyMixin):
            # Check that the system contains a perturbable molecule.
            if (
                system.nPerturbableMolecules() == 0
                and system.nDecoupledMolecules() == 0
            ):
                raise ValueError(
                    "'BioSimSpace.Protocol.FreeEnergy' requires a "
                    "perturbable molecule!"
                )

            # Check that the perturbation type is supported..
            if self._protocol.getPerturbationType() not in [
                "full",
                "release_restraint",
            ]:
                msg = (
                    "'BioSimSpace.Process.Gromacs' currently only supports the 'full' "
                    " and 'release_restraint' perturbation type. Please use 'BioSimSpace.Process.Somd' "
                    "for multistep perturbation types."
                )
                raise NotImplementedError(msg)

            # Check that we have multiple distance restraints if the perturbation type is 'release_restraint'.
            if self._protocol.getPerturbationType() == "release_restraint":
                if not self._restraint:
                    raise ValueError(
                        "'BioSimSpace.Process.Gromacs' requires a "
                        "restraint for the 'release_restraint' perturbation type!"
                    )
                if self._restraint._restraint_type != "multiple_distance":
                    raise ValueError(
                        "'BioSimSpace.Process.Gromacs' requires a "
                        "multiple distance restraint for the 'release_restraint' perturbation type!"
                    )

        else:
            # Check for perturbable molecules and convert to the chosen end state.
            system = self._checkPerturbable(system)

        # Convert the water model topology so that it matches the GROMACS naming convention.
        system._set_water_topology("GROMACS", property_map=self._property_map)

        # Check whether the system contains periodic box information.
        # For now, we'll not attempt to generate a box if the system property
        # is missing. If no box is present, we'll assume a non-periodic simulation.
        if "space" in system._sire_object.propertyKeys():
            has_box = True
        else:
            _warnings.warn("No simulation box found. Assuming gas phase simulation.")
            has_box = False

        # Deal with PBC.
        if not has_box or not self._has_water:
            # Create a 999.9 nm periodic box and apply to the system.
            space = _SireVol.PeriodicBox(_SireMaths.Vector(9999, 9999, 9999))
            system._sire_object.setProperty(
                self._property_map.get("space", "space"), space
            )

        # GRO87 coordinate files.
        if coord_file is not None:
            file = _os.path.splitext(coord_file)[0]
            _IO.saveMolecules(
                file,
                system,
                "gro87",
                match_waters=False,
                property_map=self._property_map,
            )

        # GRO87 reference files.
        if ref_file is not None:
            file = _os.path.splitext(ref_file)[0]
            _IO.saveMolecules(
                file,
                system,
                "gro87",
                match_waters=False,
                property_map=self._property_map,
            )

        # TOP file.
        if topol_file is not None:
            file = _os.path.splitext(topol_file)[0]
            _IO.saveMolecules(
                file,
                system,
                "grotop",
                match_waters=False,
                property_map=self._property_map,
            )
            self._apply_ABFE_restraint()

    def _apply_ABFE_restraint(self):
        # Write the restraint to the topology file
        if self._restraint:
            with open(self._top_file, "a") as f:
                f.write("\n")
                f.write(
                    self._restraint.toString(
                        engine="GROMACS",
                        perturbation_type=self._protocol.getPerturbationType(),
                        restraint_lambda="restraint"
                        in self._protocol.getLambda(type="series"),
                    )
                )

    def _generate_config(self):
        """Generate GROMACS configuration file strings."""

        # Clear the existing configuration list.
        self._config = []

        config_options = {}
        if not isinstance(self._protocol, _Protocol.Minimisation):
            # Set the random number seed.
            if self._is_seeded:
                seed = self._seed
            else:
                seed = -1
            config_options["ld-seed"] = seed

        if isinstance(self._protocol, _Protocol.Equilibration):
            if self._checkpoint_file is not None:
                config_options["continuation"] = "yes"

        if isinstance(self._protocol, _Protocol._PositionRestraintMixin):
            # Add any position restraints.
            self._add_position_restraints(config_options)

        # Add configuration variables for a metadynamics simulation.
        if isinstance(self._protocol, _Protocol.Metadynamics):
            # Create the PLUMED input file and copy auxiliary files to the working directory.
            self._plumed = _Plumed(str(self._work_dir))
            plumed_config, auxiliary_files = self._plumed.createConfig(
                self._system, self._protocol, self._property_map
            )
            self._setPlumedConfig(plumed_config)
            if auxiliary_files is not None:
                for file in auxiliary_files:
                    file_name = _os.path.basename(file)
                    _shutil.copyfile(file, self._work_dir + f"/{file_name}")
            self._input_files.append(self._plumed_config_file)

            # Expose the PLUMED specific member functions.
            setattr(self, "getPlumedConfig", self._getPlumedConfig)
            setattr(self, "getPlumedConfigFile", self._getPlumedConfigFile)
            setattr(self, "setPlumedConfig", self._setPlumedConfig)
            setattr(self, "getFreeEnergy", self._getFreeEnergy)
            setattr(self, "getCollectiveVariable", self._getCollectiveVariable)
            setattr(self, "sampleConfigurations", self._sampleConfigurations)
            setattr(self, "getTime", self._getTime)

        # Add configuration variables for a steered molecular dynamics protocol.
        elif isinstance(self._protocol, _Protocol.Steering):
            # Create the PLUMED input file and copy auxiliary files to the working directory.
            self._plumed = _Plumed(str(self._work_dir))
            plumed_config, auxiliary_files = self._plumed.createConfig(
                self._system, self._protocol, self._property_map
            )
            self._setPlumedConfig(plumed_config)
            if auxiliary_files is not None:
                for file in auxiliary_files:
                    file_name = _os.path.basename(file)
                    _shutil.copyfile(file, self._work_dir + f"/{file_name}")
            self._input_files.append(self._plumed_config_file)

            # Expose the PLUMED specific member functions.
            setattr(self, "getPlumedConfig", self._getPlumedConfig)
            setattr(self, "getPlumedConfigFile", self._getPlumedConfigFile)
            setattr(self, "setPlumedConfig", self._setPlumedConfig)
            setattr(self, "getCollectiveVariable", self._getCollectiveVariable)
            setattr(self, "getTime", self._getTime)

        # Set the configuration.
        if not isinstance(self._protocol, _Protocol.Dummy):
            config = _Protocol.ConfigFactory(self._system, self._protocol)
            pert_type = (
                self._protocol._perturbation_type
                if isinstance(self._protocol, _Protocol._FreeEnergyMixin)
                else None
            )
            self._apply_ABFE_restraint()
            self.addToConfig(
                config.generateGromacsConfig(
                    extra_options={**config_options, **self._extra_options},
                    extra_lines=self._extra_lines,
                    restraint=self._restraint,
                    perturbation_type=pert_type,
                )
            )

            # Flag that this isn't a custom protocol.
            self._protocol._setCustomised(False)

    def _generate_args(self):
        """Generate the dictionary of command-line arguments."""

        # Clear the existing arguments.
        self.clearArgs()

        # Add the default arguments.
        self.setArg("mdrun", True)  # Use mdrun.
        self.setArg("-deffnm", self._name)  # Output file prefix.
        self.setArg("-c", self._crd_file)  # Output out coordinate file.

        # Metadynamics and steered MD arguments.
        if isinstance(self._protocol, (_Protocol.Metadynamics, _Protocol.Steering)):
            self.setArg("-plumed", "plumed.dat")

    @staticmethod
    def _generate_binary_run_file(
        mdp_file,
        gro_file,
        top_file,
        ref_file,
        tpr_file,
        exe,
        checkpoint_file=None,
        ignore_warnings=False,
        show_errors=True,
    ):
        """
        Use grommp to generate the binary run input file.

        Parameters
        ----------

        mdp_file : str
            The path to the input mdp file.

        gro_file : str
            The path to the input coordinate file.

        top_file : str
            The path to the input topology file.

        ref_file : str
            The path to the input reference coordinate file to be used for
            position restraints.

        tpr_file : str
            The path to the output binary run file.

        exe : str
            The path to the GROMACS executable.

        checkpoint_file : str
            The path to a checkpoint file from a previous run. This can be used
            to continue an existing simulation. Currently we only support the
            use of checkpoint files for Equilibration protocols.

        ignore_warnings : bool
            Whether to ignore warnings when generating the binary run file
            with 'gmx grompp'. By default, these warnings are elevated to
            errors and will halt the program.

        show_errors : bool
            Whether to show warning/error messages when generating the binary
            run file.
        """

        if not isinstance(mdp_file, str):
            raise ValueError("'mdp_file' must be of type 'str'.")
        if not _os.path.isfile(mdp_file):
            raise IOError(f"'mdp_file' doesn't exist: '{mdp_file}'")

        if not isinstance(gro_file, str):
            raise ValueError("'gro_file' must be of type 'str'.")
        if not _os.path.isfile(gro_file):
            raise IOError(f"'gro_file' doesn't exist: '{gro_file}'")

        if not isinstance(top_file, str):
            raise ValueError("'top_file' must be of type 'str'.")
        if not _os.path.isfile(top_file):
            raise IOError(f"'top_file' doesn't exist: '{top_file}'")

        if not isinstance(ref_file, str):
            raise ValueError("'ref_file' must be of type 'str'.")
        if not _os.path.isfile(ref_file):
            raise IOError(f"'ref_file' doesn't exist: '{ref_file}'")

        if not isinstance(tpr_file, str):
            raise ValueError("'tpr_file' must be of type 'str'.")

        if not isinstance(exe, str):
            raise ValueError("'exe' must be of type 'str'.")
        if not _os.path.isfile(exe):
            raise IOError(f"'exe' doesn't exist: '{exe}'")

        if checkpoint_file is not None:
            if not isinstance(checkpoint_file, str):
                raise ValueError("'checkpoint_file' must be of type 'str'.")
            if not _os.path.isfile(checkpoint_file):
                raise IOError(f"'checkpoint_file' doesn't exist: '{checkpoint_file}'")

        if not isinstance(ignore_warnings, bool):
            raise ValueError("'ignore_warnings' must be of type 'bool'")

        if not isinstance(show_errors, bool):
            raise ValueError("'show_errors' must be of type 'bool'")

        # Create the name of the output mdp file.
        mdp_out = (
            _os.path.dirname(mdp_file)
            + "/%s.out.mdp" % _os.path.basename(mdp_file).split(".")[0]
        )

        # Use grompp to generate the portable binary run input file.
        if checkpoint_file is not None:
            command = "%s grompp -f %s -po %s -c %s -p %s -r %s -t %s -o %s" % (
                exe,
                mdp_file,
                mdp_out,
                gro_file,
                top_file,
                ref_file,
                checkpoint_file,
                tpr_file,
            )
        else:
            command = "%s grompp -f %s -po %s -c %s -p %s -r %s -o %s" % (
                exe,
                mdp_file,
                mdp_out,
                gro_file,
                top_file,
                ref_file,
                tpr_file,
            )

        # Warnings don't trigger an error.
        if ignore_warnings:
            command += " --maxwarn 9999"

        # Run the command.
        proc = _subprocess.run(
            _Utils.command_split(command),
            shell=False,
            text=True,
            stdout=_subprocess.PIPE,
            stderr=_subprocess.PIPE,
        )

        # Check that grompp ran successfully.
        if proc.returncode != 0:
            # Handle errors and warnings.
            if show_errors:
                # Capture errors and warnings from the grompp output.
                errors = []
                warnings = []
                is_error = False
                is_warn = False
                lines = proc.stderr.split("\n")
                for line in lines:
                    line = line.strip()
                    if line[0:5] == "ERROR" or is_error:
                        if line == "":
                            is_error = False
                            continue
                        errors.append(line)
                        is_error = True
                    elif line[0:7] == "WARNING" or is_warn:
                        if line == "":
                            is_warn = False
                            continue
                        warnings.append(line)
                        is_warn = True

                error_string = "\n  ".join(errors)
                warning_string = "\n  ".join(warnings)

                exception_string = "Unable to generate GROMACS binary run input file.\n"
                if len(errors) > 0:
                    exception_string += (
                        "\n'gmx grompp' reported the following errors:\n"
                        + f"{error_string}\n"
                    )
                if len(warnings) > 0:
                    exception_string += (
                        "\n'gmx grompp' reported the following warnings:\n"
                        + f"{warning_string}\n"
                        + "\nUse 'ignore_warnings' to ignore warnings."
                    )

                raise RuntimeError(exception_string)

            else:
                raise RuntimeError(
                    "Unable to generate GROMACS binary run input file. "
                    "Use 'show_errors=True' to display errors/warnings."
                )

    def addToConfig(self, config):
        """
        Add a string to the configuration list.

        Parameters
        ----------

        config : str, [str]
            A configuration string, a list of configuration strings, or a
            path to a configuration file.
        """

        # Call the base class method.
        super().addToConfig(config)

        # Use grompp to generate the portable binary run input file.
        self._generate_binary_run_file(
            self._config_file,
            self._gro_file,
            self._top_file,
            self._ref_file,
            self._tpr_file,
            self._exe,
            checkpoint_file=self._checkpoint_file,
            ignore_warnings=self._ignore_warnings,
            show_errors=self._show_errors,
        )

    def resetConfig(self):
        """Reset the configuration parameters."""
        self._generate_config()

        # Use grompp to generate the portable binary run input file.
        self._generate_binary_run_file(
            self._config_file,
            self._gro_file,
            self._top_file,
            self._ref_file,
            self._tpr_file,
            self._exe,
            checkpoint_file=self._checkpoint_file,
            ignore_warnings=self._ignore_warnings,
            show_errors=self._show_errors,
        )

    def setConfig(self, config):
        """
        Set the list of configuration file strings.

        Parameters
        ----------

        config : str, [str]
            The list of configuration strings, or a path to a configuration
            file.
        """

        # Call the base class method.
        super().setConfig(config)

        # Use grompp to generate the portable binary run input file.
        self._generate_binary_run_file(
            self._config_file,
            self._gro_file,
            self._top_file,
            self._ref_file,
            self._tpr_file,
            self._exe,
            checkpoint_file=self._checkpoint_file,
            ignore_warnings=self._ignore_warnings,
            show_errors=self._show_errors,
        )

    def start(self):
        """
        Start the GROMACS process.

        Returns
        -------

        process : :class:`Process.Gromacs <BioSimSpace.Process.Gromacs>`
            A handle to the GROMACS process.
        """

        # The process is currently queued.
        if self.isQueued():
            return

        # Process is already running.
        if self._process is not None:
            if self._process.isRunning():
                return

        # Clear any existing output.
        self._clear_output()

        # Run the process in the working directory.
        with _Utils.cd(self._work_dir):
            # Create the arguments string list.
            args = self.getArgStringList()

            # Write the command-line process to a README.txt file.
            with open("README.txt", "w") as f:
                # Set the command-line string.
                self._command = "%s " % self._exe + self.getArgString()

                # Write the command to file.
                f.write("# GROMACS was run with the following command:\n")
                f.write("%s\n" % self._command)

            # Start the timer.
            self._timer = _timeit.default_timer()

            # Start the simulation.
            self._process = _SireBase.Process.run(
                self._exe, args, "%s.out" % self._name, "%s.out" % self._name
            )

            # For historical reasons (console message aggregation with MPI), Gromacs
            # writes the majority of its output to stderr. For user convenience, we
            # redirect all output to stdout, and place a message in the stderr file
            # to highlight this.
            with open(self._stderr_file, "w") as f:
                f.write("All output has been redirected to the stdout stream!\n")

        return self

    def getSystem(self, block="AUTO"):
        """
        Get the latest molecular system.

        Parameters
        ----------

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The latest molecular system.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()
            block = True

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        if block is True and not self.isError():
            return self._getFinalFrame()
        else:
            # Minimisation trajectories have a single frame, i.e. the final state.
            if isinstance(self._protocol, _Protocol.Minimisation):
                time = 0 * _Units.Time.nanosecond
            # Get the current simulation time.
            else:
                time = self.getTime()

            # Grab the most recent frame from the trajectory file.
            return self._getFrame(time)

    def getCurrentSystem(self):
        """
        Get the latest molecular system.

        Returns
        -------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The latest molecular system.
        """
        return self.getSystem(block=False)

    def getTrajectory(self, backend="AUTO", block="AUTO"):
        """
        Return a trajectory object.

        Parameters
        ----------

        backend : str
            The backend to use for trajectory parsing. To see supported backends,
            run BioSimSpace.Trajectory.backends(). Using "AUTO" will try each in
            sequence.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        trajectory : :class:`System <BioSimSpace.Trajectory.Trajectory>`
            The latest trajectory object.
        """

        if not isinstance(backend, str):
            raise TypeError("'backend' must be of type 'str'")

        if not isinstance(block, (bool, str)):
            raise TypeError("'block' must be of type 'bool' or 'str'")

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        try:
            # Locate the trajectory file.
            traj_file = self._find_trajectory_file()

            if traj_file is None:
                return None
            else:
                self._traj_file = traj_file

            return _Trajectory.Trajectory(process=self, backend=backend)

        except:
            return None

    def getFrame(self, index):
        """
        Return a specific trajectory frame.

        Parameters
        ----------

        index : int
            The index of the frame.

        Returns
        -------

        frame : :class:`System <BioSimSpace._SireWrappers.System>`
            The System object of the corresponding frame.
        """

        if not type(index) is int:
            raise TypeError("'index' must be of type 'int'")

        max_index = (
            int(
                (self._protocol.getRunTime() / self._protocol.getTimeStep())
                / self._protocol.getRestartInterval()
            )
            - 1
        )

        if index < 0 or index > max_index:
            raise ValueError(f"'index' must be in range [0, {max_index}].")

        try:
            time = (
                index
                * self._protocol.getRestartInterval()
                * self._protocol.getTimeStep()
            )

            with _warnings.catch_warnings():
                system = self._getFrame(time)

            return self._getFrame(time)

        except:
            return None

    def getRecord(self, record, time_series=False, unit=None, block="AUTO"):
        """
        Get a record from the stdout dictionary.

        Parameters
        ----------

        record : str
            The record key.

        time_series : bool
            Whether to return a list of time series records.

        unit : :class:`Unit <BioSimSpace.Units>`
            The unit to convert the record to.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        record : :class:`Type <BioSimSpace.Types>`
            The matching record.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        self._update_energy_dict()
        return self._get_energy_record(record, time_series, unit)

    def getCurrentRecord(self, record, time_series=False, unit=None):
        """
        Get a current record from the stdout dictionary.

        Parameters
        ----------

        record : str
            The record key.

        time_series : bool
            Whether to return a list of time series records.

        unit : :class:`Unit <BioSimSpace.Units>`
            The unit to convert the record to.

        Returns
        -------

        record : :class:`Type <BioSimSpace.Types>`
            The matching record.
        """
        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        self._update_energy_dict()
        return self._get_energy_record(record, time_series, unit)

    def getRecords(self, block="AUTO"):
        """
        Return the dictionary of stdout time-series records.

        Parameters
        ----------

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        records : :class:`MultiDict <BioSimSpace.Process._process._MultiDict>`
           The dictionary of time-series records.
        """
        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        self._update_energy_dict()
        return self._energy_dict.copy()

    def getCurrentRecords(self):
        """
        Return the current dictionary of stdout time-series records.

        Parameters
        ----------

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        records : :class:`MultiDict <BioSimSpace.Process._process._MultiDict>`
           The dictionary of time-series records.
        """
        return self.getRecords(block=False)

    def getTime(self, time_series=False, block="AUTO"):
        """
        Get the simulation time.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The current simulation time in nanoseconds.
        """

        return self.getRecord("TIME", time_series, _Units.Time.picosecond, block)

    def getCurrentTime(self, time_series=False):
        """
        Get the current simulation time.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The current simulation time in nanoseconds.
        """
        return self.getTime(time_series, block=False)

    def getStep(self, time_series=False, block="AUTO"):
        """
        Get the number of integration steps.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        step : int
            The current number of integration steps.

        Notes
        -----
        The step is calculated based on
        :meth:`~BioSimSpace.Process.Gromacs.getTime` and
        :meth:`~BioSimSpace.Protocol.getTimeStep`.
        """
        records = self.getRecord("TIME", time_series, _Units.Time.picosecond, block)
        time_step = self._protocol.getTimeStep()
        if isinstance(records, list):
            return [record / time_step for record in records]
        else:
            return records / time_step

    def getCurrentStep(self, time_series=False):
        """
        Get the current number of integration steps.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        step : int
            The current number of integration steps.
        """
        return self.getStep(time_series, block=False)

    def getBondEnergy(self, time_series=False, block="AUTO"):
        """
        Get the bond energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The bond energy.
        """
        return self.getRecord("BOND", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentBondEnergy(self, time_series=False):
        """
        Get the current bond energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The bond energy.
        """
        return self.getBondEnergy(time_series, block=False)

    def getAngleEnergy(self, time_series=False, block="AUTO"):
        """
        Get the angle energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The angle energy.
        """
        return self.getRecord("ANGLE", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentAngleEnergy(self, time_series=False):
        """
        Get the current angle energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The angle energy.
        """
        return self.getAngleEnergy(time_series, block=False)

    def getDihedralEnergy(self, time_series=False, block="AUTO"):
        """
        Get the total dihedral energy (proper + improper).

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The total dihedral energy.
        """
        # Get the proper and improper energies.
        proper = self.getRecord(
            "PROPERDIH", time_series, _Units.Energy.kj_per_mol, block
        )
        improper = self.getRecord(
            "IMPROPERDIH", time_series, _Units.Energy.kj_per_mol, block
        )

        # No records.
        if proper is None and improper is None:
            return None
        elif proper is None:
            return improper
        elif improper is None:
            return proper
        else:
            if time_series:
                return [x + y for x, y in zip(proper, improper)]
            else:
                return proper + improper

    def getCurrentDihedralEnergy(self, time_series=False):
        """
        Get the current total dihedral energy (proper + improper).

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The dihedral energy.
        """
        return self.getDihedralEnergy(time_series, block=False)

    def getProperEnergy(self, time_series=False, block="AUTO"):
        """
        Get the proper dihedral energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The proper dihedral energy.
        """
        return self.getRecord("PROPERDIH", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentProperEnergy(self, time_series=False):
        """
        Get the current proper dihedral energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The proper dihedral energy.
        """
        return self.getProperEnergy(time_series, block=False)

    def getImproperEnergy(self, time_series=False, block="AUTO"):
        """
        Get the improper energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The improper energy.
        """
        return self.getRecord(
            "IMPROPERDIH", time_series, _Units.Energy.kj_per_mol, block
        )

    def getCurrentImproperEnergy(self, time_series=False):
        """
        Get the current improper energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The improper energy.
        """
        return self.getImproperEnergy(time_series, block=False)

    def getLennardJones14(self, time_series=False, block="AUTO"):
        """
        Get the Lennard-Jones energy between atoms 1 and 4.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The Lennard-Jones energy.
        """
        return self.getRecord("LJ14", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentLennardJones14(self, time_series=False):
        """
        Get the current Lennard-Jones energy between atoms 1 and 4.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The Lennard-Jones energy.
        """
        return self.getLennardJones14(time_series, block=False)

    def getLennardJonesSR(self, time_series=False, block="AUTO"):
        """
        Get the short-range Lennard-Jones energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The short-range Lennard-Jones energy.
        """
        return self.getRecord("LJSR", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentLennardJonesSR(self, time_series=False):
        """
        Get the current short-range Lennard-Jones energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The Lennard-Jones energy.
        """
        return self.getLennardJonesSR(time_series, block=False)

    def getCoulomb14(self, time_series=False, block="AUTO"):
        """
        Get the Coulomb energy between atoms 1 and 4.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The Coulomb energy.
        """
        return self.getRecord("COULOMB14", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentCoulomb14(self, time_series=False):
        """
        Get the current Coulomb energy between atoms 1 and 4.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The Coulomb energy.
        """
        return self.getCoulomb14(time_series, block=False)

    def getCoulombSR(self, time_series=False, block="AUTO"):
        """
        Get the short-range Coulomb energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The Coulomb energy.
        """
        return self.getRecord("COULOMBSR", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentCoulombSR(self, time_series=False):
        """
        Get the current short-range Coulomb energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The Coulomb energy.
        """
        return self.getCoulombSR(time_series, block=False)

    def getCoulombReciprocal(self, time_series=False, block="AUTO"):
        """
        Get the reciprocal space Coulomb energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The Coulomb energy.
        """
        return self.getRecord("COULRECIP", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentCoulombReciprocal(self, time_series=False):
        """
        Get the current reciprocal space Coulomb energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The Coulomb energy.
        """
        return self.getCoulombReciprocal(time_series, block=False)

    def getDispersionCorrection(self, time_series=False, block="AUTO"):
        """
        Get the dispersion correction.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The dispersion correction.
        """
        return self.getRecord(
            "DISPERCORR", time_series, _Units.Energy.kj_per_mol, block
        )

    def getCurrentDispersionCorrection(self, time_series=False):
        """
        Get the current dispersion correction.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The dispersion correction.
        """
        return self.getDispersionCorrection(time_series, block=False)

    def getRestraintEnergy(self, time_series=False, block="AUTO"):
        """
        Get the position restraint energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The dispersion correction.
        """
        return self.getRecord(
            "POSITIONREST", time_series, _Units.Energy.kj_per_mol, block
        )

    def getCurrentRestraintEnergy(self, time_series=False):
        """
        Get the current position restraint energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The dispersion correction.
        """
        return self.getRestraintEnergy(time_series, block=False)

    def getPotentialEnergy(self, time_series=False, block="AUTO"):
        """
        Get the potential energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The potential energy.
        """
        return self.getRecord("POTENTIAL", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentPotentialEnergy(self, time_series=False):
        """
        Get the current potential energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The potential energy.
        """
        return self.getPotentialEnergy(time_series, block=False)

    def getKineticEnergy(self, time_series=False, block="AUTO"):
        """
        Get the kinetic energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The kinetic energy.
        """
        return self.getRecord("KINETICEN", time_series, _Units.Energy.kj_per_mol, block)

    def getCurrentKineticEnergy(self, time_series=False):
        """
        Get the current kinetic energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The kinetic energy.
        """
        return self.getKineticEnergy(time_series, block=False)

    def getTotalEnergy(self, time_series=False, block="AUTO"):
        """
        Get the total energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The total energy.
        """
        return self.getRecord(
            "TOTALENERGY", time_series, _Units.Energy.kj_per_mol, block
        )

    def getCurrentTotalEnergy(self, time_series=False):
        """
        Get the current total energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The total energy.
        """
        return self.getTotalEnergy(time_series, block=False)

    def getConservedEnergy(self, time_series=False, block="AUTO"):
        """
        Get the conserved energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The conserved energy.
        """
        return self.getRecord(
            "CONSERVEDEN", time_series, _Units.Energy.kj_per_mol, block
        )

    def getCurrentConservedEnergy(self, time_series=False):
        """
        Get the current conserved energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The conserved energy.
        """
        return self.getConservedEnergy(time_series, block=False)

    def getTemperature(self, time_series=False, block="AUTO"):
        """
        Get the temperature.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The temperature.
        """
        return self.getRecord(
            "TEMPERATURE", time_series, _Units.Temperature.kelvin, block
        )

    def getCurrentTemperature(self, time_series=False):
        """
        Get the current temperature.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The current temperature.
        """
        return self.getTemperature(time_series, block=False)

    def getPressure(self, time_series=False, block="AUTO"):
        """
        Get the pressure.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure.
        """
        return self.getRecord("PRESSURE", time_series, _Units.Pressure.bar, block)

    def getCurrentPressure(self, time_series=False):
        """
        Get the current pressure.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The current pressure.
        """
        return self.getPressure(time_series, block=False)

    def getPressureDC(self, time_series=False, block="AUTO"):
        """
        Get the DC pressure.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The DC pressure.
        """
        return self.getRecord("PRESDC", time_series, _Units.Pressure.bar, block)

    def getCurrentPressureDC(self, time_series=False):
        """
        Get the current DC pressure.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The current pressure.
        """
        return self.getPressureDC(time_series, block=False)

    def getConstraintRMSD(self, time_series=False, block="AUTO"):
        """
        Get the RMSD of the constrained atoms.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        length : :class:`Length <BioSimSpace.Types.Length>`
            The constrained RMSD.
        """
        return self.getRecord("CONSTRRMSD", time_series, None, block)

    def getCurrentConstraintRMSD(self, time_series=False):
        """
        Get the current RMSD of the constrained atoms.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        length : :class:`Length <BioSimSpace.Types.Length>`
            The current constrained RMSD.
        """
        return self.getConstraintRMSD(time_series, block=False)

    def getVolume(self, time_series=False, block="AUTO"):
        """Get the volume.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        volume : :class:`Volume <BioSimSpace.Types.Volume>`
           The volume.
        """
        return self.getRecord("VOLUME", time_series, _Units.Volume.angstrom3, block)

    def getCurrentVolume(self, time_series=False):
        """Get the current volume.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        volume : :class:`Volume <BioSimSpace.Types.Volume>`
           The volume.
        """
        return self.getVolume(time_series, block=False)

    def stdout(self, n=10):
        """
        Print the last n lines of the stdout buffer.

        Parameters
        ----------

        n : int
            The number of lines to print.
        """

        # Note that thermodynamic records, e.g. energy, pressure, temperature,
        # are redirected to a log file.

        # Ensure that the number of lines is positive.
        if n < 0:
            raise ValueError("The number of lines must be positive!")

        # Append any new lines to the stdout list.
        for line in _pygtail.Pygtail(self._stdout_file):
            self._stdout.append(line.rstrip())

        # Get the current number of lines.
        num_lines = len(self._stdout)

        # Set the line from which to start printing.
        if num_lines < n:
            start = 0
        else:
            start = num_lines - n

        # Print the lines.
        for x in range(start, num_lines):
            print(self._stdout[x])

    def _add_position_restraints(self, config_options):
        """
        Helper function to add position restraints.

        Parameters
        ----------

        config_options : dict
            The dictionary of configuration options.
        """

        # Get the restraint type.
        restraint = self._protocol.getRestraint()

        if restraint is not None or self._system.getAlchemicalIon():
            # Get the force constant in units of kJ_per_mol/nanometer**2
            force_constant = self._protocol.getForceConstant()._sire_unit
            force_constant = force_constant.to(
                _SireUnits.kJ_per_mol / _SireUnits.nanometer2
            )

            # Scale reference coordinates with the scaling matrix of the pressure coupling.
            config_options["refcoord-scaling"] = "com"

            # Copy the user property map.
            property_map = self._property_map.copy()

            # Parse the topology in serial to ensure that molecules are
            # ordered correctly. Don't sort based on name.
            property_map["parallel"] = _SireBase.wrap(False)
            property_map["sort"] = _SireBase.wrap(False)

            # Create a copy of the system.
            system = self._system.copy()

            # Convert to the lambda = 0 state if this is a perturbable system and this
            # isn't a free energy protocol.
            if not isinstance(self._protocol, _Protocol._FreeEnergyMixin):
                system = self._checkPerturbable(system)

            # Convert the water model topology so that it matches the GROMACS naming convention.
            system._set_water_topology("GROMACS")

            # Create a GROMACS topology object.
            top = _SireIO.GroTop(system._sire_object, property_map)

            # Get the top file as a list of lines.
            top_lines = top.lines()

            # List of 'moleculetype' record indices.
            moltypes_top_idx = []

            # Store the line index for the start of each 'moleculetype' record.
            for idx, line in enumerate(top_lines):
                if "[ moleculetype ]" in line or "[ system ]" in line:
                    moltypes_top_idx.append(idx)

            # Create a dictionary to store the indices of the molecules in the
            # system for each GROMACS molecule type and the inverse.
            moltypes_sys_idx = {}
            sys_idx_moltypes = {}

            # Convert the topology to a GROMACS system.
            gro_system = top.groSystem()

            # Initialise the dictionary for each type.
            for mol_type in gro_system.uniqueTypes():
                moltypes_sys_idx[mol_type] = []

            # Now loop over each molecule and store the indices of the molecules
            # in the system that are of each type as well as a mapping from the
            # molecule index to the GROMACS molecule type.
            for idx, mol_type in enumerate(gro_system):
                moltypes_sys_idx[mol_type].append(idx)
                sys_idx_moltypes[idx] = mol_type

            if self._system.getAlchemicalIon():
                biggest_mol_idx = max(
                    range(system.nMolecules()), key=lambda x: system[x].nAtoms()
                )
            else:
                biggest_mol_idx = -1

            # A keyword restraint.
            if isinstance(restraint, str) or self._system.getAlchemicalIon():
                # The number of restraint files.
                num_restraint = 1

                # Loop over all of the molecule types and create a position
                # restraint file for each.
                for mol_type_idx, (mol_type, mol_idxs) in enumerate(
                    moltypes_sys_idx.items()
                ):
                    # Initialise a list of restrained atom indices.
                    restrained_atoms = []

                    # Loop over each molecule in the system that matches this
                    # type and append any atoms matching the restraint.
                    for idx, mol_idx in enumerate(mol_idxs):
                        # Get the indices of any restrained atoms in this molecule,
                        # making sure that indices are relative to the molecule.
                        if restraint is not None:
                            atom_idxs = self._system.getRestraintAtoms(
                                restraint,
                                mol_index=mol_idx,
                                is_absolute=False,
                                allow_zero_matches=True,
                            )
                        else:
                            atom_idxs = []

                        if self._system.getMolecule(mol_idx).isAlchemicalIon():
                            alch_ion = self._system.getMolecule(mol_idx).getAtoms()
                            alch_idx = alch_ion[0].index()
                            if alch_idx != 0 or len(alch_ion) != 1:
                                # The alchemical ions should only contain 1 atom
                                # and the relative index should thus be 0.
                                raise ValueError(
                                    f"{self._system.getMolecule(mol_idx)} is marked as an alchemical ion but has more than 1 atom."
                                )
                            else:
                                atom_idxs.append(alch_idx)

                        if mol_idx == biggest_mol_idx:
                            # Only triggered when there is alchemical ion present.
                            # The biggest_mol_idx is -1 when there is no alchemical ion.
                            protein_com_idx = self._system.getMolecule(
                                mol_idx
                            ).getCOMIdx()
                            if protein_com_idx not in atom_idxs:
                                atom_idxs.append(protein_com_idx)

                        # Store the atom index if it hasn't already been recorded.
                        for atom_idx in atom_idxs:
                            if not atom_idx in restrained_atoms:
                                restrained_atoms.append(atom_idx)

                    # Write the position restraint file for this molecule.
                    if len(restrained_atoms) > 0:
                        # Create the file names.
                        include_file = "posre_%04d.itp" % num_restraint
                        restraint_file = "%s/%s" % (self._work_dir, include_file)

                        with open(restraint_file, "w") as file:
                            # Write the header.
                            file.write("[ position_restraints ]\n")
                            file.write(";  i funct       fcx        fcy        fcz\n")

                            # Write restraints for each atom.
                            for atom_idx in restrained_atoms:
                                file.write(
                                    f"{atom_idx+1:4}    1       {force_constant}       {force_constant}       {force_constant}\n"
                                )

                        # Work out the offset.
                        offset = num_restraint - 1

                        # Include the position restraint file in the correct place within
                        # the topology file. We put the additional include directive at the
                        # end of the block so we move to the line before the next moleculetype
                        # record.
                        new_top_lines = top_lines[
                            : moltypes_top_idx[mol_type_idx + 1] + offset - 1
                        ]

                        # Append the additional information.
                        new_top_lines.append('#include "%s"' % include_file)
                        new_top_lines.append("")

                        # Now extend with the remainder of the file.
                        new_top_lines.extend(
                            top_lines[moltypes_top_idx[mol_type_idx + 1] + offset :]
                        )

                        # Overwrite the topology file lines.
                        top_lines = new_top_lines

                        # Increment the number of restraint files.
                        num_restraint += 1

                        # Append the restraint file to the list of autogenerated inputs.
                        self._input_files.append(restraint_file)

                # Write the updated topology to file.
                with open(self._top_file, "w") as file:
                    for line in top_lines:
                        file.write("%s\n" % line)

            # A user-defined list of atoms indices.
            else:
                # Create an empty multi-dict for each molecule type.
                mol_atoms = {}
                for mol_type in gro_system.uniqueTypes():
                    mol_atoms[mol_type] = []

                # Now work out which MolNum corresponds to each atom in the restraint.
                for idx in restraint:
                    try:
                        # Get the molecule index and relative atom index.
                        mol_idx, atom_idx = self._system._getRelativeIndices(idx)

                        # Get the type associated with this molecule index.
                        mol_type = sys_idx_moltypes[mol_idx]

                        # Append this atom if it's not already been recorded.
                        if not atom_idx in mol_atoms[mol_type]:
                            mol_atoms[mol_type].append(atom_idx)

                    except Exception as e:
                        msg = "Unable to find restrained atom in the system?"
                        if _isVerbose():
                            raise ValueError(msg) from e
                        else:
                            raise ValueError(msg) from None

                # The number of restraint files.
                num_restraint = 1

                # Loop over all of the molecule types and create a position
                # restraint file for each.
                for mol_type_idx, (mol_type, atom_idxs) in enumerate(mol_atoms.items()):
                    # Write the position restraint file for this molecule.
                    if len(atom_idxs) > 0:
                        # Create the file names.
                        include_file = "posre_%04d.itp" % num_restraint
                        restraint_file = "%s/%s" % (self._work_dir, include_file)

                        with open(restraint_file, "w") as file:
                            # Write the header.
                            file.write("[ position_restraints ]\n")
                            file.write(";  i funct       fcx        fcy        fcz\n")

                            # Write restraints for each atom.
                            for atom_idx in atom_idxs:
                                file.write(
                                    f"{atom_idx+1:4}    1       {force_constant}       {force_constant}       {force_constant}\n"
                                )

                        # Work out the offset.
                        offset = num_restraint - 1

                        # Include the position restraint file in the correct place within
                        # the topology file. We put the additional include directive at the
                        # end of the block so we move to the line before the next moleculetype
                        # record.
                        new_top_lines = top_lines[
                            : moltypes_top_idx[mol_type_idx + 1] + offset - 1
                        ]

                        # Append the additional information.
                        new_top_lines.append('#include "%s"' % include_file)
                        new_top_lines.append("")

                        # Now extend with the remainder of the file.
                        new_top_lines.extend(
                            top_lines[moltypes_top_idx[mol_type_idx + 1] + offset :]
                        )

                        # Overwrite the topology file lines.
                        top_lines = new_top_lines

                        # Increment the number of restraint files.
                        num_restraint += 1

                        # Append the restraint file to the list of autogenerated inputs.
                        self._input_files.append(restraint_file)

                # Write the updated topology to file.
                with open(self._top_file, "w") as file:
                    for line in top_lines:
                        file.write("%s\n" % line)

    def _initialise_energy_dict(self):
        # Grab the available energy terms
        command = f"{self._exe} energy -f {self._eng_file}"
        proc = _subprocess.run(
            _Utils.command_split(command),
            input="0",
            stdout=_subprocess.PIPE,
            stderr=_subprocess.PIPE,
            encoding="utf-8",
        )
        err = proc.stderr
        keys = self._parse_energy_terms(err)
        # We need to stored the original key as the one in the
        # self._energy_dict will be the sanitised keys.
        self._energy_keys = keys
        self._energy_dict["TIME"] = []
        for key in keys:
            # Skip surface tension records, since there is no appropriate general unit.
            if key != "#Surf*SurfTen":
                self._energy_dict[self._sanitise_energy_term(key)] = []

    @staticmethod
    def _parse_energy_terms(text):
        """Parse the output from gmx energy output to get the energy terms in
        the edr file. Example output look like:

        #                 :-) GROMACS - gmx energy, 2022.2-conda_forge (-:
        # Command line:
        #   gmx energy -f energy.edr
        # Opened prod.edr as single precision energy file
        # Select the terms you want from the following list by
        # selecting either (part of) the name or the number or a combination.
        # End your selection with an empty line or a zero.
        # -------------------------------------------------------------------
        #   1  Harmonic-Pot.    2  Angle            3  U-B              4  Proper-Dih.
        # -------------------------------------------------------

        Parameters
        ----------

        text : str
            The output string from the gmx energy

        Returns
        -------

        list
            A list of the string energy terms in the edr file.

        Notes
        -----
        The order that the key is stored is very important as the order that
        the energy term is stored in the xvg file will obey this order. In this
        case, the energy will be stored in the order of Harmonic-Pot., Angle,
        U-B, Proper-Dih.. Note that this order is absolute and will not be
        changed by the input to `gmx energy`.
        """
        sections = text.split("---")
        # Remove the empty sections
        sections = [section for section in sections if section]
        # Concatenate the lines
        section = sections[1].replace("\n", "")
        terms = section.split()
        # Remove the possible '-' from the separation line
        terms = [term for term in terms if term != "-"]
        # Check if the index order is correct
        indexes = [int(term) for term in terms[::2]]
        energy_names = terms[1::2]
        length_nomatch = len(indexes) != len(energy_names)
        # -1 as the index is 1-based.
        index_nomatch = (_np.arange(len(indexes)) != _np.array(indexes) - 1).any()
        if length_nomatch or index_nomatch:
            raise ValueError(f"Cannot parse the energy terms in the {edr_file} file.")
        else:
            return energy_names

    @staticmethod
    def _parse_energy_units(text):
        """Extract the energy unit from the output. Example outputs are:

        # Statistics over 15000001 steps [ 0.0000 through 30000.0000 ps ], 53 data sets
        # All statistics are over 15001 points (frames)
        # Energy                      Average   Err.Est.       RMSD  Tot-Drift
        # -------------------------------------------------------------------------------
        # Harmonic Pot.               0.34534      0.009   0.491375 -0.0225767  (kJ/mol)
        # Angle                       16159.4        4.7    202.124   -35.3327  (kJ/mol)
        # U-B                         81396.5        3.7    431.584    9.91613  (kJ/mol)
        # Proper Dih.                 71615.4         24    250.353    -101.74  (kJ/mol)

        Parameters
        ----------
        text : str
            Output text with term name and units.

        Returns
        -------
        list
            A list of the energy units of type :mod:`~BioSimSpace.Types._GeneralUnit`.

        Notes
        -----
        The order that the energy unit is printed will obey the order obtained
        from :meth:`~BioSimSpace.Process.Gromacs._parse_energy_terms`.
        """
        section = text.split("---")[-1]
        lines = section.split("\n")
        units = [
            _Units.Time.picosecond,
        ]
        for line in lines:
            terms = line.split()
            if len(terms) > 1 and terms[0] != "#Surf*SurfTen":
                unit = terms[-1][1:-1]
                if unit == "K":
                    units.append(_Units.Temperature.kelvin)
                elif unit == "kJ/mol":
                    units.append(_Units.Energy.kj_per_mol)
                elif unit == "bar":
                    units.append(_Units.Pressure.bar)
                elif unit == "":
                    units.append(_Units.Length.nanometer)
                elif unit == "nm":
                    units.append(_Units.Length.nanometer)
                elif unit == "nm^3":
                    units.append(_Units.Volume.nanometer3)
                elif unit == "bar nm":
                    units.append(_Units.Pressure.bar * _Units.Length.nanometer)
                elif unit == "nm/ps":
                    units.append(_Units.Length.nanometer / _Units.Time.picosecond)
                elif unit == "kg/m^3":
                    units.append(_Types._GeneralUnit("kg/m3"))
                else:
                    units.append(1.0)
                    _warnings.warn(
                        f"Unit {unit} cannot be parsed, recording the unit as unitless."
                    )
        return units

    @staticmethod
    def _sanitise_energy_term(key):
        """Format the energy term names to compile with the BioSimSpace
        standard.

        Parameters
        ----------
        key : str
            The original name of the energy term.

        Returns
        -------
        str
            The formatted name of the energy term.
        """
        # Convert to upper case.
        key = key.upper()

        # Strip whitespace and newlines from beginning and end.
        key = key.strip()

        # Remove whitespace.
        key = key.replace(" ", "")

        # Remove periods.
        key = key.replace(".", "")

        # Remove hyphens.
        key = key.replace("-", "")

        # Remove parentheses.
        key = key.replace("(", "")
        key = key.replace(")", "")

        # Remove instances of BAR.
        key = key.replace("BAR", "")
        return key

    def _update_energy_dict(self, initialise=False):
        if initialise or len(self._energy_dict) == 0:
            self._initialise_energy_dict()

        keys = self._energy_keys

        with _TemporaryDirectory() as tmpdirname:
            temp_dir = _Path(tmpdirname)
            output_file = temp_dir / "energy.xvg"
            command = f"{self._exe} energy -f {self._eng_file} -o {output_file}"
            proc = _subprocess.run(
                _Utils.command_split(command),
                input="\n".join(keys),
                # The order that the input keys are generated is irrelavent.
                # The order the energy term will be printed obeys
                # :meth:`~BioSimSpace.Process.Gromacs._parse_energy_terms`
                stdout=_subprocess.PIPE,
                stderr=_subprocess.PIPE,
                encoding="utf-8",
            )
            out = proc.stdout
            results = _np.loadtxt(output_file, comments=["@", "#"])
            units = self._parse_energy_units(out)

            if len(units) != len(list(self._energy_dict)):
                raise ValueError(
                    "The number of energy units does not match the "
                    "number of energy terms."
                )

        for i, key in enumerate(self._energy_dict):
            if len(results.shape) == 1:
                # Account for the case of single point energy
                self._energy_dict[key] = [results[i] * units[i]]
            else:
                self._energy_dict[key] = [result * units[i] for result in results[:, i]]

    def _get_energy_record(self, key, time_series=False, unit=None):
        """Helper function to get a stdout record from the dictionary.

        Parameters
        ----------

        key : str
            The record key.

        time_series : bool
            Whether to return a time series of records.

        unit : BioSimSpace.Types._type.Type
            The unit to convert the record to.

        Returns
        -------

        record :
            The matching stdout record.
        """

        # No data!
        if len(self._energy_dict) == 0:
            return None

        if not isinstance(time_series, bool):
            _warnings.warn("Non-boolean time-series flag. Defaulting to False!")
            time_series = False

        # Validate the unit.
        if unit is not None:
            if not isinstance(unit, _Type):
                raise TypeError("'unit' must be of type 'BioSimSpace.Types'")

        # Return the list of dictionary values.
        if time_series:
            try:
                if unit is None:
                    return [
                        x._to_default_unit().value() for x in self._energy_dict[key]
                    ]
                else:
                    return [x._to_default_unit() for x in self._energy_dict[key]]

            except KeyError:
                return None

        # Return the most recent dictionary value.
        else:
            try:
                if unit is None:
                    return self._energy_dict[key][-1]._to_default_unit().value()
                else:
                    return self._energy_dict[key][-1]._to_default_unit()

            except KeyError:
                return None

    def _getFinalFrame(self):
        """
        Get the frame from the GRO file generated at the end of the
        simulation.

        Returns
        -------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system from the final frame.
        """
        # Grab the last frame from the GRO file.
        with _Utils.cd(self._work_dir):
            # Do we need to get coordinates for the lambda=1 state.
            if "is_lambda1" in self._property_map:
                is_lambda1 = True
            else:
                is_lambda1 = False

            # Locate the coordinate file.
            if not _os.path.isfile(self._crd_file):
                _warnings.warn(
                    "Invalid coordinate file! "
                    "%s gro file not found." % (self._crd_file)
                )
                return None

            # Read the frame file.
            new_system = _IO.readMolecules(
                [self._crd_file, self._top_file], property_map=self._property_map
            )

            # Create a copy of the existing system object.
            old_system = self._system.copy()

            # Update the coordinates and velocities and return a mapping between
            # the molecule indices in the two systems.
            sire_system, mapping = _SireIO.updateCoordinatesAndVelocities(
                old_system._sire_object,
                new_system._sire_object,
                self._mapping,
                is_lambda1,
                self._property_map,
                self._property_map,
            )

            # Update the underlying Sire object.
            old_system._sire_object = sire_system

            # Store the mapping between the MolIdx in both systems so we don't
            # need to recompute it next time.
            self._mapping = mapping

            # Get the "space" property name from the property map.
            space_prop = self._property_map.get("space", "space")

            # Update the box information in the original system. Only do this if
            # the original system contains space information, since it will have
            # been added in order to run vacuum simulations.
            if (
                space_prop in old_system._sire_object.propertyKeys()
                and space_prop in new_system._sire_object.propertyKeys()
            ):
                # Get the original space.
                box = old_system._sire_object.property("space")

                # Only update the box if the space is periodic.
                if box.isPeriodic():
                    box = new_system._sire_object.property("space")
                    old_system._sire_object.setProperty(
                        self._property_map.get("space", "space"), box
                    )

            # If this is a vacuum simulation, then translate the centre of mass
            # of the system back to the origin.
            if not space_prop in old_system._sire_object.propertyKeys():
                com = old_system._getCenterOfMass()
                old_system.translate([-x for x in com])

            return old_system

    def _getFrame(self, time):
        """
        Get the trajectory frame closest to a specific time value.

        Parameters
        ----------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The time value.

        Returns
        -------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system from the closest trajectory frame.
        """

        if not isinstance(time, _Types.Time):
            raise TypeError("'time' must be of type 'BioSimSpace.Types.Time'")

        # Grab the last frame from the current trajectory file.
        try:
            with _Utils.cd(self._work_dir):
                # Do we need to get coordinates for the lambda=1 state.
                if "is_lambda1" in self._property_map:
                    is_lambda1 = True
                else:
                    is_lambda1 = False

                # Locate the trajectory file.
                traj_file = self._find_trajectory_file()

                if traj_file is None:
                    return None
                else:
                    self._traj_file = traj_file

                # Use trjconv to get the frame closest to the current simulation time.
                command = "%s trjconv -f %s -s %s -dump %f -pbc mol -o frame.gro" % (
                    self._exe,
                    self._traj_file,
                    self._tpr_file,
                    time.picoseconds().value(),
                )

                # Run the command as a pipeline.
                proc_echo = _subprocess.Popen(
                    ["echo", "0"], shell=False, stdout=_subprocess.PIPE
                )
                proc = _subprocess.Popen(
                    _Utils.command_split(command),
                    shell=False,
                    stdin=proc_echo.stdout,
                    stdout=_subprocess.PIPE,
                    stderr=_subprocess.PIPE,
                )
                proc.wait()
                proc_echo.stdout.close()

                # Read the frame file.
                new_system = _IO.readMolecules(
                    ["frame.gro", self._top_file], property_map=self._property_map
                )

                # Delete the frame file.
                _os.remove("frame.gro")

                # Create a copy of the existing system object.
                old_system = self._system.copy()

                # Update the coordinates and velocities and return a mapping between
                # the molecule indices in the two systems.
                sire_system, mapping = _SireIO.updateCoordinatesAndVelocities(
                    old_system._sire_object,
                    new_system._sire_object,
                    self._mapping,
                    is_lambda1,
                    self._property_map,
                    self._property_map,
                )

                # Update the underlying Sire object.
                old_system._sire_object = sire_system

                # Store the mapping between the MolIdx in both systems so we don't
                # need to recompute it next time.
                self._mapping = mapping

                # Get the "space" property name from the property map.
                space_prop = self._property_map.get("space", "space")

                # Update the box information in the original system. Only do this if
                # the original system contains space information, since it will have
                # been added in order to run vacuum simulations.
                if (
                    space_prop in old_system._sire_object.propertyKeys()
                    and space_prop in new_system._sire_object.propertyKeys()
                ):
                    # Get the original space.
                    box = old_system._sire_object.property("space")

                    # Only update the box if the space is periodic.
                    if box.isPeriodic():
                        box = new_system._sire_object.property("space")
                        if box.isPeriodic():
                            old_system._sire_object.setProperty(
                                self._property_map.get("space", "space"), box
                            )

                # If this is a vacuum simulation, then translate the centre of mass
                # of the system back to the origin.
                if not space_prop in old_system._sire_object.propertyKeys():
                    com = new_system._getCenterOfMass()
                    old_system.translate([-x for x in com])

                return old_system

        except:
            _warnings.warn(
                "Failed to extract trajectory frame with trjconv. "
                "Try running 'getSystem' again."
            )
            frame = "%s/frame.gro" % self._work_dir
            if _os.path.isfile(frame):
                _os.remove(frame)
            return None

    def _find_trajectory_file(self):
        """
        Helper function to find the trajectory file associated with the
        process.

        Returns
        -------

        traj_file : str
            The path to the trajectory file.
        """

        # Check that the current trajectory file is found.
        if not _os.path.isfile(self._traj_file):
            # If not, first check for any trr extension.
            traj_file = _glob.glob("%s/*.trr" % self._work_dir)

            # Store the number of trr files.
            num_trr = len(traj_file)

            # Only accept if a single trajectory file is present.
            if num_trr == 1:
                return traj_file[0]
            else:
                # Now check for any xtc files.
                traj_file = _glob.glob("%s/*.xtc" % self._work_dir)

                if len(traj_file) == 1:
                    return traj_file[0]
                else:
                    _warnings.warn(
                        "Invalid trajectory file! "
                        "%d trr files found, %d xtc files found."
                        % (num_trr, len(traj_file))
                    )
                    return None
        else:
            return self._traj_file

    def _saveMetric(
        self, filename="metric.parquet", u_nk="u_nk.parquet", dHdl="dHdl.parquet"
    ):
        """
        Helper function to save the simulation metrics to `filename`, which is a
        pandas dataframe that can be loaded with `pd.read_parquet`. if the protocol
        is Free Energy protocol, the dHdl and the u_nk data will be saved in the
        same parquet format as well.
        """

        _assert_imported(_alchemlyb)

        self._update_energy_dict(initialise=True)
        datadict_keys = [
            ("Time (ps)", _Units.Time.picosecond, "getTime"),
            (
                "PotentialEnergy (kJ/mol)",
                _Units.Energy.kj_per_mol,
                "getPotentialEnergy",
            ),
        ]
        if not isinstance(self._protocol, _Protocol.Minimisation):
            datadict_keys.extend(
                [
                    ("Volume (nm^3)", _Units.Volume.nanometer3, "getVolume"),
                    ("Pressure (bar)", _Units.Pressure.bar, "getPressure"),
                    (
                        "Temperature (kelvin)",
                        _Units.Temperature.kelvin,
                        "getTemperature",
                    ),
                ]
            )
        df = self._convert_datadict_keys(datadict_keys)
        df.to_parquet(path=f"{self.workDir()}/{filename}", index=True)
        if isinstance(self._protocol, _Protocol.FreeEnergy):
            energy = _extract(
                f"{self.workDir()}/{self._name}.xvg",
                T=self._protocol.getTemperature() / _Units.Temperature.kelvin,
            )
            with _warnings.catch_warnings():
                _warnings.filterwarnings(
                    "ignore", message="The DataFrame has column names of mixed type."
                )
                if "u_nk" in energy:
                    energy["u_nk"].to_parquet(
                        path=f"{self.workDir()}/{u_nk}", index=True
                    )
                if "dHdl" in energy:
                    energy["dHdl"].to_parquet(
                        path=f"{self.workDir()}/{dHdl}", index=True
                    )


def _is_minimisation(config):
    """
    Helper function to check whether a custom configuration
    is a minimisation.

    Parameters
    ----------

    config : [str]
        A list of configuration strings.

    Returns
    -------

    is_minimisation : bool
        Whether this is a minimisation configuration.
    """

    for line in config:
        # Convert to lower-case and remove any whitespace.
        line = line.lower().replace(" ", "")

        # Check for integrators used for minimisation.
        if (
            "integrator=steep" in line
            or "integrator=cg" in line
            or "integrator=l-bfgs" in line
        ):
            return True

    return False


def _is_vacuum(config):
    """
    Helper function to check whether a configuration corresponds to a
    vacuum simulation.

    Parameters
    ----------

    config : [str]
        A list of configuration strings.

    Returns
    -------

    is_vacuum : bool
        Whether this is (likely) a vacuum configuration.
    """

    # Join all of the config strings together.
    config = "".join(config)

    # Convert to lower-case and remove any whitespace.
    config = config.lower().replace(" ", "")

    # Check for likely indicators of a vacuum simulation.
    # These can be adapted as we think of more, or options
    # change.
    if (
        "pbc=no" in config
        and "nstlist=0" in config
        and "rlist=0" in config
        and "rvdw=0" in config
        and "rcoulomb=0" in config
    ):
        return True
    else:
        return False
