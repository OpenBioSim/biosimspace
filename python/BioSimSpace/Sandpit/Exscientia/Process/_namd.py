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

"""Functionality for running simulations using NAMD."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Namd"]

from .._Utils import _try_import

import math as _math
import os as _os

_pygtail = _try_import("pygtail")
import timeit as _timeit
import warnings as _warnings

from sire.legacy import Base as _SireBase
from sire.legacy import IO as _SireIO
from sire.legacy import Mol as _SireMol
from sire.legacy.Maths import Vector as _Vector

from .. import _isVerbose
from .._Exceptions import IncompatibleError as _IncompatibleError
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import System as _System
from ..Types._type import Type as _Type

from .. import IO as _IO
from .. import Protocol as _Protocol
from .. import Trajectory as _Trajectory
from .. import Units as _Units
from .. import _Utils

from . import _process


class Namd(_process.Process):
    """A class for running simulations using NAMD."""

    def __init__(
        self,
        system,
        protocol,
        exe=None,
        name="namd",
        work_dir=None,
        seed=None,
        property_map={},
    ):
        """
        Constructor.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.

        protocol : :class:`Protocol <BioSimSpace.Protocol>`
            The protocol for the NAMD process.

        exe : str
            The full path to the NAMD executable.

        name : str
            The name of the process.

        work_dir :
            The working directory for the process.

        seed : int
            A random number seed.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(
            system, protocol, name, work_dir, seed, property_map=property_map
        )

        # Set the package name.
        self._package_name = "NAMD"

        # This process can generate trajectory data.
        self._has_trajectory = True

        # If the path to the executable wasn't specified, then search
        # for it in $PATH.
        if exe is None:
            try:
                self._exe = _SireBase.findExe("namd2").absoluteFilePath()
            except:
                raise _MissingSoftwareError(
                    "'BioSimSpace.Process.Namd' is not supported. "
                    "Please install NAMD (http://www.ks.uiuc.edu/Research/namd)."
                ) from None
        else:
            # Make sure executable exists.
            if _os.path.isfile(exe):
                self._exe = exe
            else:
                raise IOError("NAMD executable doesn't exist: '%s'" % exe)

        # Initialise the stdout dictionary and title header.
        self._stdout_dict = _process._MultiDict()
        self._stdout_title = None

        # The names of the input files.
        self._psf_file = "%s/%s.psf" % (self._work_dir, name)
        self._top_file = "%s/%s.pdb" % (self._work_dir, name)
        self._param_file = "%s/%s.params" % (self._work_dir, name)
        self._velocity_file = None
        self._restraint_file = None

        # The name of the trajectory file.
        self._traj_file = "%s/%s_out.dcd" % (self._work_dir, name)

        # Set the path for the NAMD configuration file.
        self._config_file = "%s/%s.cfg" % (self._work_dir, name)

        # Create the list of input files.
        self._input_files = [
            self._config_file,
            self._psf_file,
            self._top_file,
            self._param_file,
        ]

        # Now set up the working directory for the process.
        self._setup()

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...

        # Create a copy of the system.
        system = self._system.copy()

        # Check for perturbable molecules and convert to the chosen end state.
        system = self._checkPerturbable(system)

        # PSF and parameter files.
        try:
            file = _os.path.splitext(self._psf_file)[0]
            _IO.saveMolecules(file, system, "psf", property_map=self._property_map)
        except Exception as e:
            msg = "Failed to write system to 'CHARMMPSF' format."
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # PDB file.
        try:
            pdb = _SireIO.PDB2(system._sire_object, self._property_map)
            pdb.writeToFile(self._top_file)
        except Exception as e:
            msg = "Failed to write system to 'PDB' format."
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Try to write a PDB "velocity" restart file.
        # The file will only be generated if all atoms in the system have
        # a "velocity" property.

        # First generate a name for the velocity file.
        velocity_file = _os.path.splitext(self._top_file)[0] + ".vel"

        # Write the velocity file.
        has_velocities = pdb.writeVelocityFile(velocity_file)

        # If a file was written, store the name of the file and update the
        # list of input files.
        if has_velocities:
            self._velocity_file = velocity_file
            self._input_files.append(self._velocity_file)

        # NAMD requires donor, acceptor, and non-bonded exclusion record entries
        # in the PSF file. We check that these are present and append blank
        # records if not. The records are actually redundant (they have no
        # affect on the MD) so could be stripped (or zeroed) by the CharmmPSF
        # parser.
        has_donors = False
        has_acceptors = False
        has_non_bonded = False

        # When converting forcefields, improper terms may be represented as
        # dihedrals (cosine impropers). As such, there will be no improper
        # records in the PSF file.
        has_impropers = False

        # Open the PSF file for reading.
        with open(self._psf_file) as file:
            # Loop over all lines.
            for line in file:
                # There are improper records.
                if "!NIMPHI" in line:
                    has_impropers = True

                # There are donor records.
                elif "!NDON" in line:
                    has_donors = True

                # There are acceptor records.
                elif "NACC" in line:
                    has_acceptors = True

                # There are non-bonded exclusion records.
                elif "NNB" in line:
                    has_non_bonded = True

        # Append empty improper record.
        if not has_impropers:
            file = open(self._psf_file, "a")
            file.write("\n%8d !NIMPHI: impropers\n" % 0)
            file.close()

        # Append empty donor record.
        if not has_donors:
            file = open(self._psf_file, "a")
            file.write("\n%8d !NDON: donors\n" % 0)
            file.close()

        # Append empty acceptor record.
        if not has_acceptors:
            file = open(self._psf_file, "a")
            file.write("\n%8d !NACC: acceptors\n" % 0)
            file.close()

        # Append empty non-bonded exclusion record.
        if not has_non_bonded:
            file = open(self._psf_file, "a")
            file.write("\n%8d !NNB: excluded\n" % 0)
            file.close()

        # Generate the NAMD configuration file.
        if isinstance(self._protocol, _Protocol.Custom):
            self.setConfig(self._protocol.getConfig())
        else:
            self._generate_config()
        self.writeConfig(self._config_file)

        # Return the list of input files.
        return self._input_files

    def _generate_config(self):
        """Generate NAMD configuration file strings."""

        # Clear the existing configuration list.
        self._config = []

        # Flag that the system doesn't contain a box.
        has_box = False

        prop = self._property_map.get("space", "space")

        # Check whether the system contains periodic box information.
        if prop in self._system._sire_object.propertyKeys():
            try:
                box = self._system._sire_object.property(prop)

                # Flag that we have found a periodic box.
                has_box = box.isPeriodic()
            except Exception:
                box = None

        # Check whether the system contains periodic box information.
        if has_box:
            # Periodic box.
            try:
                box_size = self._system._sire_object.property(prop).dimensions()
                v0 = _Vector(box_size.x().value(), 0, 0)
                v1 = _Vector(0, box_size.y().value(), 0)
                v2 = _Vector(0, 0, box_size.z().value())

            # TriclinicBox.
            except:
                box = self._system._sire_object.property(prop)
                v0 = box.vector0()
                v1 = box.vector1()
                v2 = box.vector2()

            # Work out the minimum box size.
            box_size = min(v0.magnitude(), v1.magnitude(), v2.magnitude())

            # Convert vectors to tuples.
            v0 = tuple(x.value() for x in v0)
            v1 = tuple(x.value() for x in v1)
            v2 = tuple(x.value() for x in v2)

            # Since the box is translationally invariant, we set the cell
            # origin to be the average of the atomic coordinates. This
            # ensures a consistent wrapping for coordinates in the NAMD
            # output files.
            origin = tuple(x.value() for x in self._system._getAABox().center())

        # No box information. Assume this is a gas phase simulation.
        else:
            _warnings.warn("No simulation box found. Assuming gas phase simulation.")
            has_box = False

        prop = self._property_map.get("param_format", "param_format")

        # Check whether the system contains parameter format information.
        if prop in self._system._sire_object.propertyKeys():
            # Get the parameter format.
            if self._system._sire_object.property(prop).toString() == "CHARMM":
                is_charmm = True
            else:
                is_charmm = False

        # No parameter format information. Assume these are CHARMM parameters.
        else:
            _warnings.warn(
                "No parameter format information found. Assuming CHARMM parameters."
            )
            is_charmm = True

        # Append generic configuration variables.

        # Topology.
        self.addToConfig("structure             %s" % _os.path.basename(self._psf_file))
        self.addToConfig("coordinates           %s" % _os.path.basename(self._top_file))

        # Velocities.
        if self._velocity_file is not None:
            self.addToConfig(
                "velocities            %s" % _os.path.basename(self._velocity_file)
            )

        # Parameters.
        if is_charmm:
            self.addToConfig("paraTypeCharmm        on")
        self.addToConfig(
            "parameters            %s" % _os.path.basename(self._param_file)
        )

        # Random number seed.
        if self._is_seeded:
            self.addToConfig("seed                  %d" % self._seed)

        # Exclusion policy.
        self.addToConfig("exclude               scaled1-4")

        # Non-bonded potential parameters.

        # Gas phase.
        if not has_box or not self._has_water:
            self.addToConfig("cutoff                999.")
            self.addToConfig("zeroMomentum          yes")
            self.addToConfig("switching             off")

        # Solvated.
        else:
            # Only use a cutoff if the box is large enough.
            if box_size > 26:
                self.addToConfig("cutoff                12.")
                self.addToConfig("pairlistdist          14.")
                self.addToConfig("switching             on")
                self.addToConfig("switchdist            10.")

            # Load the XSC file.
            self.addToConfig("extendedSystem        %s.xsc" % self._name)
            # We force the cell origin to be located at the system's centre
            # of geometry. This ensures a consistent periodic wrapping for
            # all NAMD output.
            self.addToConfig("cellOrigin            %.3f   %.3f   %.3f" % origin)

            # Add the cell vectors.
            self.addToConfig("cellBasisVector1      %.3f   %.3f   %.3f" % v0)
            self.addToConfig("cellBasisVector2      %.3f   %.3f   %.3f" % v1)
            self.addToConfig("cellBasisVector3      %.3f   %.3f   %.3f" % v2)

            # Wrap all molecular coordinates to the periodic box.
            self.addToConfig("wrapAll               on")

            # Periodic electrostatics.
            self.addToConfig("PME                   yes")
            self.addToConfig("PMEGridSpacing        1.")

        # Output file parameters.
        self.addToConfig("outputName            %s_out" % self._name)
        self.addToConfig("binaryOutput          no")
        self.addToConfig("binaryRestart         no")

        # Add configuration variables for a minimisation simulation.
        if isinstance(self._protocol, _Protocol.Minimisation):
            # Output frequency.
            self.addToConfig("restartfreq           500")
            self.addToConfig("xstFreq               500")

            # Printing frequency.
            self.addToConfig("outputEnergies        100")
            self.addToConfig("outputTiming          1000")

            self.addToConfig("temperature           300")

            # Work out the number of steps. This must be a multiple of
            # stepspercycle, which is set the default of 20.
            steps = 20 * _math.ceil(self._protocol.getSteps() / 20)
            self.addToConfig("minimize              %d" % steps)

        # Add configuration variables for an equilibration simulation.
        elif isinstance(self._protocol, _Protocol.Equilibration):
            # Work out the number of integration steps.
            steps = _math.ceil(
                self._protocol.getRunTime() / self._protocol.getTimeStep()
            )

            # Get the report and restart intervals.
            report_interval = self._protocol.getReportInterval()
            restart_interval = self._protocol.getRestartInterval()

            # Cap the intervals at the total number of steps.
            if report_interval > steps:
                report_interval = steps
            if restart_interval > steps:
                restart_interval = steps

            # Output frequency.
            self.addToConfig("restartfreq           %d" % restart_interval)
            self.addToConfig("xstFreq               %d" % restart_interval)

            # Printing frequency.
            self.addToConfig("outputEnergies        %d" % report_interval)
            self.addToConfig("outputTiming          1000")

            # Set the Tcl temperature variable.
            if self._protocol.isConstantTemp():
                self.addToConfig(
                    "set temperature       %.2f"
                    % self._protocol.getStartTemperature().kelvin().value()
                )
            else:
                self.addToConfig(
                    "set temperature       %.2f"
                    % self._protocol.getEndTemperature().kelvin().value()
                )
            self.addToConfig("temperature           $temperature")

            # Integrator parameters.
            self.addToConfig(
                "timestep              %.2f"
                % self._protocol.getTimeStep().femtoseconds().value()
            )
            self.addToConfig("rigidBonds            all")
            self.addToConfig("nonbondedFreq         1")
            self.addToConfig("fullElectFrequency    2")

            # Constant temperature control.
            self.addToConfig("langevin              on")
            self.addToConfig("langevinDamping       1.")
            self.addToConfig("langevinTemp          $temperature")
            self.addToConfig("langevinHydrogen      no")

            # Constant pressure control.
            if self._protocol.getPressure() is not None:
                self.addToConfig("langevinPiston        on")
                self.addToConfig(
                    "langevinPistonTarget  %.5f"
                    % self._protocol.getPressure().bar().value()
                )
                self.addToConfig("langevinPistonPeriod  100.")
                self.addToConfig("langevinPistonDecay   50.")
                self.addToConfig("langevinPistonTemp    $temperature")
                self.addToConfig("useGroupPressure      yes")
                self.addToConfig("useFlexibleCell       no")
                self.addToConfig("useConstantArea       no")

            # Restrain the backbone.
            restraint = self._protocol.getRestraint()
            if restraint is not None:
                # Create a restrained system.
                restrained = self._createRestrainedSystem(self._system, restraint)

                # Create a PDB object, mapping the "occupancy" property to "restrained".
                prop = self._property_map.get("occupancy", "occupancy")

                try:
                    p = _SireIO.PDB2(restrained._sire_object, {prop: "restrained"})

                    # File name for the restraint file.
                    self._restraint_file = "%s/%s.restrained" % (
                        self._work_dir,
                        self._name,
                    )

                    # Write the PDB file.
                    p.writeToFile(self._restraint_file)

                except:
                    _warnings.warn(
                        "Failed to add restraints to PDB file. "
                        "Perhaps there are no atoms matching the restraint?"
                    )

                # Update the configuration file.
                self.addToConfig("constraints           yes")
                self.addToConfig("consref               %s.restrained" % self._name)
                self.addToConfig("conskfile             %s.restrained" % self._name)
                self.addToConfig("conskcol              O")

            # Heating/cooling simulation.
            if not self._protocol.isConstantTemp():
                # Work out temperature step size (assuming a unit increment).
                denom = abs(
                    self._protocol.getEndTemperature().kelvin().value()
                    - self._protocol.getStartTemperature().kelvin().value()
                )
                freq = _math.floor(steps / denom)

                self.addToConfig("reassignFreq          %d" % freq)
                self.addToConfig(
                    "reassignTemp          %.2f"
                    % self._protocol.getStartTemperature().kelvin().value()
                )
                self.addToConfig("reassignIncr          1.")
                self.addToConfig(
                    "reassignHold          %.2f"
                    % self._protocol.getEndTemperature().kelvin().value()
                )

            # Trajectory output frequency.
            self.addToConfig("DCDfreq               %d" % restart_interval)

            # Run the simulation.
            self.addToConfig("run                   %d" % steps)

        # Add configuration variables for a production simulation.
        elif isinstance(self._protocol, _Protocol.Production):
            # Work out the number of integration steps.
            steps = _math.ceil(
                self._protocol.getRunTime() / self._protocol.getTimeStep()
            )

            # Get the report and restart intervals.
            report_interval = self._protocol.getReportInterval()
            restart_interval = self._protocol.getRestartInterval()

            # Cap the intervals at the total number of steps.
            if report_interval > steps:
                report_interval = steps
            if restart_interval > steps:
                restart_interval = steps

            # Output frequency.
            self.addToConfig("restartfreq           %d" % restart_interval)
            self.addToConfig("xstFreq               %d" % restart_interval)

            # Printing frequency.
            self.addToConfig("outputEnergies        %d" % report_interval)
            self.addToConfig("outputTiming          1000")

            # Set the Tcl temperature variable.
            self.addToConfig(
                "set temperature       %.2f"
                % self._protocol.getTemperature().kelvin().value()
            )
            self.addToConfig("temperature           $temperature")

            # Integrator parameters.
            self.addToConfig(
                "timestep              %.2f"
                % self._protocol.getTimeStep().femtoseconds().value()
            )
            if self._protocol.getFirstStep() != 0:
                self.addToConfig(
                    "firsttimestep         %d" % self._protocol.getFirstStep()
                )
            self.addToConfig("rigidBonds            all")
            self.addToConfig("nonbondedFreq         1")
            self.addToConfig("fullElectFrequency    2")

            # Constant temperature control.
            self.addToConfig("langevin              on")
            self.addToConfig("langevinDamping       1.")
            self.addToConfig("langevinTemp          $temperature")
            self.addToConfig("langevinHydrogen      no")

            # Constant pressure control.
            if self._protocol.getPressure() is not None:
                self.addToConfig("langevinPiston        on")
                self.addToConfig(
                    "langevinPistonTarget  %.5f"
                    % self._protocol.getPressure().bar().value()
                )
                self.addToConfig("langevinPistonPeriod  100.")
                self.addToConfig("langevinPistonDecay   50.")
                self.addToConfig("langevinPistonTemp    $temperature")
                self.addToConfig("useGroupPressure      yes")
                self.addToConfig("useFlexibleCell       no")
                self.addToConfig("useConstantArea       no")

            # Work out number of steps needed to exceed desired running time.
            steps = _math.ceil(
                self._protocol.getRunTime() / self._protocol.getTimeStep()
            )

            # Trajectory output frequency.
            self.addToConfig("DCDfreq               %d" % restart_interval)

            # Run the simulation.
            self.addToConfig("run                   %d" % steps)

        else:
            raise _IncompatibleError(
                "Unsupported protocol: '%s'" % self._protocol.__class__.__name__
            )

        # Flag that this isn't a custom protocol.
        self._protocol._setCustomised(False)

    def start(self):
        """
        Start the NAMD process.

        Returns
        -------

        process : :class:`Process.Namd <BioSimSpace.Process.Namd>`
            The process object.
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
            # Write the command-line process to a README.txt file.
            with open("README.txt", "w") as file:
                # Set the command-line string.
                self._command = "%s %s.cfg" % (self._exe, self._name)

                # Write the command to file.
                file.write("# NAMD was run with the following command:\n")
                file.write("%s\n" % self._command)

            # Start the timer.
            self._timer = _timeit.default_timer()

            # Start the simulation.
            self._process = _SireBase.Process.run(
                self._exe,
                "%s.cfg" % self._name,
                "%s.out" % self._name,
                "%s.err" % self._name,
            )

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

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        # Read the PDB coordinate file and construct a parameterised molecular
        # system using the original PSF and param files.

        has_coor = False

        # First check for final configuration.
        if _os.path.isfile("%s/%s_out.coor" % (self._work_dir, self._name)):
            coor_file = "%s/%s_out.coor" % (self._work_dir, self._name)
            has_coor = True

        # Otherwise check for a restart file.
        elif _os.path.isfile("%s/%s_out.restart.coor" % (self._work_dir, self._name)):
            coor_file = "%s/%s_out.restart.coor" % (self._work_dir, self._name)
            has_coor = True

        # Try to find an XSC file.

        has_xsc = False

        # First check for final XSC file.
        if _os.path.isfile("%s/%s_out.xsc" % (self._work_dir, self._name)):
            xsc_file = "%s/%s_out.xsc" % (self._work_dir, self._name)
            has_xsc = True

        # Otherwise check for a restart XSC file.
        elif _os.path.isfile("%s/%s_out.restart.xsc" % (self._work_dir, self._name)):
            xsc_file = "%s/%s_out.restart.xsc" % (self._work_dir, self._name)
            has_xsc = True

        # We found a coordinate file.
        if has_coor:
            # List of files.
            files = [coor_file, self._psf_file, self._param_file]

            # Add the box information.
            if has_xsc:
                files.append(xsc_file)

            # Create and return the molecular system.
            try:
                # Do we need to get coordinates for the lambda=1 state.
                if "is_lambda1" in self._property_map:
                    is_lambda1 = True
                else:
                    is_lambda1 = False

                # Load the restart file.
                new_system = _System(
                    _SireIO.MoleculeParser.read(files, self._property_map)
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

                # Update the box information in the original system.
                if "space" in new_system._sire_object.propertyKeys():
                    box = new_system._sire_object.property("space")
                    if box.isPeriodic():
                        old_system._sire_object.setProperty(
                            self._property_map.get("space", "space"), box
                        )

                return old_system

            except:
                return None

        else:
            return None

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

        trajectory : :class:`Trajectory <BioSimSpace.Trajectory>`
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
            # Do we need to get coordinates for the lambda=1 state.
            if "is_lambda1" in self._property_map:
                is_lambda1 = True
            else:
                is_lambda1 = False

            # Load the latest trajectory frame.
            new_system = _Trajectory.getFrame(self._traj_file, self._top_file, index)

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

            # Update the box information in the original system.
            if "space" in new_system._sire_object.propertyKeys():
                box = new_system._sire_object.property("space")
                old_system._sire_object.setProperty(
                    self._property_map.get("space", "space"), box
                )

            return old_system

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

        self.stdout(0)
        return self._get_stdout_record(record, time_series, unit)

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

        self.stdout(0)
        return self._get_stdout_record(record, time_series, unit)

    def getRecords(self, block="AUTO"):
        """
        Return the dictionary of stdout time-series records.

        Parameters
        ----------

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        records : dict
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

        self.stdout(0)
        return self._stdout_dict.copy()

    def getCurrentRecords(self):
        """
        Return the current dictionary of stdout time-series records.

        Parameters
        ----------

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        records : BioSimSpace.Process._process._MultiDict
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

        if isinstance(self._protocol, _Protocol.Minimisation):
            return None

        else:
            # Get the list of time steps.
            time_steps = self.getRecord("TS", time_series, None, block)

            # Convert the time step to the default unit.
            timestep = self._protocol.getTimeStep()._to_default_unit()

            # Multiply by the integration time step.
            if time_steps is not None:
                if time_series:
                    return [x * timestep for x in time_steps]
                else:
                    return timestep * time_steps

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
        """
        return self.getRecord("TS", time_series, None, block)

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
        return self.getRecord("BOND", time_series, _Units.Energy.kcal_per_mol, block)

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
        return self.getRecord("ANGLE", time_series, _Units.Energy.kcal_per_mol, block)

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
        proper = self.getRecord("DIHED", time_series, _Units.Energy.kcal_per_mol, block)
        improper = self.getRecord(
            "IMPRP", time_series, _Units.Energy.kcal_per_mol, block
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
        return self.getRecord("DIHED", time_series, _Units.Energy.kcal_per_mol, block)

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
        return self.getRecord("IMPRP", time_series, _Units.Energy.kcal_per_mol, block)

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

    def getElectrostaticEnergy(self, time_series=False, block="AUTO"):
        """
        Get the electrostatic energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The electrostatic energy.
        """
        return self.getRecord("ELECT", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentElectrostaticEnergy(self, time_series=False):
        """
        Get the current electrostatic energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The electrostatic energy.
        """
        return self.getElectrostaticEnergy(time_series, block=False)

    def getVanDerWaalsEnergy(self, time_series=False, block="AUTO"):
        """
        Get the Van der Vaals energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The Van der Vaals energy.
        """
        return self.getRecord("VDW", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentVanDerWaalsEnergy(self, time_series=False):
        """
        Get the current Van der Waals energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The Van der Vaals energy.
        """
        return self.getVanDerWaalsEnergy(time_series, block=False)

    def getBoundaryEnergy(self, time_series=False, block="AUTO"):
        """
        Get the boundary energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The boundary energy.
        """
        return self.getRecord(
            "BOUNDARY", time_series, _Units.Energy.kcal_per_mol, block
        )

    def getCurrentBoundaryEnergy(self, time_series=False):
        """
        Get the current boundary energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The boundary energy.
        """
        return self.getBoundaryEnergy(time_series, block=False)

    def getMiscEnergy(self, time_series=False, block="AUTO"):
        """
        Get the external energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The external energy.
        """
        return self.getRecord("MISC", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentMiscEnergy(self, time_series=False):
        """
        Get the current external energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The external energy.
        """
        return self.getMiscEnergy(time_series, block=False)

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
        return self.getRecord("KINETIC", time_series, _Units.Energy.kcal_per_mol, block)

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
            The current kinetic energy.
        """
        return self.getKineticEnergy(time_series, block=False)

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
        return self.getRecord(
            "POTENTIAL", time_series, _Units.Energy.kcal_per_mol, block
        )

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
        return self.getRecord("TOTAL", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentTotalEnergy(self, time_series=False):
        """
        Get the current potential energy.

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

    def getTotal2Energy(self, time_series=False, block="AUTO"):
        """
        Get the total energy. (Better KE conservation.).

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
        return self.getRecord("TOTAL2", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentTotal2Energy(self, time_series=False):
        """
        Get the current total energy. (Better KE conservation.).

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The total energy.
        """
        return self.getTotal2Energy(time_series, block=False)

    def getTotal3Energy(self, time_series=False, block="AUTO"):
        """
        Get the total energy. (Smaller short-time fluctuations.).

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
        return self.getRecord("TOTAL3", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentTotal3Energy(self, time_series=False):
        """
        Get the total energy. (Smaller short-time fluctuations.).

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The total energy.
        """
        return self.getTotal3Energy(time_series, block=False)

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
        return self.getRecord("TEMP", time_series, _Units.Temperature.kelvin, block)

    def getCurrentTemperature(self, time_series=False):
        """
        Get the temperature.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The temperature.
        """
        return self.getTemperature(time_series, block=False)

    def getTemperatureAverage(self, time_series=False, block="AUTO"):
        """
        Get the average temperature.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The average temperature.
        """
        return self.getRecord("TEMPAVG", time_series, _Units.Temperature.kelvin, block)

    def getCurrentTemperatureAverage(self, time_series=False):
        """
        Get the current average temperature.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The average temperature.
        """
        return self.getTemperatureAverage(time_series, block=False)

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
            The pressure.
        """
        return self.getPressure(time_series, block=False)

    def getPressureAverage(self, time_series=False, block="AUTO"):
        """
        Get the average pressure.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The average pressure.
        """
        return self.getRecord("PRESSAVG", time_series, _Units.Pressure.bar, block)

    def getCurrentPressureAverage(self, time_series=False):
        """
        Get the current average pressure.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The average pressure.
        """
        return self.getPressureAverage(time_series, block=False)

    def getGPressure(self, time_series=False, block="AUTO"):
        """
        Get the pressure. (Hydrogens incorporated into bonded atoms.).

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
        return self.getRecord("GPRESSURE", time_series, _Units.Pressure.bar, block)

    def getCurrentGPressure(self, time_series=False):
        """
        Get the current pressure. (Hydrogens incorporated into bonded atoms.).

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure.
        """
        return self.getGPressure(time_series, block=False)

    def getGPressureAverage(self, time_series=False, block="AUTO"):
        """
        Get the average pressure. (Hydrogens incorporated into bonded atoms.).

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The average pressure.
        """
        return self.getRecord("GPRESSAVG", time_series, _Units.Pressure.bar, block)

    def getCurrentGPressureAverage(self, time_series=False):
        """
        Get the current average pressure. (Hydrogens incorporated into bonded atoms.).

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The average pressure.
        """
        return self.getGPressureAverage(time_series, block=False)

    def getVolume(self, time_series=False, block="AUTO"):
        """
        Get the volume.

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
        """
        Get the current volume.

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

    def eta(self):
        """
        Get the estimated time for the process to finish (in minutes).

        Returns
        -------

        eta : :class:`Time <BioSimSpace.Types.Time>`
            The estimated remaining time in minutes.
        """

        # Make sure the list of stdout records is up to date.
        # Print the last zero lines, i.e. no output.
        self.stdout(0)

        # Now search backwards through the list to find the last TIMING record.
        for _, record in reversed(list(enumerate(self._stdout))):
            # Split the record using whitespace.
            data = record.split()

            # We've found a TIMING record.
            if len(data) > 0:
                if data[0] == "TIMING:":
                    # Try to find the "hours" record.
                    # If found, return the entry precedeing it.
                    try:
                        return (
                            float(data[data.index("hours") - 1]) * 60
                        ) * _Units.Time.minutes

                    # No record found.
                    except ValueError:
                        return None

    def stdout(self, n=10):
        """
        Print the last n lines of the stdout buffer.

        Parameters
        ----------

        n : int
            The number of lines to print.
        """

        # Ensure that the number of lines is positive.
        if n < 0:
            raise ValueError("The number of lines must be positive!")

        # Append any new lines to the stdout list.
        for line in _pygtail.Pygtail(self._stdout_file):
            self._stdout.append(line.rstrip())

            # Split the record using whitespace.
            data = self._stdout[-1].split()

            # Make sure there is at least one record.
            if len(data) > 0:
                # Store the updated energy title.
                if data[0] == "ETITLE:":
                    self._stdout_title = data[1:]

                # This is an energy record.
                elif data[0] == "ENERGY:":
                    # Extract the data.
                    stdout_data = data[1:]

                    # Add the records to the dictionary.
                    if len(stdout_data) == len(self._stdout_title):
                        for title, data in zip(self._stdout_title, stdout_data):
                            self._stdout_dict[title] = data

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

    def _createRestrainedSystem(self, system, restraint):
        """
        Restrain protein backbone atoms.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.

        restraint : str, [int]
            The type of restraint.

        Returns
        -------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system with an added 'restrained' property.
        """

        # Get the force constant value in the default units. This is
        # the same as used by NAMD, i.e. kcal_per_mol/angstrom**2
        force_constant = self._protocol.getForceConstant().value()

        # Copy the original system.
        s = system.copy()

        # Keyword restraint.
        if isinstance(restraint, str):
            # Loop over all molecules by number.
            for x, mol in enumerate(s):
                # Get the indices of the restrained atoms for this molecule.
                atoms = s.getRestraintAtoms(
                    restraint, x, is_absolute=False, allow_zero_matches=True
                )

                # Extract the molecule and make it editable.
                edit_mol = mol._sire_object.edit()

                # First set all restraint force constants to 0, i.e. the restraint
                # will be ignored.
                for atom in edit_mol.atoms():
                    edit_mol = (
                        edit_mol.atom(atom.index())
                        .setProperty("restrained", 0.0)
                        .molecule()
                    )

                # Now apply restraints to the selected atoms.
                for idx in atoms:
                    edit_mol = (
                        edit_mol.atom(_SireMol.AtomIdx(idx))
                        .setProperty("restrained", 10.0)
                        .molecule()
                    )

                # Update the system.
                s._sire_object.update(edit_mol.commit())

        # A user-defined list of atoms.
        elif isinstance(restraint, (list, tuple)):
            # Create an empty multi dict for each MolNum.
            mol_atoms = {}
            for num in s._mol_nums:
                mol_atoms[num] = []

            # Now work out which MolNum corresponds to each atom in the restraint.
            for idx in restraint:
                try:
                    mol_idx, atom_idx = s._getRelativeIndices(idx)
                    mol_num = s._mol_nums[mol_idx]
                    atom_idx = _SireMol.AtomIdx(atom_idx)
                    mol_atoms[mol_num].append(atom_idx)
                except Exception as e:
                    msg = "Unable to find restrained atom in the system?"
                    if _isVerbose():
                        raise ValueError(msg) from e
                    else:
                        raise ValueError(msg) from None

            # Now loop over the multi-dict.
            for num, idxs in mol_atoms.items():
                # Extract the molecule and make it editable.
                edit_mol = s._sire_object[num].edit()

                # First set all restraints to zero.
                for atom in edit_mol.atoms():
                    edit_mol = (
                        edit_mol.atom(atom.index())
                        .setProperty("restrained", 0.0)
                        .molecule()
                    )

                # Now apply restraints to the selected atoms.
                # TODO: The fixed-width PDB format means that the force
                # constant can't exceed 999 kcal_per_mol/angstrom**2 when
                # written in the standard 2dp floating point format. We
                # could warn when the value is too large, or write as an
                # integer instead. (This latter would require tweaking the
                # PDB parser.
                for idx in idxs:
                    edit_mol = (
                        edit_mol.atom(idx)
                        .setProperty("restrained", force_constant)
                        .molecule()
                    )

                # Update the system.
                s._sire_object.update(edit_mol.commit())

        # Return the new system.
        return s

    def _get_stdout_record(self, key, time_series=False, unit=None):
        """
        Helper function to get a stdout record from the dictionary.

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
        if len(self._stdout_dict) == 0:
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
                if key == "TS":
                    return [int(x) for x in self._stdout_dict[key]]
                else:
                    if unit is None:
                        return [float(x) for x in self._stdout_dict[key]]
                    else:
                        return [
                            (float(x) * unit)._to_default_unit()
                            for x in self._stdout_dict[key]
                        ]

            except KeyError:
                return None

        # Return the most recent dictionary value.
        else:
            try:
                if key == "TS":
                    return int(self._stdout_dict[key][-1])
                else:
                    if unit is None:
                        return float(self._stdout_dict[key][-1])
                    else:
                        return (
                            float(self._stdout_dict[key][-1]) * unit
                        )._to_default_unit()

            except KeyError:
                return None
