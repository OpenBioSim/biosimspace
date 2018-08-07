######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2018
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
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
Functionality for running simulations using NAMD.
Author: Lester Hedges <lester.hedges@gmail.com>
"""

import Sire as _Sire

from . import _process
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import System as _System
from ..Trajectory import Trajectory as _Trajectory

import BioSimSpace.Protocol as _Protocol
import BioSimSpace.Types._type as _Type
import BioSimSpace.Units as _Units

import math as _math
import os as _os
import pygtail as _pygtail
import timeit as _timeit
import warnings as _warnings

__all__ = ["Namd"]

class Namd(_process.Process):
    """A class for running simulations using NAMD."""

    def __init__(self, system, protocol, exe=None,
            name="namd", work_dir=None, seed=None, map={}):
        """Constructor.

           Positional arguments
           --------------------

           system : BioSimSpace._SireWrappers.System
               The molecular system.

           protocol : BioSimSpace.Protocol
               The protocol for the AMBER process.


           Keyword arguments
           -----------------

           exe : str
               The full path to the AMBER executable.

           name : str
               The name of the process.

           work_dir :
               The working directory for the process.

           seed : int
               A random number seed.

           map : dict
               A dictionary that maps system "properties" to their user defined
               values. This allows the user to refer to properties with their
               own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(system, protocol, name, work_dir, seed, map)

        # Set the package name.
        self._package_name = "NAMD"

        # This process can generate trajectory data.
        self._has_trajectory = True

        # If the path to the executable wasn't specified, then search
        # for it in $PATH.
        if exe is None:
            try:
                self._exe = _Sire.Base.findExe("namd2").absoluteFilePath()
            except:
                raise _MissingSoftwareError("'BioSimSpace.Process.Namd' is not supported. "
                    + "Please install NAMD (http://www.ks.uiuc.edu/Research/namd).") from None
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

        # Set the path for the NAMD configuration file.
        self._config_file = "%s/%s.namd" % (self._work_dir, name)

        # Create the list of input files.
        self._input_files = [self._config_file, self._psf_file, self._top_file, self._param_file]

        # Now set up the working directory for the process.
        self._setup()

    def _setup(self):
        """Setup the input files and working directory ready for simulation."""

        # Create the input files...

        # PSF and parameter files.
        try:
            psf = _Sire.IO.CharmmPSF(self._system)
            psf.writeToFile(self._psf_file)
        except:
            raise IOError("Failed to write system to 'CHARMMPSF' format.") from None

        # PDB file.
        try:
            pdb = _Sire.IO.PDB2(self._system)
            pdb.writeToFile(self._top_file)
        except:
            raise IOError("Failed to write system to 'PDB' format.") from None

        # Try to write a PDB "velocity" restart file.
        # The file will only be generated if all atoms in self._system have
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
                if "!NDON" in line:
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
        if not has_acceptors:
            file = open(self._psf_file, "a")
            file.write("\n%8d !NNB: excluded\n" % 0)
            file.close()

        # Generate the NAMD configuration file.
        if type(self._protocol) is _Protocol.Custom:
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

        if "space" in self._map:
            prop = self._map["space"]
        else:
            prop = "space"

        # Check whether the system contains periodic box information.
        if prop in self._system.propertyKeys():
            # Flag that we have found a box.
            has_box = True

            # Get the box size.
            box_size = self._system.property(prop).dimensions()

            # Since the box is translationally invariant, we set the cell
            # origin to be the average of the atomic coordinates. This
            # ensures a consistent wrapping for coordinates in the  NAMD
            # output files.
            origin = tuple(_process._getAABox(self._system).center())

        # No box information. Assume this is a gas phase simulation.
        else:
            _warnings.warn("No simulation box found. Assuming gas phase simulation.")
            has_box = False

        if "param_format" in self._map:
            prop = self._map["param_format"]
        else:
            prop = "param_format"

        # Check whether the system contains parameter format information.
        if prop in self._system.propertyKeys():
            # Get the parameter format.
            if self._system.property(prop).toString() == "CHARMM":
                is_charmm = True
            else:
                is_charmm = False

        # No parameter format information. Assume these are CHARMM parameters.
        else:
            _warnings.warn("No parameter format information found. Assuming CHARMM parameters.")
            is_charmm = True

        # Append generic configuration variables.

        # Topology.
        self.addToConfig("structure             %s" % _os.path.basename(self._psf_file))
        self.addToConfig("coordinates           %s" % _os.path.basename(self._top_file))

        # Velocities.
        if self._velocity_file is not None:
            self.addToConfig("velocities            %s" % _os.path.basename(self._velocity_file))

        # Parameters.
        if is_charmm:
            self.addToConfig("paraTypeCharmm        on")
        self.addToConfig("parameters            %s" % _os.path.basename(self._param_file))

        # Random number seed.
        if self._is_seeded:
            self.addToConfig("seed                  %d" % self._seed)

        # Exclusion policy.
        self.addToConfig("exclude               scaled1-4")

        # Non-bonded potential parameters.

        # Gas phase.
        if not has_box:
            self.addToConfig("cutoff                999.")
            self.addToConfig("zeroMomentum          yes")
            self.addToConfig("switching             off")

        # Solvated.
        else:
            # Only use a cutoff if the box is large enough.
            if min(box_size) > 26:
                self.addToConfig("cutoff                12.")
                self.addToConfig("pairlistdist          14.")
                self.addToConfig("switching             on")
                self.addToConfig("switchdist            10.")

            # Load the XSC file.
            self.addToConfig("extendedSystem        %s.xsc" % self._name)
            # We force the cell origin to be located at the system's centre
            # of geometry. This ensures a consistent periodic wrapping for
            # all NAMD output.
            self.addToConfig("cellOrigin            %.1f   %.1f   %.1f" % origin)

            # Wrap all molecular coordinates to the periodic box.
            self.addToConfig("wrapAll               on")

            # Periodic electrostatics.
            self.addToConfig("PME                   yes")
            self.addToConfig("PMEGridSpacing        1.")

        # Output file parameters.
        self.addToConfig("outputName            %s_out" % self._name)
        self.addToConfig("binaryOutput          no")
        self.addToConfig("binaryRestart         no")

        # Output frequency.
        self.addToConfig("restartfreq           500")
        self.addToConfig("xstFreq               500")

        # Printing frequency.
        self.addToConfig("outputEnergies        100")
        self.addToConfig("outputTiming          1000")

        # Add configuration variables for a minimisation simulation.
        if type(self._protocol) is _Protocol.Minimisation:
            self.addToConfig("temperature           300")

            # Work out the number of steps. This must be a multiple of
            # stepspercycle, which is set the default of 20.
            steps = 20 * _math.ceil(self._protocol.getSteps() / 20)
            self.addToConfig("minimize              %d" % steps)

        # Add configuration variables for an equilibration simulation.
        if type(self._protocol) is _Protocol.Equilibration:
            # Set the Tcl temperature variable.
            if self._protocol.isConstantTemp():
                self.addToConfig("set temperature       %.2f"
                    % self._protocol.getStartTemperature().kelvin().magnitude())
            else:
                self.addToConfig("set temperature       %.2f"
                    % self._protocol.getEndTemperature().kelvin().magnitude())
            self.addToConfig("temperature           $temperature")

            # Integrator parameters.
            self.addToConfig("timestep              %.2f"
                % self._protocol.getTimeStep().femtoseconds().magnitude())
            self.addToConfig("rigidBonds            all")
            self.addToConfig("nonbondedFreq         1")
            self.addToConfig("fullElectFrequency    2")

            # Constant temperature control.
            self.addToConfig("langevin              on")
            self.addToConfig("langevinDamping       1.")
            self.addToConfig("langevinTemp          $temperature")
            self.addToConfig("langevinHydrogen      no")

            # Constant pressure control.
            if not self._protocol.isConstantTemp():
                self.addToConfig("langevinPiston        on")
                self.addToConfig("langevinPistonTarget  1.01325")
                self.addToConfig("langevinPistonPeriod  100.")
                self.addToConfig("langevinPistonDecay   50.")
                self.addToConfig("langevinPistonTemp    $temperature")
                self.addToConfig("useGroupPressure      yes")
                self.addToConfig("useFlexibleCell       no")
                self.addToConfig("useConstantArea       no")

            # Restrain the backbone.
            if self._protocol.isRestrained():
                # Create a restrained system.
                restrained = _process._restrain_backbone(self._system)

                # Create a PDB object, mapping the "occupancy" property to "restrained".
                if "occupancy" in self._map:
                    prop = self._map["occupancy"]
                else:
                    prop = "occupancy"

                try:
                    p = PDB2(restrained, {prop : "restrained"})

                    # File name for the restraint file.
                    self._restraint_file = "%s/%s.restrained" % (self._work_dir, self._name)

                    # Write the PDB file.
                    p.writeToFile(self._restraint_file)
                except:
                    raise IOError("Failed to add restraints to 'PDB' file.") from None

                # Update the configuration file.
                self.addToConfig("fixedAtoms            yes")
                self.addToConfig("fixedAtomsFile        %s.restrained" % self._name)

            # Work out number of steps needed to exceed desired running time,
            # rounded up to the nearest 20.
            steps = _math.ceil(self._protocol.getRunTime().nanoseconds().magnitude() /
                               self._protocol.getTimeStep().nanoseconds().magnitude())
            steps = 20 * _math.ceil(steps / 20)

            # Heating/cooling simulation.
            if not self._protocol.isConstantTemp():
                # Work out temperature step size (assuming a unit increment).
                denom = abs(self._protocol.getEndTemperature().kelvin().magnitude() -
                            self._protocol.getStartTemperature().kelvin().magnitude())
                freq = _math.floor(steps / denom)

                self.addToConfig("reassignFreq          %d" % freq)
                self.addToConfig("reassignTemp          %.2f"
                    % self._protocol.getStartTemperature().kelvin().magnitude())
                self.addToConfig("reassignIncr          1.")
                self.addToConfig("reassignHold          %.2f"
                    % self._protocol.getEndTemperature().kelvin().magnitude())

            # Trajectory output frequency.
            self.addToConfig("DCDfreq               %d" % _math.floor(steps / self._protocol.getFrames()))

            # Run the simulation.
            self.addToConfig("run                   %d" % steps)

        # Add configuration variables for a production simulation.
        elif type(self._protocol) is _Protocol.Production:
            # Set the Tcl temperature variable.
            self.addToConfig("set temperature       %.2f"
                % self._protocol.getTemperature().kelvin().magnitude())
            self.addToConfig("temperature           $temperature")

            # Integrator parameters.
            self.addToConfig("timestep              %.2f"
                % self._protocol.getTimeStep().femtoseconds().magnitude())
            if self._protocol.getFirstStep() != 0:
                self.addToConfig("firsttimestep         %d" % self._protocol.getFirstStep())
            self.addToConfig("rigidBonds            all")
            self.addToConfig("nonbondedFreq         1")
            self.addToConfig("fullElectFrequency    2")

            # Constant temperature control.
            self.addToConfig("langevin              on")
            self.addToConfig("langevinDamping       1.")
            self.addToConfig("langevinTemp          $temperature")
            self.addToConfig("langevinHydrogen      no")

            # Constant pressure control.
            if self._protocol.getEnsemble() == "NPT":
                self.addToConfig("langevinPiston        on")
                self.addToConfig("langevinPistonTarget  1.01325")
                self.addToConfig("langevinPistonPeriod  100.")
                self.addToConfig("langevinPistonDecay   50.")
                self.addToConfig("langevinPistonTemp    $temperature")
                self.addToConfig("useGroupPressure      yes")
                self.addToConfig("useFlexibleCell       no")
                self.addToConfig("useConstantArea       no")

            # Work out number of steps needed to exceed desired running time,
            # rounded up to the nearest 20.
            steps = _math.ceil(self._protocol.getRunTime().nanoseconds().magnitude() /
                               self._protocol.getTimeStep().nanoseconds().magnitude())
            steps = 20 * _math.ceil(steps / 20)

            # Trajectory output frequency.
            self.addToConfig("DCDfreq               %d" % _math.floor(steps / self._protocol.getFrames()))

            # Run the simulation.
            self.addToConfig("run                   %d" % steps)

    def start(self):
        """Start the NAMD process.

           Returns
           -------

           process : BioSimSpace.Process.Namd
               The process object.
        """

        # Process is already running.
        if self._process is not None:
            if self._process.isRunning():
                return

        # Clear any existing output.
        self._clear_output()

        # Store the current working directory.
        dir = _os.getcwd()

        # Change to the working directory for the process.
        # This avoids problems with relative paths.
        _os.chdir(self._work_dir)

        # Write the command-line process to a README.txt file.
        with open("README.txt", "w") as file:

            # Set the command-line string.
            self._command = "%s %s.namd" % (self._exe, self._name)

            # Write the command to file.
            file.write("# NAMD was run with the following command:\n")
            file.write("%s\n" % self._command)

        # Start the timer.
        self._timer = _timeit.default_timer()

        # Start the simulation.
        self._process = _Sire.Base.Process.run(self._exe,
            "%s.namd" % self._name, "%s.out" % self._name, "%s.err" % self._name)

        # Change back to the original working directory.
        _os.chdir(dir)

        return self

    def getSystem(self, block="AUTO"):
        """Get the latest molecular configuration as a Sire system.

           Keyword arguments
           -----------------

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           system : BioSimSpace._SireWrappers.System
               The latest molecular system.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

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
            files = [ coor_file, self._psf_file, self._param_file ]

            # Add the box information.
            if has_xsc:
                files.append(xsc_file)

            # Create and return the molecular system.
            try:
                return _System(_Sire.IO.MoleculeParser.read(files))
            except:
                return None

        else:
            return None

    def getCurrentSystem(self):
        """Get the latest molecular configuration as a Sire system.

           Returns
           -------

           system : BioSimSpace._SireWrappers.System
               The latest molecular system.
        """
        return self.getSystem(block=False)

    def getTrajectory(self, block="AUTO"):
        """Return a trajectory object.

           Keyword arguments
           -----------------

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           trajectory : mdtraj.Trajectory, MDAnalysis.Universe
        """
        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        return _Trajectory(process=self)

    def getRecord(self, record, time_series=False, unit=None, block="AUTO"):
        """Get a record from the stdout dictionary.

           Positional arguments
           --------------------

           record : str
               The record key.


           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           unit : BioSimSpace.Types.Type
               The unit to convert the record to.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           record :
               The matching record.
        """

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        self.stdout(0)
        return self._get_stdout_record(record, time_series, unit)

    def getCurrentRecord(self, record, time_series=False, unit=None):
        """Get a current record from the stdout dictionary.

           Positional arguments
           --------------------

           record : str
               The record key.


           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           unit : BioSimSpace.Types.Type
               The unit to convert the record to.

           Returns
           -------

           record :
               The matching record.
        """
        self.stdout(0)
        return self._get_stdout_record(record, time_series, unit)

    def getRecords(self, block="AUTO"):
        """Return the dictionary of stdout time-series records.

           Keyword arguments:

           block       -- Whether to block until the process has finished running.
        """
        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        return self._stdout_dict.copy()

    def getCurrentRecords(self):
        """Return the current dictionary of stdout time-series records.

           Keyword arguments
           -----------------

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           records : BioSimSpace.Process._process._MultiDict
              The dictionary of time-series records.
        """
        return getRecords(block=False)

    def getTime(self, time_series=False, block="AUTO"):
        """Get the time (in nanoseconds).

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           time : BioSimSpace.Types.Time
               The current simulation time in nanoseconds.
        """

        if type(self._protocol) is _Protocol.Minimisation:
            return None

        else:
            # Get the list of time steps.
            time_steps = self.getRecord("TS", time_series, None, block)

            # Convert the time step to nanoseconds.
            timestep = self._protocol.getTimeStep().nanoseconds()

            # Multiply by the integration time step.
            if time_steps is not None:
                if time_series:
                    return [x * timestep for x in time_steps]
                else:
                    return timestep * time_steps

    def getCurrentTime(self, time_series=False):
        """Get the current time (in nanoseconds).

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           time : BioSimSpace.Types.Time
               The current simulation time in nanoseconds.
        """
        return self.getTime(time_series, block=False)

    def getStep(self, time_series=False, block="AUTO"):
        """Get the number of integration steps.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           step : int
               The current number of integration steps.
        """
        return self.getRecord("TS", time_series, block)

    def getCurrentStep(self, time_series=False):
        """Get the current number of integration steps.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           step : int
               The current number of integration steps.
        """
        return self.getStep(time_series, block=False)

    def getBondEnergy(self, time_series=False, block="AUTO"):
        """Get the bond energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The bond energy.
        """
        return self.getRecord("BOND", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentBondEnergy(self, time_series=False):
        """Get the current bond energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The bond energy.
        """
        return self.getBondEnergy(time_series, block=False)

    def getAngleEnergy(self, time_series=False, block="AUTO"):
        """Get the angle energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The angle energy.
        """
        return self.getRecord("ANGLE", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentAngleEnergy(self, time_series=False):
        """Get the current angle energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The angle energy.
        """
        return self.getAngleEnergy(time_series, block=False)

    def getDihedralEnergy(self, time_series=False, block="AUTO"):
        """Get the dihedral energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The dihedral energy.
        """
        return self.getRecord("DIHED", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentDihedralEnergy(self, time_series=False):
        """Get the current dihedral energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The dihedral energy.
        """
        return self.getDihedralEnergy(time_series, block=False)

    def getImproperEnergy(self, time_series=False, block="AUTO"):
        """Get the improper energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The improper energy.
        """
        return self.getRecord("IMPRP", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentImproperEnergy(self, time_series=False):
        """Get the current improper energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The improper energy.
        """
        return self.getImproperEnergy(time_series, block=False)

    def getElectrostaticEnergy(self, time_series=False, block="AUTO"):
        """Get the electrostatic energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The electrostatic energy.
        """
        return self.getRecord("ELECT", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentElectrostaticEnergy(self, time_series=False):
        """Get the current electrostatic energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The electrostatic energy.
        """
        return self.getElectrostaticEnergy(time_series, block=False)

    def getVanDerWaalsEnergy(self, time_series=False, block="AUTO"):
        """Get the Van der Vaals energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The Van der Vaals energy.
        """
        return self.getRecord("VDW", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentVanDerWaalsEnergy(self, time_series=False):
        """Get the current Van der Waals energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The Van der Vaals energy.
        """
        return self.getVanDerWaalsEnergy(time_series, block=False)

    def getBoundaryEnergy(self, time_series=False, block="AUTO"):
        """Get the boundary energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The boundary energy.
        """
        return self.getRecord("BOUNDARY", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentBoundaryEnergy(self, time_series=False):
        """Get the current boundary energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The boundary energy.
        """
        return self.getBoundaryEnergy(time_series, block=False)

    def getMiscEnergy(self, time_series=False, block="AUTO"):
        """Get the external energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The external energy.
        """
        return self.getRecord("MISC", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentMiscEnergy(self, time_series=False):
        """Get the current external energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The external energy.
        """
        return self.getMiscEnergy(time_series, block=False)

    def getKineticEnergy(self, time_series=False, block="AUTO"):
        """Get the kinetic energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
           block       -- Whether to block until the process has finished running.
        """
        return self.getRecord("KINETIC", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentKineticEnergy(self, time_series=False):
        """Get the current kinetic energy.

           Keyword arguments:

           time_series -- Whether to return a list of time series records.
        """
        return self.getKineticEnergy(time_series, block=False)

    def getPotentialEnergy(self, time_series=False, block="AUTO"):
        """Get the potential energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The potential energy.
        """
        return self.getRecord("POTENTIAL", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentPotentialEnergy(self, time_series=False):
        """Get the current potential energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The potential energy.
        """
        return self.getPotentialEnergy(time_series, block=False)

    def getTotalEnergy(self, time_series=False, block="AUTO"):
        """Get the total energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The total energy.
        """
        return self.getRecord("TOTAL", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentTotalEnergy(self, time_series=False):
        """Get the current potential energy.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The total energy.
        """
        return self.getTotalEnergy(time_series, block=False)

    def getTotal2Energy(self, time_series=False, block="AUTO"):
        """Get the total energy. (Better KE conservation.)

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The total energy.
        """
        return self.getRecord("TOTAL2", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentTotal2Energy(self, time_series=False):
        """Get the current total energy. (Better KE conservation.)

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The total energy.
        """
        return self.getTotal2Energy(time_series, block=False)

    def getTotal3Energy(self, time_series=False, block="AUTO"):
        """Get the total energy. (Smaller short-time fluctuations.)

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The total energy.
        """
        return self.getRecord("TOTAL3", time_series, _Units.Energy.kcal_per_mol, block)

    def getCurrentTotal3Energy(self, time_series=False):
        """Get the total energy. (Smaller short-time fluctuations.)

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           energy : BioSimSpace.Types.Energy
               The total energy.
        """
        return self.getTotal3Energy(time_series, block=False)

    def getTemperature(self, time_series=False, block="AUTO"):
        """Get the temperature.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           temperature : BioSimSpace.Types.Temperature
               The temperature.
        """
        return self.getRecord("TEMP", time_series, _Units.Temperature.kelvin, block)

    def getCurrentTemperature(self, time_series=False):
        """Get the temperature.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           temperature : BioSimSpace.Types.Temperature
               The temperature.
        """
        return self.getTemperature(time_series, block=False)

    def getTemperatureAverage(self, time_series=False, block="AUTO"):
        """Get the average temperature.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           temperature : BioSimSpace.Types.Temperature
               The average temperature.
        """
        return self.getRecord("TEMPAVG", time_series, _Units.Temperature.kelvin, block)

    def getCurrentTemperatureAverage(self, time_series=False):
        """Get the current average temperature.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           temperature : BioSimSpace.Types.Temperature
               The average temperature.
        """
        return self.getTemperatureAverage(time_series, block=False)

    def getPressure(self, time_series=False, block="AUTO"):
        """Get the pressure.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           pressure : BioSimSpace.Types.Pressure
               The pressure.
        """
        return self.getRecord("PRESSURE", time_series, _Units.Pressure.bar, block)

    def getCurrentPressure(self, time_series=False):
        """Get the current pressure.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           pressure : BioSimSpace.Types.Pressure
               The pressure.
        """
        return self.getPressure(time_series, block=False)

    def getPressureAverage(self, time_series=False, block="AUTO"):
        """Get the average pressure.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           pressure : BioSimSpace.Types.Pressure
               The average pressure.
        """
        return self.getRecord("PRESSAVG", time_series, _Units.Pressure.bar, block)

    def getCurrentPressureAverage(self, time_series=False):
        """Get the current average pressure.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           pressure : BioSimSpace.Types.Pressure
               The average pressure.
        """
        return self.getPressureAverage(time_series, block=False)

    def getGPressure(self, time_series=False, block="AUTO"):
        """Get the pressure. (Hydrogens incorporated into bonded atoms.)

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           pressure : BioSimSpace.Types.Pressure
               The pressure.
        """
        return self.getRecord("GPRESSURE", time_series, _Units.Pressure.bar, block)

    def getCurrentGPressure(self, time_series=False):
        """Get the current pressure. (Hydrogens incorporated into bonded atoms.)

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           pressure : BioSimSpace.Types.Pressure
               The pressure.
        """
        return self.getGPressure(time_series, block=False)

    def getGPressureAverage(self, time_series=False, block="AUTO"):
        """Get the average pressure. (Hydrogens incorporated into bonded atoms.)

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           pressure : BioSimSpace.Types.Pressure
               The average pressure.
        """
        return self.getRecord("GPRESSAVG", time_series, _Units.Pressure.bar, block)

    def getCurrentGPressureAverage(self, time_series=False):
        """Get the current average pressure. (Hydrogens incorporated into bonded atoms.)

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           pressure : BioSimSpace.Types.Pressure
               The average pressure.
        """
        return self.getGPressureAverage(time_series, block=False)

    def getVolume(self, time_series=False, block="AUTO"):
        """Get the volume.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.

           block : bool
               Whether to block until the process has finished running.


           Returns
           -------

           volume : BioSimSpace.Types.Volume
               The volume.
        """
        return self.getRecord("VOLUME", time_series, _Units.Volume.angstrom3, block)

    def getCurrentVolume(self, time_series=False):
        """Get the current volume.

           Keyword arguments
           -----------------

           time_series : bool
               Whether to return a list of time series records.


           Returns
           -------

           volume : BioSimSpace.Types.Volume
               The volume.
        """
        return self.getVolume(time_series, block=False)

    def eta(self):
        """Get the estimated time for the process to finish (in minutes)."""

        # Make sure the list of stdout records is up to date.
        # Print the last zero lines, i.e. no output.
        self.stdout(0)

        # Now search backwards through the list to find the last TIMING record.
        for x, record in reversed(list(enumerate(self._stdout))):

            # Split the record using whitespace.
            data = record.split()

            # We've found a TIMING record.
            if len(data) > 0:
                if data[0] == "TIMING:":

                    # Try to find the "hours" record.
                    # If found, return the entry preceeding it.
                    try:
                        return (float(data[data.index("hours") - 1]) * 60) * _Units.Time.minutes

                    # No record found.
                    except ValueError:
                        return None

    def stdout(self, n=10):
        """Print the last n lines of the stdout buffer.

           Keyword arguments
           -----------------

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
                    if (len(stdout_data) == len(self._stdout_title)):
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

    def _get_stdout_record(self, key, time_series=False, unit=None):
        """Helper function to get a stdout record from the dictionary.

           Positional arguments
           --------------------

           key : str
               The record key.


           Keyword arguments
           -----------------

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
        if len(self._stdout_dict) is 0:
            return None

        if type(time_series) is not bool:
            _warnings.warn("Non-boolean time-series flag. Defaulting to False!")
            time_series = False

        # Valdate the unit.
        if unit is not None:
            if not isinstance(unit, _Type.Type):
                raise TypeError("'unit' must be of type 'BioSimSpace.Types'")

        # Return the list of dictionary values.
        if time_series:
            try:
                if key is "TS":
                    return [int(x) for x in self._stdout_dict[key]]
                else:
                    if unit is None:
                        return [float(x) for x in self._stdout_dict[key]]
                    else:
                        return [float(x) * unit for x in self._stdout_dict[key]]

            except KeyError:
                return None

        # Return the most recent dictionary value.
        else:
            try:
                if key is "TS":
                    return int(self._stdout_dict[key][-1])
                else:
                    if unit is None:
                        return float(self._stdout_dict[key][-1])
                    else:
                        return float(self._stdout_dict[key][-1]) * unit

            except KeyError:
                return None

    def _get_trajectory_files(self):
        """Get all files associated with the molecular trajectory."""

        # Name of the trajectory file.
        traj_file = "%s/%s_out.dcd" % (self._work_dir, self._name)

        # Return the trajectory and topology file.
        if _os.path.isfile(traj_file):
            return (traj_file, self._top_file)

        # No trajectory file.
        else:
            return None
