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

from ._atom_utils import _AToMUtils
import warnings as _warnings
import math as _math
from .._Exceptions import IncompatibleError as _IncompatibleError
from .. import Protocol as _Protocol
from ._openmm import OpenMM as _OpenMM


class OpenMMAToM(_OpenMM):
    """config generator functions for AToM simulations using OpenMM.
    Designed to overload the _generate_config() method of the OpenMM class
    to introduce AToM-specific methods."""

    def __init__(
        self,
        system=None,
        protocol=None,
        reference_system=None,
        exe=None,
        name="openmm",
        platform="CPU",
        work_dir=None,
        seed=None,
        property_map={},
        **kwargs,
    ):
        # Look for the is_testing flag in the kwargs.
        # Only used for calculating single point energies.
        if "_is_testing" in kwargs:
            _warnings.warn("NOW IN TESTING MODE")
            self._is_testing = kwargs["_is_testing"]
        else:
            self._is_testing = False
        super().__init__(
            system,
            protocol,
            reference_system=reference_system,
            exe=exe,
            name=name,
            platform=platform,
            work_dir=work_dir,
            seed=seed,
            property_map=property_map,
            **kwargs,
        )

    def _generate_config(self):
        if isinstance(self._protocol, _Protocol.AToMMinimisation):
            self._generate_config_minimisation()
        elif isinstance(self._protocol, _Protocol.AToMEquilibration):
            self._generate_config_equilibration()
        elif isinstance(self._protocol, _Protocol.AToMAnnealing):
            self._generate_config_annealing()
        elif isinstance(self._protocol, _Protocol.AToMProduction) and self._is_testing:
            self._generate_config_single_point_testing()
        elif isinstance(self._protocol, _Protocol.AToMProduction):
            self._generate_config_production()

    def _check_space(self):
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

        return has_box

    def _add_initialisation(self, has_box):
        # Write the OpenMM import statements.
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
        is_periodic = True
        if not has_box or not self._has_water:
            is_periodic = False
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

        return is_periodic

    def _add_pressure_check(self, pressure, temperature, is_periodic):
        # Add a Monte Carlo barostat if the simulation is at constant pressure.
        is_constant_pressure = False
        if pressure is not None:
            # Cannot use a barostat with a non-periodic system.
            if not is_periodic:
                _warnings.warn(
                    "Cannot use a barostat for a vacuum or non-periodic simulation"
                )
            else:
                is_constant_pressure = True

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

        return is_constant_pressure

    def _add_simulation_instantiation(self):
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

    def _generate_config_minimisation(self):
        util = _AToMUtils(self._protocol)
        # Clear the existing configuration list.
        self._config = []

        has_box = self._check_space()
        self._add_config_imports()
        self._add_config_monkey_patches()
        self._add_initialisation(has_box)

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
        if self._protocol._getCoreAlignment():
            alignment = util.createAlignmentForce()
            self.addToConfig("\n# Add alignment force.")
            self.addToConfig(alignment)
        if self._protocol._getCMCMRestraint():
            CMCM = util.createCOMRestraint()
            self.addToConfig("\n# Add COM restraint.")
            self.addToConfig(CMCM)

        self._add_simulation_instantiation()

        self.addToConfig(
            f"simulation.minimizeEnergy(maxIterations={self._protocol._getSteps()})"
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

    def _generate_config_equilibration(self):
        util = _AToMUtils(self._protocol)
        # Clear the existing configuration list.
        self._config = []

        has_box = self._check_space()
        self._add_config_imports()
        self._add_config_monkey_patches()
        is_periodic = self._add_initialisation(has_box)

        # Get the starting temperature and system pressure.
        temperature = self._protocol._getStartTemperature().kelvin().value()
        pressure = self._protocol._getPressure()

        is_constant_pressure = self._add_pressure_check(
            pressure, temperature, is_periodic
        )
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
        if self._protocol._getUseATMForce():
            atm = util.createATMForce(index=None)
            self.addToConfig(atm)

        if self._protocol._getCoreAlignment():
            alignment = util.createAlignmentForce()
            self.addToConfig("\n# Add alignment force.")
            self.addToConfig(alignment)
        if self._protocol._getCMCMRestraint():
            CMCM = util.createCOMRestraint()
            self.addToConfig("\n# Add COM restraint.")
            self.addToConfig(CMCM)

        # Get the integration time step from the protocol.
        timestep = self._protocol.getTimeStep().picoseconds().value()

        # Set the integrator.
        self.addToConfig("\n# Define the integrator.")
        self.addToConfig(f"integrator = LangevinMiddleIntegrator({temperature}*kelvin,")
        friction = 1 / self._protocol._getThermostatTimeConstant().picoseconds().value()
        self.addToConfig(f"                                {friction:.5f}/picosecond,")
        self.addToConfig(f"                                {timestep}*picoseconds)")
        if self._is_seeded:
            self.addToConfig(f"integrator.setRandomNumberSeed({self._seed})")

        # Add the platform information.
        self._add_config_platform()

        self._add_simulation_instantiation()

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
                    self._protocol._getEndTemperature().kelvin().value()
                    - self._protocol._getStartTemperature().kelvin().value()
                ) / temp_cycles

                self.addToConfig(f"start_temperature = {temperature}")
                self.addToConfig(f"for x in range(0, {temp_cycles}):")
                self.addToConfig(f"    temperature = {temperature} + x*{delta_temp}")
                self.addToConfig(f"    integrator.setTemperature(temperature*kelvin)")
                if is_constant_pressure:
                    self.addToConfig(
                        f"    barostat.setDefaultTemperature(temperature*kelvin)"
                    )
                self.addToConfig("    simulation.step(100)")
            else:
                # Work out the temperature change per step.
                delta_temp = (
                    self._protocol._getEndTemperature().kelvin().value()
                    - self._protocol._getStartTemperature().kelvin().value()
                ) / steps

                self.addToConfig(f"start_temperature = {temperature}")
                self.addToConfig(f"for x in range(0, {steps}):")
                self.addToConfig(f"    temperature = {temperature} + x*{delta_temp}")
                self.addToConfig(f"    integrator.setTemperature(temperature*kelvin)")
                if is_constant_pressure:
                    self.addToConfig(
                        f"    barostat.setDefaultTemperature(temperature*kelvin)"
                    )
                self.addToConfig("    simulation.step(1)")

    def _generate_config_annealing(self):
        self._protocol._set_current_index(0)
        util = _AToMUtils(self._protocol)
        # Clear the existing configuration list.
        self._config = []

        has_box = self._check_space()

        # Add standard openMM config
        self.addToConfig("from glob import glob")
        self.addToConfig("import math")
        self.addToConfig("import os")
        self.addToConfig("import shutil")
        self._add_config_imports()
        self._add_config_monkey_patches()

        is_periodic = self._add_initialisation(has_box)

        # Get the starting temperature and system pressure.
        temperature = self._protocol._getTemperature().kelvin().value()
        pressure = self._protocol._getPressure()

        is_constant_pressure = self._add_pressure_check(
            pressure, temperature, is_periodic
        )

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

        # Use utils to create AToM-specific forces
        # Atom force is the only window-dependent force
        disp = util.createDisplacement()
        self.addToConfig(disp)
        self.addToConfig("\n# Add AToM Force.")
        self.addToConfig(util.createATMForce(self._protocol._get_window_index()))
        if self._protocol._getCoreAlignment():
            alignment = util.createAlignmentForce()
            self.addToConfig("\n# Add alignment force.")
            self.addToConfig(alignment)

        if self._protocol._getCMCMRestraint():
            CMCM = util.createCOMRestraint()
            self.addToConfig("\n# Add COM restraint.")
            self.addToConfig(CMCM)

        # Get the integration time step from the protocol.
        timestep = self._protocol.getTimeStep().picoseconds().value()

        # Set the integrator.
        self.addToConfig("\n# Define the integrator.")
        self.addToConfig(f"integrator = LangevinMiddleIntegrator({temperature}*kelvin,")
        friction = 1 / self._protocol._getThermostatTimeConstant().picoseconds().value()
        self.addToConfig(f"                                {friction:.5f}/picosecond,")
        self.addToConfig(f"                                {timestep}*picoseconds)")
        if self._is_seeded:
            self.addToConfig(f"integrator.setRandomNumberSeed({self._seed})")

        # Add the platform information.
        self._add_config_platform()

        self._add_simulation_instantiation()

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

    def _generate_config_production(self):
        self._protocol._set_current_index(0)
        analysis_method = self._protocol._getAnalysisMethod()
        util = _AToMUtils(self._protocol)
        # Clear the existing configuration list.
        self._config = []

        has_box = self._check_space()

        # TODO: check extra_options, extra_lines and property_map
        if self._protocol._get_window_index() is None:
            raise _IncompatibleError(
                "AToM protocol requires the current window index to be set."
            )

        # Write the OpenMM import statements.

        self.addToConfig("import pandas as pd")
        self.addToConfig("import numpy as np")
        self.addToConfig("from glob import glob")
        self.addToConfig("import math")
        self.addToConfig("import os")
        self.addToConfig("import shutil")
        self._add_config_imports()
        self._add_config_monkey_patches()
        self.addToConfig("\n")
        if analysis_method == "UWHAM" or analysis_method == "both":
            self.addToConfig(util.createSoftcorePertE())
        # Add standard openMM config

        is_periodic = self._add_initialisation(has_box)
        # Get the starting temperature and system pressure.
        temperature = self._protocol._getTemperature().kelvin().value()
        pressure = self._protocol._getPressure()

        is_constant_pressure = self._add_pressure_check(
            pressure, temperature, is_periodic
        )

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

        # Use utils to create AToM-specific forces
        # Atom force is the only window-dependent force
        disp = util.createDisplacement()
        self.addToConfig(disp)
        self.addToConfig("\n# Add AToM Force.")
        self.addToConfig(util.createATMForce(self._protocol._get_window_index()))
        if self._protocol._getCoreAlignment():
            alignment = util.createAlignmentForce()
            self.addToConfig("\n# Add alignment force.")
            self.addToConfig(alignment)

        if self._protocol._getCMCMRestraint():
            CMCM = util.createCOMRestraint()
            self.addToConfig("\n# Add COM restraint.")
            self.addToConfig(CMCM)

        # Get the integration time step from the protocol.
        timestep = self._protocol.getTimeStep().picoseconds().value()

        # Set the integrator.
        self.addToConfig("\n# Define the integrator.")
        self.addToConfig(f"integrator = LangevinMiddleIntegrator({temperature}*kelvin,")
        friction = 1 / self._protocol._getThermostatTimeConstant().picoseconds().value()
        self.addToConfig(f"                                {friction:.5f}/picosecond,")
        self.addToConfig(f"                                {timestep}*picoseconds)")
        if self._is_seeded:
            self.addToConfig(f"integrator.setRandomNumberSeed({self._seed})")

        # Add the platform information.
        self._add_config_platform()

        self._add_simulation_instantiation()

        # Set initial velocities from temperature distribution.
        self.addToConfig("\n# Setting initial system velocities.")
        self.addToConfig(
            f"simulation.context.setVelocitiesToTemperature({temperature})"
        )

        # Check for a restart file and load the simulation state.
        is_restart, _ = self._add_config_restart()

        # NOTE: The restarting logic here is different to previous openMM classes
        # It doesn't use the steps value from the restart function, instead
        # the number of steps is worked out at runtime within the openmm script
        # this means that restarting either by using the biosimspace runner
        # OR by running the openmm script directly will work the same.
        step = 0

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

        # Get the report and restart intervals.
        report_interval = self._protocol.getReportInterval()
        restart_interval = self._protocol.getRestartInterval()

        # Cap the intervals at the total number of steps.
        if report_interval > steps:
            report_interval = steps
        if restart_interval > steps:
            restart_interval = steps

        # Work out the total simulation time in picoseconds.
        run_time = steps * timestep

        # Work out the number of cycles in 100 picosecond intervals.
        cycles = _math.ceil(run_time / (report_interval * timestep))

        # Work out the number of steps per cycle.
        steps_per_cycle = int(steps / cycles)

        self.addToConfig(
            util.createRestartLogic(
                total_cycles=cycles, steps_per_cycle=steps_per_cycle
            )
        )
        # Inform user that a restart was loaded.
        self.addToConfig("\n# Print restart information.")
        self.addToConfig("if is_restart:")
        self.addToConfig(f"    steps = {total_steps}")
        self.addToConfig("    percent_complete = 100 * (step / steps)")
        self.addToConfig("    print('Loaded state from an existing simulation.')")
        self.addToConfig("    print(f'Simulation is {percent_complete}% complete.')")
        self.addToConfig("    print(f'running an additional {numcycles} cycles')")

        # Add the reporters.
        self.addToConfig("\n# Add reporters.")
        self._add_config_reporters(
            state_interval=report_interval,
            traj_interval=restart_interval,
            is_restart=is_restart,
        )

        self.addToConfig(f"\ntemperature = {temperature}")
        if analysis_method == "UWHAM":
            # Now run the simulation.
            self.addToConfig(
                util.createSoftcorePertELoop(
                    name=self._name,
                    steps_per_cycle=steps_per_cycle,
                    report_interval=report_interval,
                    timestep=timestep,
                )
            )
        elif analysis_method == "both":
            direction = self._protocol._getDirection()
            inflex = 0
            for i in range(len(direction) - 1):
                if direction[i] != direction[i + 1]:
                    inflex = i + 1
                    break
            # Now run the simulation.
            self.addToConfig(
                util.createReportingBoth(
                    name=self._name,
                    steps_per_cycle=steps_per_cycle,
                    timestep=timestep,
                    inflex_point=inflex,
                )
            )
        else:
            direction = self._protocol._getDirection()
            inflex = 0
            for i in range(len(direction) - 1):
                if direction[i] != direction[i + 1]:
                    inflex = i + 1
                    break
            self.addToConfig(
                util.createLoopWithReporting(
                    name=self._name,
                    steps_per_cycle=steps_per_cycle,
                    report_interval=report_interval,
                    timestep=timestep,
                    inflex_point=inflex,
                )
            )

    def _generate_config_single_point_testing(self):
        # Designed as a hidden method - uses a production protocol to
        # calculate single point energies for each lambda window
        # quite hacky, but not designed to be exposed to the user anyway
        self._protocol._set_current_index(0)
        if not isinstance(self._protocol, _Protocol.AToMProduction):
            raise _IncompatibleError(
                "Single point testing requires an AToMProduction protocol."
            )
        util = _AToMUtils(self._protocol)
        # Clear the existing configuration list.
        self._config = []

        has_box = self._check_space()

        # TODO: check extra_options, extra_lines and property_map
        if self._protocol._get_window_index() is None:
            raise _IncompatibleError(
                "AToM protocol requires the current window index to be set."
            )

        # Write the OpenMM import statements.

        self.addToConfig("import pandas as pd")
        self.addToConfig("import numpy as np")
        self.addToConfig("from glob import glob")
        self.addToConfig("import math")
        self.addToConfig("import os")
        self.addToConfig("import shutil")
        self._add_config_imports()
        self._add_config_monkey_patches()
        self.addToConfig("\n")
        # Add standard openMM config

        is_periodic = self._add_initialisation(has_box)
        # Get the starting temperature and system pressure.
        temperature = self._protocol._getTemperature().kelvin().value()
        pressure = self._protocol._getPressure()

        is_constant_pressure = self._add_pressure_check(
            pressure, temperature, is_periodic
        )

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
            frc = util.create_flat_bottom_restraint(restrained_atoms, force_group=5)
            self.addToConfig(frc)

        # Use utils to create AToM-specific forces
        # Atom force is the only window-dependent force
        disp = util.createDisplacement()
        self.addToConfig(disp)
        self.addToConfig("\n# Add AToM Force.")
        self.addToConfig(
            util.createATMForce(self._protocol._get_window_index(), force_group=10)
        )
        if self._protocol._getCoreAlignment():
            alignment = util.createAlignmentForce(force_group=[6, 7, 8])
            self.addToConfig("\n# Add alignment force.")
            self.addToConfig(alignment)

        if self._protocol._getCMCMRestraint():
            CMCM = util.createCOMRestraint(force_group=9)
            self.addToConfig("\n# Add COM restraint.")
            self.addToConfig(CMCM)

        # Get the integration time step from the protocol.
        timestep = self._protocol.getTimeStep().picoseconds().value()

        # Set the integrator.
        self.addToConfig("\n# Define the integrator.")
        self.addToConfig(f"integrator = LangevinMiddleIntegrator({temperature}*kelvin,")
        friction = 1 / self._protocol._getThermostatTimeConstant().picoseconds().value()
        self.addToConfig(f"                                {friction:.5f}/picosecond,")
        self.addToConfig(f"                                {timestep}*picoseconds)")
        if self._is_seeded:
            self.addToConfig(f"integrator.setRandomNumberSeed({self._seed})")

        # Add the platform information.
        self._add_config_platform()

        self._add_simulation_instantiation()

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
        # reading in the directions from the protocol, find the index at which direction changes
        direction = self._protocol._getDirection()
        inflex = 0
        for i in range(len(direction) - 1):
            if direction[i] != direction[i + 1]:
                inflex = i + 1
                break
        self.addToConfig(
            util.createSinglePointTest(
                inflex,
                self._name,
                atm_force_group=10,
                position_restraint_force_group=5,
                alignment_force_groups=[6, 7, 8],
                com_force_group=9,
            )
        )
