######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2025
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

"""Functionality for generating configuration files for SOMD."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Somd"]

import math as _math
import warnings as _warnings

from .. import Protocol as _Protocol
from ..Protocol._free_energy_mixin import _FreeEnergyMixin
from ..Protocol._position_restraint_mixin import _PositionRestraintMixin

from ._config import Config as _Config


class Somd(_Config):
    """A class for generating configuration files for SOMD."""

    def __init__(self, system, protocol, property_map={}):
        """
        Constructor.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.

        protocol : :class:`Protocol <BioSimSpace.Protocol>`
            The protocol for the process.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(system, protocol, property_map=property_map)

    def createConfig(self, extra_options={}, extra_lines=[]):
        """
        Create the list of configuration strings.

        extra_options : dict
            A dictionary containing extra options. Overrides the defaults generated
            by the protocol.

        extra_lines : [str]
            A list of extra lines to put at the end of the configuration file.

        Returns
        -------

        config : [str]
            The list of SOMD format configuration strings.
        """

        # Validate input.

        if not isinstance(extra_options, dict):
            raise TypeError("'extra_options' must be of type 'dict'.")
        else:
            keys = extra_options.keys()
            if not all(isinstance(k, str) for k in keys):
                raise TypeError("Keys of 'extra_options' must be of type 'str'.")

        if not isinstance(extra_lines, list):
            raise TypeError("'extra_lines' must be of type 'list'.")
        else:
            if not all(isinstance(line, str) for line in extra_lines):
                raise TypeError("Lines in 'extra_lines' must be of type 'str'.")

        # Define some miscellaneous defaults.
        protocol_dict = {"save coordinates": True}  # Save molecular coordinates.

        # Minimisation.
        if isinstance(self._protocol, _Protocol.Minimisation):
            protocol_dict["minimise"] = True
            # Maximum number of steps.
            protocol_dict["minimise maximum iterations"] = self.steps()
            # Convergence tolerance.
            protocol_dict["minimise tolerance"] = 1
            # Perform a single SOMD cycle.
            protocol_dict["ncycles"] = 1
            # Perform a single MD move.
            protocol_dict["nmoves"] = 1
        else:
            # Get the report and restart intervals.
            report_interval = self._protocol.getReportInterval()
            restart_interval = self._protocol.getRestartInterval()

            # Work out the number of cycles. We want one per nanosecond of simulation.
            if self._protocol.getRunTime().nanoseconds().value() < 2:
                ncycles = 1
                moves_per_cycle = self.steps()
            else:
                ncycles = _math.ceil(self._protocol.getRunTime().nanoseconds().value())
                moves_per_cycle = _math.ceil(self.steps() / ncycles)

            # The number of moves needs to be multiple of the report interval.
            if moves_per_cycle % report_interval != 0:
                moves_per_cycle = (
                    _math.ceil(moves_per_cycle / report_interval) * report_interval
                )

            # Work out the number of cycles per frame.
            cycles_per_frame = restart_interval / moves_per_cycle

            # Work out whether we need to adjust the buffer frequency.
            buffer_freq = 0
            if cycles_per_frame < 1:
                buffer_freq = restart_interval
                cycles_per_frame = 1
            else:
                cycles_per_frame = _math.floor(cycles_per_frame)

            # Make sure that we aren't buffering more than 1000 frames per cycle.
            if buffer_freq > 0 and moves_per_cycle / buffer_freq > 1000:
                _warnings.warn(
                    "Trajectory buffering will exceed limit. Reducing buffer frequency."
                )
                buffer_freq = moves_per_cycle / 1000

            # For free energy simulations, the buffer frequency must be an integer
            # multiple of the frequency at which free energies are written, which
            # is report interval. Round down to the closest multiple.
            if isinstance(self._protocol, _FreeEnergyMixin):
                if buffer_freq > 0:
                    buffer_freq = report_interval * _math.floor(
                        buffer_freq / report_interval
                    )

            # The number of SOMD cycles.
            protocol_dict["ncycles"] = ncycles
            # The number of moves per cycle.
            protocol_dict["nmoves"] = moves_per_cycle
            # Cycles per trajectory write.
            protocol_dict["ncycles_per_snap"] = cycles_per_frame
            # Buffering frequency.
            protocol_dict["buffered coordinates frequency"] = buffer_freq
            timestep = self._protocol.getTimeStep().femtoseconds().value()
            # Integration time step.
            protocol_dict["timestep"] = "%.2f femtosecond" % timestep

            # Use the Langevin Middle integrator if it is a 4 fs timestep
            if timestep >= 4.00:
                # Langevin middle integrator
                protocol_dict["integrator_type"] = "langevinmiddle"
            else:
                pass

        # Periodic boundary conditions.
        if self.hasWater(self._system):
            # Solvated box.
            protocol_dict["reaction field dielectric"] = "78.3"
        if not self.hasBox(self._system, self._property_map) or not self.hasWater(
            self._system
        ):
            # No periodic box.
            protocol_dict["cutoff type"] = "cutoffnonperiodic"
        else:
            # Periodic box.
            protocol_dict["cutoff type"] = "cutoffperiodic"
        # Non-bonded cut-off.
        protocol_dict["cutoff distance"] = "10 angstrom"

        # Restraints.
        if (
            isinstance(self._protocol, _PositionRestraintMixin)
            and self._protocol.getRestraint() is not None
        ):
            raise _IncompatibleError(
                "We currently don't support position restraints with SOMD."
            )

        # Pressure control.
        protocol_dict["barostat"] = False
        if not isinstance(self._protocol, _Protocol.Minimisation):
            if self._protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if self.hasBox(self._system, self._property_map) and self.hasWater(
                    self._system
                ):
                    # Enable barostat.
                    protocol_dict["barostat"] = True
                    pressure = self._protocol.getPressure().atm().value()
                    # Presure in atmosphere.
                    protocol_dict["pressure"] = "%.5f atm" % pressure
                else:
                    _warnings.warn(
                        "Cannot use a barostat for a vacuum or non-periodic simulation"
                    )

        # Temperature control.
        if not isinstance(self._protocol, _Protocol.Minimisation):
            if (
                isinstance(self._protocol, _Protocol.Equilibration)
                and not self._protocol.isConstantTemp()
            ):
                raise _IncompatibleError(
                    "SOMD only supports constant temperature equilibration."
                )

            # Turn on the thermostat.
            protocol_dict["thermostat"] = "True"
            if not isinstance(self._protocol, _Protocol.Equilibration):
                protocol_dict["temperature"] = (
                    "%.2f kelvin" % self._protocol.getTemperature().kelvin().value()
                )
            else:
                protocol_dict["temperature"] = (
                    "%.2f kelvin"
                    % self._protocol.getStartTemperature().kelvin().value()
                )

            # Friction coefficient (1 / ps). Note that the unit is written
            # to the configuration file as "picosecond". This is because SOMD
            # applies the value to OpenMM without inverting it, so has things the
            # wrong way round, i.e. it uses the inverse friction as the friction.
            protocol_dict["inverse friction"] = "{:.5f} picosecond".format(
                1 / self._protocol.getThermostatTimeConstant().picoseconds().value()
            )

        # Free energies.
        if isinstance(self._protocol, _Protocol.FreeEnergyProduction):
            # Perform a pre-minimisation.
            protocol_dict["minimise"] = True
            # Handle hydrogen perturbations.
            protocol_dict["constraint"] = "hbonds-notperturbed"
            # Write gradients at the report interval.
            protocol_dict["energy frequency"] = report_interval

            protocol = [str(x) for x in self._protocol.getLambdaValues()]
            protocol_dict["lambda array"] = ", ".join(protocol)
            # Current lambda value.
            protocol_dict["lambda_val"] = self._protocol.getLambda()
            res_num = (
                self._system.search("perturbable")
                .residues()[0]
                ._sire_object.number()
                .value()
            )
            # Perturbed residue number.
            protocol_dict["perturbed residue number"] = res_num

        # Put everything together in a line-by-line format.
        total_dict = {**protocol_dict, **extra_options}
        total_lines = [
            f"{k} = {v}" for k, v in total_dict.items() if v is not None
        ] + extra_lines

        return total_lines
