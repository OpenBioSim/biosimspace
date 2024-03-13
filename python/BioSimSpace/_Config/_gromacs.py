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

"""Functionality for generating configuration files for AMBER."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Gromacs"]

import math as _math
import warnings as _warnings

from .. import Protocol as _Protocol
from ..Protocol._free_energy_mixin import _FreeEnergyMixin
from ..Protocol._position_restraint_mixin import _PositionRestraintMixin

from ._config import Config as _Config


class Gromacs(_Config):
    """A class for generating configuration files for GROMACS."""

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

    def createConfig(self, version=None, extra_options={}, extra_lines=[]):
        """
        Create the list of configuration strings.

        version : float
            The GROMACS version.

        extra_options : dict
            A dictionary containing extra options. Overrides the defaults generated
            by the protocol.

        extra_lines : [str]
            A list of extra lines to put at the end of the configuration file.

        Returns
        -------

        config : [str]
            The list of AMBER format configuration strings.
        """

        # Validate input.

        if version and not isinstance(version, float):
            raise TypeError("'version' must be of type 'float'.")

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

        # Make sure the report interval is a multiple of nstcalcenergy.
        if isinstance(self._protocol, _FreeEnergyMixin):
            nstcalcenergy = 250
        else:
            nstcalcenergy = 100
        report_interval = self.reportInterval()
        if report_interval % nstcalcenergy != 0:
            report_interval = nstcalcenergy * _math.ceil(
                report_interval / nstcalcenergy
            )

        # Define some miscellaneous defaults.
        protocol_dict = {
            # Interval between writing to the log file.
            "nstlog": report_interval,
            # Interval between writing to the energy file.
            "nstenergy": report_interval,
            # Interval between writing to the trajectory file.
            "nstxout-compressed": self.restartInterval(),
        }

        # Minimisation.
        if isinstance(self._protocol, _Protocol.Minimisation):
            protocol_dict["integrator"] = "steep"
        else:
            # Timestep in picoseconds
            timestep = self._protocol.getTimeStep().picoseconds().value()
            # Integration time step.
            protocol_dict["dt"] = f"{timestep:.3f}"
        # Number of integration steps.
        protocol_dict["nsteps"] = self.steps()

        # Constraints.
        if not isinstance(self._protocol, _Protocol.Minimisation):
            # Rigid bonded hydrogens.
            protocol_dict["constraints"] = "h-bonds"
            # Linear constraint solver.
            protocol_dict["constraint-algorithm"] = "LINCS"

        # Periodic boundary conditions.

        # Simulate a fully periodic box.
        protocol_dict["pbc"] = "xyz"
        # Use Verlet pair lists.
        protocol_dict["cutoff-scheme"] = "Verlet"
        if self.hasBox(self._system, self._property_map) and self.hasWater(
            self._system
        ):
            # Use a grid to search for neighbours.
            protocol_dict["ns-type"] = "grid"
            # Rebuild neighbour list every 20 steps.
            protocol_dict["nstlist"] = "20"
            # Set short-range cutoff.
            protocol_dict["rlist"] = "0.8"
            # Set van der Waals cutoff.
            protocol_dict["rvdw"] = "0.8"
            # Set Coulomb cutoff.
            protocol_dict["rcoulomb"] = "0.8"
            # Fast smooth Particle-Mesh Ewald.
            protocol_dict["coulombtype"] = "PME"
            # Dispersion corrections for energy and pressure.
            protocol_dict["DispCorr"] = "EnerPres"
        else:
            # Perform vacuum simulations by implementing pseudo-PBC conditions,
            # i.e. run calculation in a near-infinite box (333.3 nm).
            # c.f.: https://pubmed.ncbi.nlm.nih.gov/29678588
            # Single neighbour list (all particles interact).
            protocol_dict["nstlist"] = "1"
            # "Infinite" short-range cutoff.
            protocol_dict["rlist"] = "333.3"
            # "Infinite" van der Waals cutoff.
            protocol_dict["rvdw"] = "333.3"
            # "Infinite" Coulomb cutoff.
            protocol_dict["rcoulomb"] = "333.3"
            # Plain cut-off.
            protocol_dict["coulombtype"] = "Cut-off"
        # Twin-range van der Waals cut-off.
        protocol_dict["vdwtype"] = "Cut-off"

        # Position restraints.
        if isinstance(self._protocol, _PositionRestraintMixin):
            # Note that constraints will be defined by the GROMACS process.
            protocol_dict["refcoord-scaling"] = "com"

        # Pressure control.
        if not isinstance(self._protocol, _Protocol.Minimisation):
            if self._protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if self.hasBox(self._system, self._property_map) and self.hasWater(
                    self._system
                ):
                    # Barostat type.
                    if version and version >= 2021:
                        protocol_dict["pcoupl"] = "c-rescale"
                    else:
                        protocol_dict["pcoupl"] = "berendsen"
                    # 1ps time constant for pressure coupling.
                    protocol_dict["tau-p"] = 1
                    # Pressure in bar.
                    protocol_dict[
                        "ref-p"
                    ] = f"{self._protocol.getPressure().bar().value():.5f}"
                    # Compressibility of water.
                    protocol_dict["compressibility"] = "4.5e-5"
                else:
                    _warnings.warn(
                        "Cannot use a barostat for a vacuum or non-periodic simulation"
                    )

        # Temperature control.
        if not isinstance(self._protocol, _Protocol.Minimisation):
            # Leap-frog molecular dynamics.
            protocol_dict["integrator"] = "md"
            # Temperature coupling using velocity rescaling with a stochastic term.
            protocol_dict["tcoupl"] = "v-rescale"
            # A single temperature group for the entire system.
            protocol_dict["tc-grps"] = "system"
            # Thermostat coupling frequency (ps).
            protocol_dict["tau-t"] = "{:.5f}".format(
                self._protocol.getThermostatTimeConstant().picoseconds().value()
            )

            if isinstance(self._protocol, _Protocol.Equilibration):
                if self._protocol.isConstantTemp():
                    temp = (
                        "%.2f" % self._protocol.getStartTemperature().kelvin().value()
                    )
                    protocol_dict["ref-t"] = temp
                    protocol_dict["gen-vel"] = "yes"
                    protocol_dict["gen-temp"] = temp
                else:
                    # still need a reference temperature for each group,
                    # even when heating/cooling
                    protocol_dict["ref-t"] = (
                        "%.2f" % self._protocol.getEndTemperature().kelvin().value()
                    )
                    # Work out the final time of the simulation.
                    timestep = self._protocol.getTimeStep().picoseconds().value()
                    end_time = _math.floor(timestep * self.steps())
                    # Single sequence of annealing points.
                    protocol_dict["annealing"] = "single"
                    # Two annealing points for "system" temperature group.
                    protocol_dict["annealing-npoints"] = 2

                    # Linearly change temperature between start and end times.
                    protocol_dict["annealing-time"] = "0 %d" % end_time
                    protocol_dict["annealing-temp"] = "%.2f %.2f" % (
                        self._protocol.getStartTemperature().kelvin().value(),
                        self._protocol.getEndTemperature().kelvin().value(),
                    )
            else:
                # Fixed temperature.
                protocol_dict["ref-t"] = (
                    "%.2f" % self._protocol.getTemperature().kelvin().value()
                )

        # Free energies.
        if isinstance(self._protocol, _FreeEnergyMixin):
            # Extract the lambda array.
            lam_vals = self._protocol.getLambdaValues()

            # Determine the index of the lambda value.
            idx = self._protocol.getLambdaIndex()

            # Free energy mode.
            protocol_dict["free-energy"] = "yes"
            # Initial lambda state.
            protocol_dict["init-lambda-state"] = idx
            # Lambda value array.
            protocol_dict["fep-lambdas"] = " ".join([str(x) for x in lam_vals])
            # All interactions on at lambda = 0
            protocol_dict["couple-lambda0"] = "vdw-q"
            # All interactions on at lambda = 1
            protocol_dict["couple-lambda1"] = "vdw-q"
            # Write all lambda values.
            protocol_dict["calc-lambda-neighbors"] = -1
            # Calculate energies every 250 steps.
            protocol_dict["nstcalcenergy"] = 250
            # Write gradients every 250 steps.
            protocol_dict["nstdhdl"] = 250

        # Put everything together in a line-by-line format.
        total_dict = {**protocol_dict, **extra_options}
        total_lines = [
            f"{k} = {v}" for k, v in total_dict.items() if v is not None
        ] + extra_lines

        return total_lines
