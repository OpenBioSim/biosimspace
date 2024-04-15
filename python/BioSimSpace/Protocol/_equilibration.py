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

"""Functionality for equilibration protocols."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Equilibration"]

import math as _math
import warnings as _warnings

from .. import Types as _Types
from .. import Units as _Units

from ._position_restraint_mixin import _PositionRestraintMixin
from ._protocol import Protocol as _Protocol


class Equilibration(_Protocol, _PositionRestraintMixin):
    """A class for storing equilibration protocols."""

    def __init__(
        self,
        timestep=_Types.Time(2, "femtosecond"),
        runtime=_Types.Time(0.2, "nanoseconds"),
        temperature_start=_Types.Temperature(300, "kelvin"),
        temperature_end=_Types.Temperature(300, "kelvin"),
        temperature=None,
        pressure=None,
        thermostat_time_constant=_Types.Time(1, "picosecond"),
        report_interval=100,
        restart_interval=500,
        restraint=None,
        force_constant=10 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
    ):
        """
        Constructor.

        Parameters
        ----------

        timestep : :class:`Time <BioSimSpace.Types.Time>`
            The integration timestep.

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The running time.

        temperature_start : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The starting temperature.

        temperature_end : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The final temperature.

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The equilibration temperature. This takes precedence of over
            the other temperatures, i.e. to run at fixed temperature.

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure. If this argument is omitted then the simulation
            is run using the NVT ensemble.

        thermostat_time_constant : :class:`Time <BioSimSpace.Types.Time>`
            Time constant for thermostat coupling.

        report_interval : int
            The frequency at which statistics are recorded. (In integration steps.)

        restart_interval : int
            The frequency at which restart configurations and trajectory
            frames are saved. (In integration steps.)

        restraint : str, [int]
            The type of restraint to perform. This should be one of the
            following options:
                "backbone"
                     Protein backbone atoms. The matching is done by a name
                     template, so is unreliable on conversion between
                     molecular file formats.
                "heavy"
                     All non-hydrogen atoms that aren't part of water
                     molecules or free ions.
                "all"
                     All atoms that aren't part of water molecules or free
                     ions.
            Alternatively, the user can pass a list of atom indices for
            more fine-grained control. If None, then no restraints are used.

        force_constant : :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`, float
            The force constant for the restraint potential. If a 'float' is
            passed, then default units of 'kcal_per_mol / angstrom**2' will
            be used.
        """

        # Call the base class constructor.
        super().__init__()

        # Set the time step.
        self.setTimeStep(timestep)

        # Set the running time.
        self.setRunTime(runtime)

        # Constant temperature equilibration.
        if temperature is not None:
            self.setStartTemperature(temperature)
            self.setEndTemperature(temperature)
            self._is_const_temp = True

        # Heating / cooling simulation.
        else:
            self._is_const_temp = False

            # Set the start temperature.
            self.setStartTemperature(temperature_start)

            # Set the final temperature.
            self.setEndTemperature(temperature_end)

            # Constant temperature simulation.
            if self._temperature_start == self._temperature_end:
                self._is_const_temp = True

        # Constant pressure simulation.
        if pressure is not None:
            self.setPressure(pressure)
        else:
            self._pressure = None

        # Set the thermostat time constant.
        self.setThermostatTimeConstant(thermostat_time_constant)

        # Set the report interval.
        self.setReportInterval(report_interval)

        # Set the restart interval.
        self.setRestartInterval(restart_interval)

        # Set the posistion restraint.
        _PositionRestraintMixin.__init__(self, restraint, force_constant)

    def _get_parm(self):
        """Return a string representation of the parameters."""
        return (
            f"timestep={self._timestep}, "
            f"runtime={self._runtime}, "
            f"temperature_start={self._temperature_start}, "
            f"temperature_end={self._temperature_end}, "
            f"pressure={self._pressure}, "
            f"thermostat_time_constant={self._thermostat_time_constant}, "
            f"report_interval={self._report_interval}, "
            f"restart_interval={self._restart_interval}, "
            + _PositionRestraintMixin._get_parm(self)
        )

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return f"<BioSimSpace.Protocol.Equilibration: {self._get_parm()}>"

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._is_customised:
            return "BioSimSpace.Protocol.Custom"
        else:
            return f"BioSimSpace.Protocol.Equilibration({self._get_parm()})"

    def __eq__(self, other):
        """Equality operator."""

        if not isinstance(other, Equilibration):
            return False

        if self._is_customised or other._is_customised:
            return False

        return (
            self._timestep == other._timestep
            and self._runtime == other._runtime
            and self._temperature_start == other._temperature_start
            and self._temperature_end == other._temperature_end
            and self._pressure == other._pressure
            and self._thermostat_time_constant == other._thermostat_time_constant
            and self._report_interval == other._report_interval
            and self._restart_interval == other._restart_interval
            and _PositionRestraintMixin.__eq__(self, other)
        )

    def getTimeStep(self):
        """
        Return the time step.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The integration time step.
        """
        return self._timestep

    def setTimeStep(self, timestep):
        """
        Set the time step.

        Parameters
        ----------

        time : str, :class:`Time <BioSimSpace.Types.Time>`
            The integration time step.
        """
        if isinstance(timestep, str):
            try:
                self._timestep = _Types.Time(timestep)
            except:
                raise ValueError("Unable to parse 'timestep' string.") from None
        elif isinstance(timestep, _Types.Time):
            self._timestep = timestep
        else:
            raise TypeError(
                "'timestep' must be of type 'str' or 'BioSimSpace.Types.Time'"
            )

    def getRunTime(self):
        """
        Return the running time.

        Returns
        -------

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The simulation run time.
        """
        return self._runtime

    def setRunTime(self, runtime):
        """
        Set the running time.

        Parameters
        ----------

        runtime : str, :class:`Time <BioSimSpace.Types.Time>`
            The simulation run time.
        """
        if isinstance(runtime, str):
            try:
                self._runtime = _Types.Time(runtime)
            except:
                raise ValueError("Unable to parse 'runtime' string.") from None
        elif isinstance(runtime, _Types.Time):
            self._runtime = runtime
        else:
            raise TypeError(
                "'runtime' must be of type 'str' or 'BioSimSpace.Types.Time'"
            )

    def getStartTemperature(self):
        """
        Return the starting temperature.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The starting temperature.
        """
        return self._temperature_start

    def setStartTemperature(self, temperature):
        """
        Set the starting temperature.

        Parameters
        ----------

        temperature : str, :class:`Temperature <BioSimSpace.Types.Temperature>`
            The starting temperature.
        """

        if isinstance(temperature, str):
            try:
                temperature = _Types.Temperature(temperature)
            except:
                raise ValueError("Unable to parse 'temperature' string.") from None
        elif not isinstance(temperature, _Types.Temperature):
            raise TypeError(
                "'temperature' must be of type 'str' or 'BioSimSpace.Types.Temperature'"
            )

        if _math.isclose(temperature.kelvin().value(), 0, rel_tol=1e-6):
            temperature._value = 0.01
        self._temperature_start = temperature

    def getEndTemperature(self):
        """
        Return the final temperature.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The final temperature.
        """
        return self._temperature_end

    def setEndTemperature(self, temperature):
        """
        Set the final temperature.

        Parameters
        ----------

        temperature : str, :class:`Temperature <BioSimSpace.Types.Temperature>`
            The final temperature.
        """
        if isinstance(temperature, str):
            try:
                temperature = _Types.Temperature(temperature)
            except:
                raise ValueError("Unable to parse 'temperature' string.") from None
        elif not isinstance(temperature, _Types.Temperature):
            raise TypeError(
                "'temperature' must be of type 'str' or 'BioSimSpace.Types.Temperature'"
            )

        if _math.isclose(temperature.kelvin().value(), 0, rel_tol=1e-6):
            temperature._value = 0.01
        self._temperature_end = temperature

    def getPressure(self):
        """
        Return the pressure.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure.
        """
        return self._pressure

    def setPressure(self, pressure):
        """
        Set the pressure.

        Parameters
        ----------

        pressure : str, :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure.
        """
        if isinstance(pressure, str):
            try:
                self._pressure = _Types.Pressure(pressure)
            except:
                raise ValueError("Unable to parse 'pressure' string.") from None
        elif isinstance(pressure, _Types.Pressure):
            self._pressure = pressure
        else:
            raise TypeError(
                "'pressure' must be of type 'str' or 'BioSimSpace.Types.Pressure'"
            )

    def getThermostatTimeConstant(self):
        """
        Return the time constant for the thermostat.

        Returns
        -------

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The time constant for the thermostat.
        """
        return self._thermostat_time_constant

    def setThermostatTimeConstant(self, thermostat_time_constant):
        """
        Set the time constant for the thermostat.

        Parameters
        ----------

        thermostat_time_constant : str, :class:`Time <BioSimSpace.Types.Time>`
            The time constant for the thermostat.
        """
        if isinstance(thermostat_time_constant, str):
            try:
                self._thermostat_time_constant = _Types.Time(thermostat_time_constant)
            except:
                raise ValueError(
                    "Unable to parse 'thermostat_time_constant' string."
                ) from None
        elif isinstance(thermostat_time_constant, _Types.Time):
            self._thermostat_time_constant = thermostat_time_constant
        else:
            raise TypeError(
                "'thermostat_time_constant' must be of type 'str' or 'BioSimSpace.Types.Time'"
            )

    def getReportInterval(self):
        """
        Return the interval between reporting statistics. (In integration steps.).

        Returns
        -------

        report_interval : int
            The number of integration steps between reporting statistics.
        """
        return self._report_interval

    def setReportInterval(self, report_interval):
        """
        Set the interval at which statistics are reported. (In integration steps.).

        Parameters
        ----------

        report_interval : int
            The number of integration steps between reporting statistics.
        """
        if not type(report_interval) is int:
            raise TypeError("'report_interval' must be of type 'int'")

        if report_interval <= 0:
            _warnings.warn("'report_interval' must be positive. Using default (100).")
            report_interval = 100

        self._report_interval = report_interval

    def getRestartInterval(self):
        """
        Return the interval between saving restart confiugrations, and/or
        trajectory frames. (In integration steps.).

        Returns
        -------

        restart_interval : int
            The number of integration steps between saving restart
            configurations and/or trajectory frames.
        """
        return self._restart_interval

    def setRestartInterval(self, restart_interval):
        """
        Set the interval between saving restart confiugrations, and/or
        trajectory frames. (In integration steps.).

        Parameters
        ----------

        restart_interval : int
            The number of integration steps between saving restart
            configurations and/or trajectory frames.
        """
        if not type(restart_interval) is int:
            raise TypeError("'restart_interval' must be of type 'int'")

        if restart_interval <= 0:
            _warnings.warn("'restart_interval' must be positive. Using default (500).")
            restart_interval = 500

        self._restart_interval = restart_interval

    def isConstantTemp(self):
        """
        Return whether the protocol has a constant temperature.

        Returns
        -------

        is_const_temp : bool
            Whether the temperature is fixed.
        """
        return self._temperature_start == self._temperature_end

    @classmethod
    def restraints(cls):
        """
        Return a list of the supported restraint keywords.

        Returns
        -------

        restraints : [str]
            A list of the supported restraint keywords.
        """
        return cls._restraints.copy()
