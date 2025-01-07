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

"""Functionality for production protocols."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Production"]

import math as _math
import warnings as _warnings

from ._position_restraint import _PositionRestraintMixin
from ._protocol import Protocol as _Protocol
from .. import Types as _Types
from .. import Units as _Units


class Production(_Protocol, _PositionRestraintMixin):
    """A class for storing production protocols."""

    def __init__(
        self,
        timestep=_Types.Time(2, "femtosecond"),
        runtime=_Types.Time(1, "nanosecond"),
        temperature=_Types.Temperature(300, "kelvin"),
        pressure=_Types.Pressure(1, "atmosphere"),
        tau_t=_Types.Time(1, "picosecond"),
        report_interval=200,
        restart_interval=1000,
        first_step=0,
        restart=False,
        restraint=None,
        force_constant=10 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
    ):
        """Constructor.

        Parameters
        ----------

        timestep : :class:`Time <BioSimSpace.Types.Time>`
            The integration timestep.

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The running time.

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The temperature.

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure. Pass pressure=None to use the NVT ensemble.

        tau_t : :class:`Time <BioSimSpace.Types.Time>`
            Time constant for thermostat coupling.

        report_interval : int
            The frequency at which statistics are recorded. (In integration steps.)

        restart_interval : int
            The frequency at which restart configurations and trajectory

        first_step : int
            The initial time step (for restart simulations).

        restart : bool
            Whether this is a continuation of a previous simulation.

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
        _Protocol.__init__(self)

        # Set the time step.
        self.setTimeStep(timestep)

        # Set the runtime.
        self.setRunTime(runtime)

        # Set the system temperature.
        self.setTemperature(temperature)

        # Set the system pressure.
        if pressure is not None:
            self.setPressure(pressure)
        else:
            self._pressure = None

        # Set the time constant for the thermostat
        self.setTauT(tau_t)

        # Set the report interval.
        self.setReportInterval(report_interval)

        # Set the restart interval.
        self.setRestartInterval(restart_interval)

        # Set the restart flag.
        self.setRestart(restart)

        # Set the first time step.
        self.setFirstStep(first_step)

        # Set the restraint.
        _PositionRestraintMixin.__init__(self, restraint, force_constant)

    def _get_parm(self):
        """Return a string representation of the parameters."""
        return (
            f"timestep={self._timestep}, "
            f"runtime={self._runtime}, "
            f"temperature={self._temperature}, "
            f"pressure={self._pressure}, "
            f"tau_t={self._tau_t}, "
            f"report_interval={self._report_interval}, "
            f"restart_interval={self._restart_interval}, "
            f"first_step={self._first_step}, "
            f"restart={self._restart}, " + _PositionRestraintMixin._get_parm(self)
        )

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return f"<BioSimSpace.Protocol.Production: {self._get_parm()}>"

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._is_customised:
            return "BioSimSpace.Protocol.Custom"
        else:
            return f"BioSimSpace.Protocol.Production({self._get_parm()})"

    def __eq__(self, other):
        """Equality operator."""

        if not isinstance(other, Production):
            return False

        if self._is_customised or other._is_customised:
            return False

        return (
            self._timestep == other._timestep
            and self._runtime == other._runtime
            and self._temperature == other._temperature
            and self._pressure == other._pressure
            and self._tau_t == other._tau_t
            and self._report_interval == other._report_interval
            and self._restart_interval == other._restart_interval
            and self._restart == other._restart
            and self._first_step == other._first_step
            and _PositionRestraintMixin.__eq__(self, other)
        )

    def getTimeStep(self):
        """
        Return the time step.

        Returns
        -------

        timestep : :class:`Time <BioSimSpace.Types.Time>`
            The integration time step.
        """
        return self._timestep

    def setTimeStep(self, timestep):
        """
        Set the time step.

        Parameters
        ----------

        timestep : str, :class:`Time <BioSimSpace.Types.Time>`
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

    def getTauT(self):
        """
        Return the time constant for the thermostat.

        Returns
        -------

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The time constant for the thermostat.
        """
        return self._tau_t

    def setTauT(self, tau_t):
        """
        Set the time constant for the thermostat.

        Parameters
        ----------

        tau_t : str, :class:`Time <BioSimSpace.Types.Time>`
            The time constant for the thermostat.
        """
        if isinstance(tau_t, str):
            try:
                self._tau_t = _Types.Time(tau_t)
            except:
                raise ValueError("Unable to parse 'tau_t' string.") from None
        elif isinstance(tau_t, _Types.Time):
            self._tau_t = tau_t
        else:
            raise TypeError("'tau_t' must be of type 'str' or 'BioSimSpace.Types.Time'")

    def getTemperature(self):
        """
        Return temperature.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The simulation temperature.
        """
        return self._temperature

    def setTemperature(self, temperature):
        """
        Set the temperature.

        Parameters
        ----------

        temperature : str, :class:`Temperature <BioSimSpace.Types.Temperature>`
            The simulation temperature.
        """
        if isinstance(temperature, str):
            try:
                self._temperature = _Types.Temperature(temperature)
            except:
                raise ValueError("Unable to parse 'temperature' string.") from None
        elif isinstance(temperature, _Types.Temperature):
            self._temperature = temperature
        else:
            raise TypeError(
                "'temperature' must be of type 'str' or 'BioSimSpace.Types.Temperature'"
            )

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
            _warnings.warn("'report_interval' must be positive. Using default (200).")
            report_interval = 200

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
            _warnings.warn("'restart_interval' must be positive. Using default (1000).")
            restart_interval = 1000

        self._restart_interval = restart_interval

    def getFirstStep(self):
        """
        Return the first time step.

        Returns
        -------

        step : int
            The first time step.
        """
        return self._first_step

    def setFirstStep(self, first_step):
        """
        Set the initial time step.

        Parameters
        ----------

        step : int
            The first time step.
        """
        if not type(first_step) is int:
            raise TypeError("'first_step' must be of type 'int'")

        if first_step < 0:
            _warnings.warn("The initial time step must be positive. Using default (0).")
            self._first_step = 0
        else:
            self._first_step = _math.ceil(first_step)

    def isRestart(self):
        """
        Return whether this restart simulation.

        Returns
        -------

        is_restart : bool
            Whether this is a restart simulation.
        """
        return self._restart

    def setRestart(self, restart):
        """
        Set the restart flag.

        Parameters
        ----------

        restart : bool
            Whether this is a restart simulation.
        """
        if isinstance(restart, bool):
            self._restart = restart
        else:
            _warnings.warn("Non-boolean restart flag. Defaulting to False!")
            self._restart = False
