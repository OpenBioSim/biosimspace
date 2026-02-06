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

"""Functionality for equilibration protocols."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Equilibration"]


from .. import Types as _Types
from .. import Units as _Units

from ._protocol import Protocol as _Protocol
from ._position_restraint import _PositionRestraintMixin


class Equilibration(_Protocol, _PositionRestraintMixin):
    """A class for storing equilibration protocols."""

    # Supported restraint keywords.
    _restraints = ["backbone", "heavy", "all", "none"]

    def __init__(
        self,
        timestep=_Types.Time(2, "femtosecond"),
        runtime=_Types.Time(0.2, "nanoseconds"),
        temperature_start=_Types.Temperature(300, "kelvin"),
        temperature_end=_Types.Temperature(300, "kelvin"),
        temperature=None,
        pressure=None,
        tau_t=_Types.Time(1, "picosecond"),
        report_interval=200,
        restart_interval=1000,
        restraint=None,
        force_constant=10 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
        restart=False,
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
            The pressure. If None, then the simulation is run using the NVT ensemble.
            If set to a value, then the simulation is run using the NPT ensemble.

        tau_t : :class:`Time <BioSimSpace.Types.Time>`
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

        restart : bool
            Whether this is a continuation of a previous simulation.
        """

        # Call the base class constructor.

        _Protocol.__init__(self)

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

        # Set the time constant for the thermostat
        self.setTauT(tau_t)

        # Set the report interval.
        self.setReportInterval(report_interval)

        # Set the restart interval.
        self.setRestartInterval(restart_interval)

        # Set the restart flag.
        self.setRestart(restart)

        # Set the restraint.
        _PositionRestraintMixin.__init__(self, restraint, force_constant)

    def _get_parm(self):
        """Return a string representation of the parameters."""
        return (
            f"timestep={self._timestep}, "
            f"runtime={self._runtime}, "
            f"temperature_start={self._temperature_start}, "
            f"temperature_end={self._temperature_end}, "
            f"pressure={self._pressure}, "
            f"tau_t={self._tau_t}, "
            f"report_interval={self._report_interval}, "
            f"restart_interval={self._restart_interval}, "
            f"restart={self._restart}, " + _PositionRestraintMixin._get_parm(self)
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
            and self._tau_t == other._tau_t
            and self._report_interval == other._report_interval
            and self._restart_interval == other._restart_interval
            and self._restart == other._restart
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
        import math as _math

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
        import math as _math

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
        import warnings as _warnings

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
        import warnings as _warnings

        if not type(restart_interval) is int:
            raise TypeError("'restart_interval' must be of type 'int'")

        if restart_interval <= 0:
            _warnings.warn("'restart_interval' must be positive. Using default (1000).")
            restart_interval = 1000

        self._restart_interval = restart_interval

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
        import warnings as _warnings

        if isinstance(restart, bool):
            self._restart = restart
        else:
            _warnings.warn("Non-boolean restart flag. Defaulting to False!")
            self._restart = False

    def getRestraint(self):
        """
        Return the type of restraint..

        Returns
        -------

        restraint : str, [int]
            The type of restraint.
        """
        return self._restraint

    def setRestraint(self, restraint):
        """
        Set the type of restraint.

        Parameters
        ----------

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
            more fine-grained control.
        """

        if isinstance(restraint, str):
            # Convert to lower case and strip whitespace.
            restraint = restraint.lower().replace(" ", "")
            if restraint not in self._restraints:
                raise ValueError(f"'restraint' must be one of: {self._restraints}")
            # Set to NoneType if equal to "none", since this makes checking
            # whether a restraint is set elsewhere much easier.
            if restraint == "none":
                restraint = None

        elif isinstance(restraint, (list, tuple)):
            if not all(type(x) is int for x in restraint):
                raise ValueError("'restraint' must be a list of 'int' types!")
            # Create a set to sort and ensure no duplicates, then convert back to a list.
            restraint = list(set(restraint))
            restraint.sort()

        else:
            raise TypeError(
                "'restraint' must be of type 'str', or a list of 'int' types."
            )

        self._restraint = restraint

    def isConstantTemp(self):
        """
        Return whether the protocol has a constant temperature.

        Returns
        -------

        is_const_temp : bool
            Whether the temperature is fixed.
        """
        return self._temperature_start == self._temperature_end
