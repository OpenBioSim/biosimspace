#####################################################################
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

"""Functionality for metadynamics protocols."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Metadynamics"]


from .. import Types as _Types
from ..Metadynamics import CollectiveVariable as _CollectiveVariable

from ._protocol import Protocol as _Protocol

# Store the collective variable base type.
_colvar_type = _CollectiveVariable._collective_variable.CollectiveVariable


class Metadynamics(_Protocol):
    """A class for storing metadynamics protocols."""

    def __init__(
        self,
        collective_variable,
        timestep=_Types.Time(2, "femtosecond"),
        runtime=_Types.Time(1, "nanosecond"),
        temperature=_Types.Temperature(300, "kelvin"),
        pressure=_Types.Pressure(1, "atmosphere"),
        thermostat_time_constant=_Types.Time(1, "picosecond"),
        hill_height=_Types.Energy(1, "kj per mol"),
        hill_frequency=1000,
        report_interval=1000,
        restart_interval=1000,
        bias_factor=None,
    ):
        """
        Constructor.

        Parameters
        ----------

        collective_variable : :class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`, \
                             [:class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`]
            The collective variable (or variables) for the simulation.

        timestep : :class:`Time <BioSimSpace.Types.Time>`
            The integration timestep.

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The running time.

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The temperature.

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure. Pass pressure=None to use the NVT ensemble.

        thermostat_time_constant : :class:`Time <BioSimSpace.Types.Time>`
            Time constant for thermostat coupling.

        hill_height : :class:`Energy <BioSimSpace.Types.Energy>`
            The height of the Gaussian hills.

        hill_frequency : int
            The frequency at which hills are deposited.

        report_interval : int
            The frequency at which statistics are recorded. (In integration steps.)

        restart_interval : int
            The frequency at which restart configurations and trajectory
            frames are saved. (In integration steps.)

        bias_factor : float
            The bias factor for well tempered metadynamics.
        """

        # Call the base class constructor.
        super().__init__()

        # Whether this is a newly created object.
        self._is_new_object = True

        # Set the collective variable.
        self.setCollectiveVariable(collective_variable)

        # Set the time step.
        self.setTimeStep(timestep)

        # Set the runtime.
        self.setRunTime(runtime)

        # Set the system temperature.
        self.setTemperature(temperature)

        # Set the report interval.
        self.setReportInterval(report_interval)

        # Set the restart interval.
        self.setRestartInterval(restart_interval)

        # Set the system pressure.
        if pressure is not None:
            self.setPressure(pressure)
        else:
            self._pressure = None

        # Set the thermostat time constant.
        self.setThermostatTimeConstant(thermostat_time_constant)

        # Set the hill parameters: height, frequency.
        self.setHillHeight(hill_height)
        self.setHillFrequency(hill_frequency)

        # Set the bias factor for well tempered metadynamics.
        if bias_factor is not None:
            self.setBiasFactor(bias_factor)
        else:
            self._bias_factor = None

        # Flag that the object has been created.
        self._is_new_object = False

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            string = "<BioSimSpace.Protocol.Metadynamics: "
            string += "collective_variable=%s" % self._collective_variable
            string += ", timestep=%s" % self._timestep
            string += ", runtime=%s" % self._runtime
            string += ", temperature=%s" % self._temperature
            if self._pressure is not None:
                string += ", pressure=%s" % self._pressure
            string += ", hill_height=%s" % self._hill_height
            string += ", hill_frequency=%s" % self._hill_frequency
            if self._bias_factor is not None:
                string += ", bias_factor=%s" % self._bias_factor
            string += ", report_interval=%d" % self._report_interval
            string += ", restart_interval=%d" % self._restart_interval
            string += ">"

            return string

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return self.__str__()

    def __eq__(self, other):
        """Equality operator."""

        if not isinstance(other, Metadynamics):
            return False

        if self._is_customised or other._is_customised:
            return False

        return (
            self._collective_variable == other._collective_variable
            and self._timestep == other._timestep
            and self._runtime == other._runtime
            and self._temperature == other._temperature
            and self._pressure == other._pressure
            and self._hill_height == other._hill_height
            and self._hill_frequency == other._hill_frequency
            and self._bias_factor == other._bias_factor
            and self._report_interval == other._report_interval
            and self._restart_interval == other._restart_interval
        )

    def getCollectiveVariable(self):
        """
        Return the collective variable (or variables).

        Returns
        -------

        collective_variable : [:class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`]
            The collective variable (or variables) for the simulation.
        """
        return self._collective_variable.copy()

    def setCollectiveVariable(self, collective_variable):
        """
        Set the collective variable (or variables).

        Parameters
        ----------

        collective_variable : :class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`, \
                             [:class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`]
            The collective variable (or variables) for the simulation.
        """

        # A single collective variable.
        if isinstance(collective_variable, _colvar_type):
            self._collective_variable = [collective_variable]
            return

        if isinstance(collective_variable, (list, tuple)):
            if not all(isinstance(x, _colvar_type) for x in collective_variable):
                raise TypeError(
                    "'collective_variable' must all be of type "
                    "'BioSimSpace.Metadynamics.CollectiveVariable'"
                )
        else:
            raise TypeError(
                "'collective_variable' must be of type "
                "'BioSimSpace.Metadynamics.CollectiveVariable' "
                "or a list of 'BioSimSpace.Metadynamics.CollectiveVariable' types."
            )

        # Make sure all of the collective variables are consistent. If any have
        # a grid set, then so must all other variables.

        num_grid = 0

        for colvar in collective_variable:
            if colvar.getGrid() is not None:
                num_grid += 1

        if num_grid > 0 and num_grid != len(collective_variable):
            raise ValueError(
                "If a 'grid' is desired, then all collective "
                "variables must define one."
            )

        self._collective_variable = collective_variable

        # If the object has already been created, then check that other member
        # data is consistent.
        if not self._is_new_object:
            self.setHillHeight(self._hill_height)

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

    def getHillHeight(self):
        """
        Return the height of the Gaussian hills.

        Returns
        -------

        hill_width : :class:`Energy <BioSimSpace.Types.Energy>`
            The height of the Gaussian hills.
        """
        return self._hill_height

    def setHillHeight(self, hill_height):
        """
        Set the height of the Gaussian hills.

        Parameters
        ----------

        hill_height : str, :class:`Energy <BioSimSpace.Types.Energy>`
            The hill height.
        """

        if isinstance(hill_height, str):
            try:
                hill_height = _Types.Energy(hill_height)
            except:
                raise ValueError("Unable to parse 'hill_height' string.") from None
        elif not isinstance(hill_height, _Types.Energy):
            raise TypeError(
                "'hill_height' must be of type 'str' or 'BioSimSpace.Types.Energy'"
            )

        # Check that heights is greater than zero.
        if hill_height.value() <= 0:
            raise ValueError("Cannot have 'hill_height' with value <= 0")

        self._hill_height = hill_height.kj_per_mol()

    def getHillFrequency(self):
        """
        Return the frequency at which Gaussian hills are deposited.

        Returns
        -------

        hill_frequency : int
            The frequency at which hills are deposited.
        """
        return self._hill_frequency

    def setHillFrequency(self, hill_frequency):
        """
        Set the frequency at which Gaussian hills are deposited.

        Parameters
        ----------

        hill_frequency : int
            The frequency at which hills are deposited.
        """

        try:
            hill_frequency = int(hill_frequency)
        except:
            raise TypeError("'hill_frequency' must be of type 'int'")

        if hill_frequency < 1:
            raise ValueError("'hill_frequency' must be >= 1")

        self._hill_frequency = hill_frequency

    def getBiasFactor(self):
        """
        Return the bias factor for well tempered metadynamics.

        Returns
        -------

        bias_factor : float
            The bias factor for well tempered metadynamics.
        """
        return self._bias_factor

    def setBiasFactor(self, bias_factor=None):
        """
        Set the bias factor for well tempered metadynamics.
        Call with no arguments to clear the bias factor.

        Parameters
        ----------

        bias_factor : float
            The bias factor for well tempered metadynamics.
        """

        if bias_factor is None:
            self._bias_factor = None
            return

        try:
            bias_factor = float(bias_factor)
        except:
            raise TypeError("'bias_factor' must be of type 'float'")

        if bias_factor < 1.0:
            raise ValueError("'bias_factor' must be > 1.0")

        # OpenMM crashes when the bias factor is exactly one, so add a delta.
        delta = 1e-6
        if abs(bias_factor - 1.0) < delta:
            bias_factor += delta

        self._bias_factor = bias_factor

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
        import warnings as _warnings

        if not type(restart_interval) is int:
            raise TypeError("'restart_interval' must be of type 'int'")

        if restart_interval <= 0:
            _warnings.warn("'restart_interval' must be positive. Using default (500).")
            restart_interval = 500

        self._restart_interval = restart_interval
