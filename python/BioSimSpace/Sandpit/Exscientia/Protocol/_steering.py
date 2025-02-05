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

"""Functionality for steered molecular dynamics protocols."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Steering"]

import math as _math
import os as _os

from .. import Types as _Types
from ..Metadynamics import CollectiveVariable as _CollectiveVariable
from ..Metadynamics import Restraint as _Restraint

from ._protocol import Protocol as _Protocol

# Store the collective variable base type.
_colvar_type = _CollectiveVariable._collective_variable.CollectiveVariable


class Steering(_Protocol):
    """A class for storing steered molecular dynamics protocols."""

    def __init__(
        self,
        collective_variable,
        schedule,
        restraints,
        verse="both",
        timestep=_Types.Time(2, "femtosecond"),
        runtime=_Types.Time(1, "nanosecond"),
        temperature=_Types.Temperature(300, "kelvin"),
        pressure=_Types.Pressure(1, "atmosphere"),
        report_interval=1000,
        restart_interval=1000,
        colvar_file=None,
    ):
        """
        Constructor.

        Parameters
        ----------

        collective_variable : :class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`, \
                             [:class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`]
            The collective variable (or variables) used to perform the steering.

        schedule : [:class:`Time <BioSimSpace.Types.Time>`]
            The time schedule for the steering.

        restraints : [:class:`Restraint <BioSimSpace.Metadynamics.Restraint>`] \
                    [(:class:`Restraint <BioSimSpace.Metadynamics.Restraint>`)]]
            The position of the restraint on each collective variable for
            each stage of the schedule.

        verse : str
            Whether the restraint is acting for values of the collective
            variable "larger" or "smaller" than the restraint, or acting
            on "both" sides (default).

        timestep : :class:`Time <BioSimSpace.Types.Time>`
            The integration timestep.

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The running time.

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The temperature.

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure. Pass pressure=None to use the NVT ensemble.

        restart_interval : int
            The frequency at which restart configurations and trajectory
            frames are saved. (In integration steps.)

        colvar_file : str
            The path to a COLVAR file from a previous simulation. The
            information in the file must be consistent with the
            'collective_variable' argument.
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

        # Set the verse.
        self.setVerse(verse)

        # Set the schedule.
        self.setSchedule(schedule)

        # Set the restraints.
        self.setRestraints(restraints)

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

        if colvar_file is not None:
            self.setColvarFile(colvar_file)
        else:
            self._colvar_file = None

        # Flag that the object has been created.
        self._is_new_object = False

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            string = "<BioSimSpace.Protocol.Steering: "
            string += "collective_variable=%s" % self._collective_variable
            string += ", schedule=%s" % self._schedule
            string += ", restraints=%s" % self._restraints
            string += ", verse=%s" % self._verse
            string += ", timestep=%s" % self._timestep
            string += ", runtime=%s" % self._runtime
            string += ", temperature=%s" % self._temperature
            if self._pressure is not None:
                string += ", pressure=%s" % self._pressure
            if self._colvar_file is not None:
                string += ", colvar_file=%r" % self._colvar_file
            string += ", report_interval=%d" % self._report_interval
            string += ", restart_interval=%d" % self._restart_interval
            string += ">"

            return string

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return self.__str__()

    def __eq__(self, other):
        """Equality operator."""

        if not isinstance(other, Steering):
            return False

        if self._is_customised or other._is_customised:
            return False

        return (
            self._collective_variable == other._collective_variable
            and self._schedule == other._schedule
            and self._restraints == other._restraints
            and self._verse == other._verse
            and self._timestep == other._timestep
            and self._runtime == other._runtime
            and self._temperature == other._temperature
            and self._pressure == other._pressure
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

        else:
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

            self._collective_variable = collective_variable

        # If the object has already been created, then check that other member
        # data is consistent.
        if not self._is_new_object:
            self.setRestraints(self._restraints)

    def getSchedule(self):
        """
        Return steering schedule.

        Returns
        -------

        schedule : [int]
            The schedule for the steering, i.e. the integration time steps
            at which restraints are applied/adjusted.
        """
        steps = []
        for time in self._schedule:
            steps.append(_math.floor(time / self._timestep))
        return steps

    def setSchedule(self, schedule):
        """
        Set the steering schedule.

        Parameters
        ----------

        schedule : [:class:`Time <BioSimSpace.Types.Tme>`]
            The time schedule for the steering.
        """

        if isinstance(schedule, (list, tuple)):
            if not all(isinstance(x, _Types.Time) for x in schedule):
                raise TypeError(
                    "'schedule' must all be of type " "'BioSimSpace.Types.Time'"
                )
        else:
            raise TypeError(
                "'schedule' must be a list of " "'BioSimSpace.Types.Time' types."
            )

        # Make sure the times are linearly increasing and are less than
        # the total run time.
        last_time = schedule[0]
        for x in range(1, len(schedule)):
            if last_time > self._runtime:
                raise ValueError("'schedule' values cannot exceed the 'runtime'!")
            if schedule[x] > self._runtime:
                raise ValueError("'schedule' values cannot exceed the 'runtime'!")
            if schedule[x] < last_time:
                raise ValueError("The steering 'schedule' must increase in time!")
            last_time = schedule[x]

        self._schedule = schedule

        # If the object has already been created, then check that other member
        # data is consistent.
        if not self._is_new_object:
            self.setRestraints(self._restraints)

    def getRestraints(self):
        """
        Return the restraint on each collective variable for each stage in
        the schedule.

        Returns
        -------

        restraints : [:class`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`] \
                     [(:class:`CollectiveVariable <BioSimSpace.Metadynamics.CollectiveVariable>`)]]
            The position of the restraint on each collective variable for
            each stage of the schedule.
        """
        return self._restraints.copy()

    def setRestraints(self, restraints):
        """
        Set the restraints for the steering schedule.

        Parameters
        ----------

        restraint : [:class`Restraint <BioSimSpace.Metadynamics.Restraint>`] \
                    [(:class:`Restraint <BioSimSpace.Metadynamics.Restraint>`)]]
            The position of the restraint on each collective variable for
            each stage of the schedule.
        """

        # Validate type.
        if not isinstance(restraints, (list, tuple)):
            raise TypeError(
                "'restraints' must be a list of "
                "'BioSimSpace.CollectiveVariable.Restraint' types."
            )

        # Make sure we are self-consistent with the schedule.
        if len(restraints) != len(self._schedule):
            raise ValueError(
                f"'len(restraints) != len(schedule), i.e. {len(restraints)} != {len(self._schedule)}"
            )

        # Convert all single entries to lists and validate types.
        new_restraints = []
        num_cvs = len(self._collective_variable)
        for restraint in restraints:
            # Convert tuple to list.
            if isinstance(restraint, tuple):
                restraint = list(restraint)
            # Convert any non-list type to a list.
            if not isinstance(restraint, list):
                restraint = [restraint]
            # Validate that there is a restraint for each collective variable
            # for each stage in the schedule.
            if len(restraint) != num_cvs:
                raise ValueError(
                    "Must have a restraint for each collective "
                    "variable for each stage of the schedule."
                )
            # Validate type.
            if not all(isinstance(x, _Restraint) for x in restraint):
                raise TypeError(
                    "'restraint' must all be of type "
                    "'BioSimSpace.Metadynamics.Restraint'"
                )
            # Validate that restraints are harmonic.
            if not all(r.getSlope() == 0 for r in restraint):
                raise ValueError("'restraints' can only contain harmonic restraints!")
            # Validate that the value of the restraint matches the type of each
            # collective variable for each stage in the schedule.
            if not all(
                type(r.getValue()) in self._collective_variable[x]._types
                for x, r in enumerate(restraint)
            ):
                raise ValueError(
                    "The type of value for each restraint must match the "
                    "type of the collective variable to which it corresponds."
                )

            # Append to the new list.
            new_restraints.append(restraint)

        self._restraints = new_restraints

    def getVerse(self):
        """
        Returns whether the restraint is acting for values of the collective
        variable "larger" or "smaller" than the restraint, or acting on "both"
        sides (default).

        Returns
        -------

        verse : str
            Whether the restraint is acting for values of the collective
            variable "larger" or "smaller" than the restraint, or acting
            on "both" sides (default).
        """
        return self._verse

    def setVerse(self, verse):
        """
        Set whether the restraint is acting for values of the collective
        variable "larger" or "smaller" than the restraint, or acting on "both"
        sides (default).

        Parameters
        ----------

        verse : str, [str]
            Whether the restraint is acting for values of the collective
            variable "larger" or "smaller" than the restraint, or acting
            on "both" sides (default).
        """

        # Convert string to list.
        if isinstance(verse, str):
            # Use the same verse for all collective variables.
            verse = [verse] * len(self.getCollectiveVariable())

        if isinstance(verse, (list, tuple)):
            if not all(isinstance(x, str) for x in verse):
                raise TypeError(
                    "'verse' must be of type 'str' or a list of 'str' types."
                )
            elif len(verse) != len(self.getCollectiveVariable()):
                raise ValueError(
                    "'verse' must be of same length as 'collective_variable'"
                )
        else:
            raise TypeError("'verse' must be of type 'str' or a list of 'str' types.")

        new_verse = []
        allowed = ["both", "larger", "smaller"]

        for v in verse:
            # Convert to lower case and strip whitespace.
            v = v.lower().replace(" ", "")

            if v not in allowed:
                raise ValueError(f"'verse' must be one of: {allowed}")
            else:
                new_verse.append(v)

        self._verse = new_verse

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

        # If the object has already been created, then check that other member
        # data is consistent.
        if not self._is_new_object:
            self.setSchedule(self._schedule)

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

        # If the object has already been created, then check that other member
        # data is consistent.
        if not self._is_new_object:
            self.setSchedule(self._schedule)

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

    def getColvarFile(self):
        """
        Return the path to the COLVAR file.

        Returns
        -------

        colvar_file : str
            The path to the COLVAR file.
        """
        return self._colvar_file

    def setColvarFile(self, colvar_file):
        """
        Set the location of an existing COLVAR file.

        Parameters
        ----------

        colvar_file : str
            The path to an existing COLVAR file.
        """
        if not isinstance(colvar_file, str):
            raise ValueError("'colvar_file' must be of type 'str'")

        if not _os.path.isfile(colvar_file):
            raise ValueError("'colvar_file' doesn't exist: %s" % colvar_file)

        self._colvar_file = colvar_file
