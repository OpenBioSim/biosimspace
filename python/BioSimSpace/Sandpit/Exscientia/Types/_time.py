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

"""A time type."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Time"]

from sire.legacy import Units as _SireUnits

from ._type import Type as _Type


class Time(_Type):
    """A time type."""

    # A list of the supported Sire unit names.
    _sire_units = [
        "day",
        "hour",
        "minute",
        "second",
        "millisecond",
        "nanosecond",
        "picosecond",
        "femtosecond",
    ]

    # Dictionary of allowed units.
    _supported_units = {
        "DAY": _SireUnits.day,
        "HOUR": _SireUnits.hour,
        "MINUTE": _SireUnits.minute,
        "SECOND": _SireUnits.second,
        "MILLISECOND": _SireUnits.millisecond,
        "NANOSECOND": _SireUnits.nanosecond,
        "PICOSECOND": _SireUnits.picosecond,
        "FEMTOSECOND": _SireUnits.femtosecond,
    }

    # Map unit abbreviations to the full name.
    _abbreviations = {
        "HR": "HOUR",
        "MIN": "MINUTE",
        "SEC": "SECOND",
        "MS": "MILLISECOND",
        "NS": "NANOSECOND",
        "PS": "PICOSECOND",
        "FS": "FEMTOSECOND",
    }

    # Print formatting.
    _print_format = {
        "DAY": "day",
        "HOUR": "hour",
        "MINUTE": "min",
        "SECOND": "sec",
        "MILLISECOND": "ms",
        "NANOSECOND": "ns",
        "PICOSECOND": "ps",
        "FEMTOSECOND": "fs",
    }

    # Documentation strings.
    _doc_strings = {
        "DAY": "A time in days.",
        "HOUR": "A time in hours.",
        "MINUTE": "A time in minutes.",
        "SECOND": "A time in seconds.",
        "MILLISECOND": "A time in milliseconds.",
        "NANOSECOND": "A time in nanoseconds.",
        "PICOSECOND": "A time in picoseconds.",
        "FEMTOSECOND": "A time in femtoseconds.",
    }

    # Null type unit for avoiding issue printing configargparse help.
    _default_unit = "NANOSECOND"

    # The dimension mask:
    #     Angle, Charge, Length, Mass, Quantity, Temperature, Time
    _dimensions = (0, 0, 0, 0, 0, 0, 1)

    def __init__(self, *args):
        """
        Constructor.

        ``*args`` can be a value and unit, or a string representation
        of the time, e.g. "0.2 fs".

        Parameters
        ----------

        value : float
            The value.

        unit : str
            The unit.

        string : str
            A string representation of the time.

        Examples
        --------

        Create an object representing a time of 17.3 femtoseconds then
        print the time in nanoseconds.

        >>> import BioSimSpace as BSS
        >>> time = BSS.Types.Time(17.3, "fs")
        >>> print(time.nanoseconds())

        The same as above, except passing a string representation of the
        time to the constructor.

        >>> import BioSimSpace as BSS
        >>> time = BSS.Types.Time("17.3 fs")
        >>> print(time.nanoseconds())

        The string matching is extremeley flexible, so all of the following
        would be valid arguments: "17.3 fs", "17.3 femtoseconds", "1.73e1 fs".
        """

        # Call the base class constructor.
        super().__init__(*args)

    def __str__(self):
        """Return a human readable string representation of the object."""

        abbrev = self._print_format[self._unit]
        if self._value != 1:
            if abbrev[-1] != "s":
                abbrev = abbrev + "s"

        if abs(self._value) > 1e4 or abs(self._value) < 1e-4:
            return "%.4e %s" % (self._value, abbrev)
        else:
            return "%5.4f %s" % (self._value, abbrev)

    def days(self):
        """
        Return the time in days.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The time in days.
        """
        return Time(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.day), "DAY"
        )

    def hours(self):
        """
        Return the time in hours.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The time in hours.
        """
        return Time(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.hour),
            "HOUR",
        )

    def minutes(self):
        """
        Return the time in minutes.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The time in minutes.
        """
        return Time(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.minute),
            "MINUTE",
        )

    def seconds(self):
        """
        Return the time in seconds.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The time in seconds.
        """
        return Time(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.second),
            "SECOND",
        )

    def milliseconds(self):
        """
        Return the time in milliseconds.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The time in milliseconds.
        """
        return Time(
            (self._value * self._supported_units[self._unit]).to(
                _SireUnits.millisecond
            ),
            "MILLISECOND",
        )

    def nanoseconds(self):
        """
        Return the time in nanoseconds.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The time in nanoseconds.
        """
        return Time(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.nanosecond),
            "NANOSECOND",
        )

    def picoseconds(self):
        """
        Return the time in picoseconds.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The time in picoseconds.
        """
        return Time(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.picosecond),
            "PICOSECOND",
        )

    def femtoseconds(self):
        """
        Return the time in femtoseconds.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The time in femtoseconds.
        """
        return Time(
            (self._value * self._supported_units[self._unit]).to(
                _SireUnits.femtosecond
            ),
            "FEMTOSECOND",
        )

    def _to_default_unit(self, mag=None):
        """
        Internal method to return an object of the same type in the default unit.

        Parameters
        ----------

        mag : float
            The value (optional).

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The time in the default unit of picoseconds.
        """
        if mag is None:
            return self.picoseconds()
        else:
            return Time(mag, "PICOSECOND")

    def _convert_to(self, unit):
        """
        Return the time in a different unit.

        Parameters
        ----------

        unit : str
            The unit to convert to.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The time in the specified unit.
        """
        if unit == "DAY":
            return self.days()
        elif unit == "HOUR":
            return self.hours()
        elif unit == "MINUTE":
            return self.minutes()
        elif unit == "SECOND":
            return self.seconds()
        elif unit == "MILLISECOND":
            return self.milliseconds()
        elif unit == "NANOSECOND":
            return self.nanoseconds()
        elif unit == "PICOSECOND":
            return self.picoseconds()
        elif unit == "FEMTOSECOND":
            return self.femtoseconds()
        else:
            raise ValueError(
                "Supported units are: '%s'" % list(self._supported_units.keys())
            )

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Check that the unit is supported.
        if unit in self._supported_units:
            return unit
        elif unit[:-1] in self._supported_units:
            return unit[:-1]
        elif unit in self._abbreviations:
            return self._abbreviations[unit]
        elif unit[:-1] in self._abbreviations:
            return self._abbreviations[unit[:-1]]
        else:
            raise ValueError(
                "Supported units are: '%s'" % list(self._supported_units.keys())
            )

    @staticmethod
    def _to_sire_format(unit):
        """
        Reformat the unit string so it adheres to the Sire unit formatting.

        Parameters
        ----------

        unit : str
            A string representation of the unit.

        Returns
        -------

        sire_unit : str
            The unit string in Sire compatible format.
        """

        unit = unit.replace("days", "day")
        unit = unit.replace("seconds", "second")
        unit = unit.replace("secs", "second")

        # Convert powers. (Just 2nd and third for now.)
        unit = unit.replace("day2", "(day*day)")
        unit = unit.replace("day3", "(day*day*day)")
        unit = unit.replace("day-1", "(1/day)")
        unit = unit.replace("day-2", "(1/(day*day))")
        unit = unit.replace("day-3", "(1/(day*day*day))")
        unit = unit.replace("hour2", "(hour*hour)")
        unit = unit.replace("hour3", "(hour*hour*hour)")
        unit = unit.replace("hour-1", "(1/hour)")
        unit = unit.replace("hour-2", "(1/(hour*hour))")
        unit = unit.replace("hour-3", "(1/(hour*hour*hour))")
        unit = unit.replace("minute2", "(minute*minute)")
        unit = unit.replace("minute3", "(minute*minute*minute)")
        unit = unit.replace("minute-1", "(1/minute)")
        unit = unit.replace("minute-2", "(1/(minute*minute))")
        unit = unit.replace("minute-3", "(1/(minute*minute*minute))")
        unit = unit.replace("femtosecond2", "(femtosecond*femtosecond)")
        unit = unit.replace("femtosecond3", "(femtosecond*femtosecond*femtosecond)")
        unit = unit.replace("femtosecond-1", "(1/femtosecond)")
        unit = unit.replace("femtosecond-2", "(1/(femtosecond*femtosecond))")
        unit = unit.replace(
            "femtosecond-3", "(1/(femtosecond*femtosecond*femtosecond))"
        )
        unit = unit.replace("picosecond2", "(picosecond*picosecond)")
        unit = unit.replace("picosecond3", "(picosecond*picosecond*picosecond)")
        unit = unit.replace("picosecond-1", "(1/picosecond)")
        unit = unit.replace("picosecond-2", "(1/(picosecond*picosecond))")
        unit = unit.replace("picosecond-3", "(1/(picosecond*picosecond*picosecond))")
        unit = unit.replace("nanosecond2", "(nanosecond*nanosecond)")
        unit = unit.replace("nanosecond3", "(nanosecond*nanosecond*nanosecond)")
        unit = unit.replace("nanosecond-1", "(1/nanosecond)")
        unit = unit.replace("nanosecond-2", "(1/(nanosecond*nanosecond))")
        unit = unit.replace("nanosecond-3", "(1/(nanosecond*nanosecond*nanosecond))")
        unit = unit.replace("millisecond2", "(millisecond*millisecond)")
        unit = unit.replace("millisecond3", "(millisecond*millisecond*millisecond)")
        unit = unit.replace("millisecond-1", "(1/millisecond)")
        unit = unit.replace("millisecond-2", "(1/(millisecond*millisecond))")
        unit = unit.replace(
            "millisecond-3", "(1/(millisecond*millisecond*millisecond))"
        )
        unit = unit.replace("second2", "(second*second)")
        unit = unit.replace("second3", "(second*second*second)")
        unit = unit.replace("second-1", "(1/second)")
        unit = unit.replace("second-2", "(1/(second*second))")
        unit = unit.replace("second-3", "(1/(second*second*second))")

        return unit
