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

"""An angle type."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Angle"]

from sire.legacy import Units as _SireUnits

from ._type import Type as _Type


class Angle(_Type):
    """An angle type."""

    # A list of the supported Sire unit names.
    _sire_units = ["radian", "degree"]

    # Dictionary of allowed units.
    _supported_units = {"RADIAN": _SireUnits.radian, "DEGREE": _SireUnits.degree}

    # Map unit abbreviations to the full name.
    _abbreviations = {"R": "RADIAN", "D": "DEGREE"}

    # Print format.
    _print_format = {"RADIAN": "radian", "DEGREE": "degree"}

    # Documentation strings.
    _doc_strings = {"RADIAN": "An angle in radians.", "DEGREE": "An angle in degrees."}

    # Null type unit for avoiding issue printing configargparse help.
    _default_unit = "RADIAN"

    # The dimension mask:
    #     Angle, Charge, Length, Mass, Quantity, Temperature, Time
    _dimensions = (1, 0, 0, 0, 0, 0, 0)

    def __init__(self, *args):
        """
        Constructor.

        ``*args`` can be a value and unit, or a string representation
        of the angle, e.g. "3 radians".

        Parameters
        ----------

        value : float
            The value.

        unit : str
            The unit.

        string : str
            A string representation of the angle.

        Examples
        --------

        Create an object representing an angle of 3.14 radians then
        print the length in degrees.

        >>> import BioSimSpace as BSS
        >>> length = BSS.Types.Angle(3.14, "R")
        >>> print(length.degrees())

        The same as above, except passing a string representation of the
        angle to the constructor.

        >>> import BioSimSpace as BSS
        >>> length = BSS.Types.Angle("3.14 R")
        >>> print(length.degrees())

        The string matching is extremeley flexible, so all of the following
        would be valid arguments: "3.14 R", "3.14 radians", "314e-2 Radians".
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

    def radians(self):
        """
        Return the angle in radians.

        Returns
        -------

        angle : :class:`Angle <BioSimSpace.Types.Angle>`
            The angle in radians.
        """
        return Angle(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.radian),
            "RADIAN",
        )

    def degrees(self):
        """
        Return the angle in degrees.

        Returns
        -------

        angle : :class:`Angle <BioSimSpace.Types.Angle>`
            The angle in degrees.
        """
        return Angle(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.degree),
            "DEGREE",
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

        angle : :class:`Angle <BioSimSpace.Types.Angle>`
            The length in the default unit of radians.
        """
        if mag is None:
            return self.radians()
        else:
            return Angle(mag, "RADIAN")

    def _convert_to(self, unit):
        """
        Return the angle in a different unit.

        Parameters
        ----------

        unit : str
            The unit to convert to.

        Returns
        -------

        angle : :class:`Angle <BioSimSpace.Types.Angle>`
            The angle in the specified unit.
        """
        if unit == "RADIAN":
            return self.radians()
        elif unit == "DEGREE":
            return self.degrees()
        else:
            raise ValueError(
                "Supported units are: '%s'" % list(self._supported_units.keys())
            )

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Strip any "S" characters.
        unit = unit.replace("S", "")

        # Strip "EGREE".
        unit = unit.replace("EGREE", "")

        # Strip "EG".
        unit = unit.replace("EG", "")

        # Strip "ADIAN".
        unit = unit.replace("ADIAN", "")

        # Strip "AD".
        unit = unit.replace("AD", "")

        # Check that the unit is supported.
        if unit in self._supported_units:
            return unit
        elif unit in self._abbreviations:
            return self._abbreviations[unit]
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

        # First, handle plurals and abbreviations.
        unit = unit.replace("radians", "rad")
        unit = unit.replace("radian", "rad")
        unit = unit.replace("rads", "rad")

        # Now convert back to correct format.
        unit = unit.replace("rad", "radian")

        # Convert powers. (Limited selection, for now.)
        unit = unit.replace("radian2", "(radian*radian)")
        unit = unit.replace("radian3", "(radian*radian*radian)")
        unit = unit.replace("degree2", "(degree*degree)")
        unit = unit.replace("degree3", "(degree*degree*degree)")
        unit = unit.replace("radian-1", "(1/(radian))")
        unit = unit.replace("radian-2", "(1/(radian*radian))")
        unit = unit.replace("radian-3", "(1/(radian*radian*radian))")
        unit = unit.replace("degree-1", "(1/(degree))")
        unit = unit.replace("degree-2", "(1/(degree*degree))")
        unit = unit.replace("degree-3", "(1/(degree*degree*degree))")

        return unit
