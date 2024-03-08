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

"""A length type."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Length"]

from sire.legacy import Units as _SireUnits

from ._type import Type as _Type


class Length(_Type):
    """A length type."""

    # A list of the supported Sire unit names.
    _sire_units = [
        "meter",
        "centimeter",
        "millimeter",
        "nanometer",
        "angstrom",
        "picometer",
    ]

    # Dictionary of allowed units.
    _supported_units = {
        "METER": _SireUnits.meter,
        "CENTIMETER": _SireUnits.centimeter,
        "MILLIMETER": _SireUnits.millimeter,
        "NANOMETER": _SireUnits.nanometer,
        "ANGSTROM": _SireUnits.angstrom,
        "PICOMETER": _SireUnits.picometer,
    }

    # Map unit abbreviations to the full name.
    _abbreviations = {
        "M": "METER",
        "CM": "CENTIMETER",
        "MM": "MILLIMETER",
        "NM": "NANOMETER",
        "A": "ANGSTROM",
        "PM": "PICOMETER",
    }

    # Print format.
    _print_format = {
        "METER": "m",
        "CENTIMETER": "cm",
        "MILLIMETER": "mm",
        "NANOMETER": "nm",
        "ANGSTROM": "A",
        "PICOMETER": "pm",
    }

    # Documentation strings.
    _doc_strings = {
        "METER": "A length in meters.",
        "CENTIMETER": "A length in centimeters.",
        "MILLIMETER": "A length in millimeters.",
        "NANOMETER": "A length in nanometers.",
        "ANGSTROM": "A length in Angstrom.",
        "PICOMETER": "A length in picometers.",
    }

    # Null type unit for avoiding issue printing configargparse help.
    _default_unit = "ANGSTROM"

    # The dimension mask.
    _dimensions = tuple(list(_supported_units.values())[0].dimensions())

    def __init__(self, *args):
        """
        Constructor.

        ``*args`` can be a value and unit, or a string representation
        of the length, e.g. "12 Angstrom".

        Parameters
        ----------

        value : float
            The value.

        unit : str
            The unit.

        string : str
            A string representation of the length.

        Examples
        --------

        Create an object representing a length of 148.6 Angstrom then
        print the length in nanometers.

        >>> import BioSimSpace as BSS
        >>> length = BSS.Types.Length(148.6, "A")
        >>> print(length.nanometers())

        The same as above, except passing a string representation of the
        length to the constructor.

        >>> import BioSimSpace as BSS
        >>> length = BSS.Types.Length("148.6 A")
        >>> print(length.nanometers())

        The string matching is extremeley flexible, so all of the following
        would be valid arguments: "148.6 A", "148.6 angstrom", "1.48e2 Angstrom".
        """

        # Call the base class constructor.
        super().__init__(*args)

    def __mul__(self, other):
        """Multiplication operator."""

        # Handle containers by converting each item in the container to
        # this type.
        if isinstance(other, list):
            return [self.__mul__(item) for item in other]
        if isinstance(other, tuple):
            return tuple([self.__mul__(item) for item in other])

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Multiplication by float.
        if isinstance(other, float):
            mag = self._value * other
            return Length(mag, self._unit)

        # Multiplication by another Length.
        elif isinstance(other, Length):
            mag = self.angstroms().value() * other.angstroms().value()
            return _Area(mag, "A2")

        # Multiplication by an Area.
        elif isinstance(other, _Area):
            mag = self.angstroms().value() * other.angstroms2().value()
            return _Volume(mag, "A3")

        # Multiplication by another type.
        elif isinstance(other, _Type):
            from ._general_unit import GeneralUnit as _GeneralUnit

            return _GeneralUnit(self._to_sire_unit() * other._to_sire_unit())

        # Multiplication by a string.
        elif isinstance(other, str):
            try:
                length = Length(other)
                return self * length
            except:
                try:
                    area = _Area(other)
                    return self * area
                except:
                    raise ValueError(
                        "Could not convert the string to a 'BioSimSpace.Length' "
                        "or 'BioSimSpace.Area' type."
                    )
        else:
            raise TypeError(
                "unsupported operand type(s) for *: '%s' and '%s'"
                % (type(self), type(other))
            )

    def __rmul__(self, other):
        """Multiplication operator."""

        # Multiplication is commutative: a*b = b*a
        return self.__mul__(other)

    def meters(self):
        """
        Return the length in meters.

        Returns
        -------

        length : :class:`Length <BioSimSpace.Types.Length>`
            The length in meters.
        """
        return Length(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.meter),
            "METER",
        )

    def centimeters(self):
        """
        Return the length in centimeters.

        Returns
        -------

        length : :class:`Length <BioSimSpace.Types.Length>`
            The length in centimeters.
        """
        return Length(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.centimeter),
            "CENTIMETER",
        )

    def millimeters(self):
        """
        Return the length in millimeters.

        Returns
        -------

        length : :class:`Length <BioSimSpace.Types.Length>`
            The length in millimeters.
        """
        return Length(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.millimeter),
            "MILLIMETER",
        )

    def nanometers(self):
        """
        Return the length in nanometers.

        Returns
        -------

        length : :class:`Length <BioSimSpace.Types.Length>`
            The length in nanometers.
        """
        return Length(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.nanometer),
            "NANOMETER",
        )

    def angstroms(self):
        """
        Return the length in angstroms.

        Returns
        -------

        length : :class:`Length <BioSimSpace.Types.Length>`
            The length in angstrom.
        """
        return Length(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.angstrom),
            "ANGSTROM",
        )

    def picometers(self):
        """
        Return the length in picometers.

        Returns
        -------

        length : :class:`Length <BioSimSpace.Types.Length>`
            The length in picometers.
        """
        return Length(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.picometer),
            "PICOMETER",
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

        length : :class:`Length <BioSimSpace.Types.Length>`
            The length in the default unit of Angstrom.
        """
        if mag is None:
            return self.angstroms()
        else:
            return Length(mag, "ANGSTROM")

    def _convert_to(self, unit):
        """
        Return the length in a different unit.

        Parameters
        ----------

        unit : str
            The unit to convert to.

        Returns
        -------

        length : :class:`Length <BioSimSpace.Types.Length>`
            The length in the specified unit.
        """
        if unit == "METER":
            return self.meters()
        elif unit == "CENTIMETER":
            return self.centimeters()
        elif unit == "MILLIMETER":
            return self.millimeters()
        elif unit == "NANOMETER":
            return self.nanometers()
        elif unit == "ANGSTROM":
            return self.angstroms()
        elif unit == "PICOMETER":
            return self.picometers()
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

        # Fix for ANGSTROM (since it contains an "S").
        if unit[0:3] == "ANG":
            unit = "ANGS" + unit[3:]

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

        unit = unit.replace("nm", "nanometer")
        unit = unit.replace("pm", "picometer")

        # Convert powers.
        unit = unit.replace("angstrom-1", "(1/angstrom)")
        unit = unit.replace("picometer-1", "(1/picometer)")
        unit = unit.replace("nanometer-1", "(1/nanometer)")
        unit = unit.replace("millimeter-1", "(1/millimeter)")
        unit = unit.replace("centimeter-1", "(1/centimeter)")
        unit = unit.replace("meter-1", "(1/meter)")

        return unit


# Import at bottom of module to avoid circular dependency.
from ._area import Area as _Area
from ._volume import Volume as _Volume
