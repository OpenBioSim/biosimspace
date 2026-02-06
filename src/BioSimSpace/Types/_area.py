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

"""An area type."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Area"]

from sire.legacy import Units as _SireUnits

from ._type import Type as _Type


class Area(_Type):
    """An area type."""

    # A list of the supported Sire unit names.
    _sire_units = ["meter2", "nanometer2", "angstrom2", "picometer2"]

    # Dictionary of allowed units.
    _supported_units = {
        "METER2": _SireUnits.meter2,
        "NANOMETER2": _SireUnits.nanometer2,
        "ANGSTROM2": _SireUnits.angstrom2,
        "PICOMETER2": _SireUnits.picometer2,
    }

    # Map unit abbreviations to the full name.
    _abbreviations = {
        "M2": "METER2",
        "NM2": "NANOMETER2",
        "A2": "ANGSTROM2",
        "PM2": "PICOMETER2",
    }

    # Print format.
    _print_format = {
        "METER2": "m^2",
        "NANOMETER2": "nm^2",
        "ANGSTROM2": "A^2",
        "PICOMETER2": "pm^2",
    }

    # Documentation strings.
    _doc_strings = {
        "METER2": "An area in square meters.",
        "NANOMETER2": "An area in square nanometers.",
        "ANGSTROM2": "An area in square Angstrom.",
        "PICOMETER2": "An area in square picometers.",
    }

    # Null type unit for avoiding issue printing configargparse help.
    _default_unit = "ANGSTROM2"

    # The dimension mask.
    _dimensions = tuple(list(_supported_units.values())[0].dimensions())

    def __init__(self, *args):
        """
        Constructor.

        ``*args`` can be a value and unit, or a string representation
        of the area, e.g. "30 nm^2".

        Parameters
        ----------

        value : float
            The value.

        unit : str
            The unit.

        string : str
            A string representation of the area.

        Examples
        --------

        Create an object representing an area of 30 square nanometers then
        print the area in square Angstrom.

        >>> import BioSimSpace as BSS
        >>> area = BSS.Types.Area(30, "nm^2")
        >>> print(area.angstroms2())

        The same as above, except passing a string representation of the
        area to the constructor.

        >>> import BioSimSpace as BSS
        >>> area = BSS.Types.Area("30 nm^2")
        >>> print(area.angstroms2())

        The string matching is extremeley flexible, so all of the following
        would be valid arguments: "30 nm^2", "30 square nanometers",
        "30 nanometers squared".
        """

        # Call the base class constructor.
        super().__init__(*args)

    def __mul__(self, other):
        """Multiplication operator."""
        from ._length import Length as _Length
        from ._volume import Volume as _Volume

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
            return Area(mag, self._unit)

        # Multiplication by a Length.
        elif isinstance(other, _Length):
            mag = self.angstroms2().value() * other.angstroms().value()
            return _Volume(mag, "A3")

        # Multiplication by another type.
        elif isinstance(other, _Type):
            from ._general_unit import GeneralUnit as _GeneralUnit

            return _GeneralUnit(self._to_sire_unit() * other._to_sire_unit())

        # Multiplication by a string.
        elif isinstance(other, str):
            try:
                length = _Length(other)
                return self * length
            except:
                raise ValueError(
                    "Could not convert the string to a 'BioSimSpace.Types.Length'"
                )

        else:
            raise TypeError(
                "unsupported operand type(s) for *: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def __rmul__(self, other):
        """Multiplication operator."""

        # Multiplication is commutative: a*b = b*a
        return self.__mul__(other)

    def __truediv__(self, other):
        """Division operator."""
        from ._length import Length as _Length

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Float division.
        if isinstance(other, float):
            mag = self._value / other
            return Area(mag, self._unit)

        # Division by another Area.
        elif isinstance(other, Area):
            return self.angstroms2().value() / other.angstroms2().value()

        # Division by a Length.
        elif isinstance(other, _Length):
            mag = self.angstroms2().value() / other.angstroms().value()
            return _Length(mag, "A")

        # Division by another type.
        elif isinstance(other, _Type):
            from ._general_unit import GeneralUnit as _GeneralUnit

            return _GeneralUnit(self._to_sire_unit() / other._to_sire_unit())

        # Division by a string.
        elif isinstance(other, str):
            try:
                length = _Length(other)
                return self / length
            except:
                try:
                    area = Area(other)
                    return self / area
                except:
                    raise ValueError(
                        "Could not convert the string to a "
                        "'BioSimSpace.Types.Length' or a 'BioSimSpace.Types.Area'."
                    )

        else:
            raise TypeError(
                "unsupported operand type(s) for /: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def meters2(self):
        """
        Return the area in square meters.

        Returns
        -------

        area : :class:`Area <BioSimSpace.Types.Area>`
            The area in square meters.
        """
        return Area(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.meter2),
            "METER2",
        )

    def nanometers2(self):
        """
        Return the area in square nanometers.

        Returns
        -------

        area : :class:`Area <BioSimSpace.Types.Area>`
            The area in square nanometers.
        """
        return Area(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.nanometer2),
            "NANOMETER2",
        )

    def angstroms2(self):
        """
        Return the area in square Angstrom.

        Returns
        -------

        area : :class:`Area <BioSimSpace.Types.Area>`
            The area in square Angstrom.
        """
        return Area(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.angstrom2),
            "ANGSTROM2",
        )

    def picometers2(self):
        """
        Return the area in square picometers.

        Returns
        -------

        area : :class:`Area <BioSimSpace.Types.Area>`
            The area in square picometers.
        """
        return Area(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.picometer2),
            "PICOMETER2",
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

        area : :class:`Area <BioSimSpace.Types.Area>`
            The area in the default unit of square Angstrom.
        """
        if mag is None:
            return self.angstroms2()
        else:
            return Area(mag, "ANGSTROM2")

    def _convert_to(self, unit):
        """
        Return the area in a different unit.

        Parameters
        ----------

        unit : str
            The unit to convert to.

        Returns
        -------

        area : :class:`Area <BioSimSpace.Types.Area>`
            The area in the specified unit.
        """
        if unit == "METER2":
            return self.meters2()
        elif unit == "NANOMETER2":
            return self.nanometers2()
        elif unit == "ANGSTROM2":
            return self.angstroms2()
        elif unit == "PICOMETER2":
            return self.picometers2()
        else:
            raise ValueError(
                "Supported units are: '%s'" % list(self._supported_units.keys())
            )

    @classmethod
    def _validate_unit(cls, unit):
        """Validate that the unit is supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Replace any occurrence of squared with 2.
        unit = unit.replace("SQUARED", "2")
        unit = unit.replace("SQUARE", "2")

        # Strip "^" character.
        unit = unit.replace("^", "")

        # Strip "**" characters.
        unit = unit.replace("**", "")

        # Strip any "S" characters.
        unit = unit.replace("S", "")

        # Fix for ANGSTROM (since it contains an "S").
        if unit[0:3] == "ANG":
            unit = "ANGS" + unit[3:]

        # Make sure that the "2" character appears last. This allows the user
        # to write, e.g. "square nm" or "nm squared".
        index = unit.find("2")
        if index != -1:
            unit = unit[0:index] + unit[index + 1 :] + "2"

        # Check that the unit is supported.
        if unit in cls._supported_units:
            return unit
        elif unit in cls._abbreviations:
            return cls._abbreviations[unit]
        else:
            raise ValueError(
                "Supported units are: '%s'" % list(cls._supported_units.keys())
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

        unit = unit.replace("angstroms", "angstrom")
        unit = unit.replace("meters", "meter")
        unit = unit.replace("nm", "nanometer")
        unit = unit.replace("pm", "picometer")

        # Convert powers.
        unit = unit.replace("angstrom-2", "(1/angstrom2)")
        unit = unit.replace("picometer-2", "(1/picometer2)")
        unit = unit.replace("nanometer-2", "(1/nanometer2)")
        unit = unit.replace("meter-2", "(1/meter2)")

        return unit
