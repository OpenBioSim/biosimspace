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

"""A temperature type."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Temperature"]

from sire.legacy import Units as _SireUnits

from ._type import Type as _Type


class Temperature(_Type):
    """A temperature type."""

    # A list of the supported Sire unit names.
    _sire_units = ["kelvin", "celsius", "fahrenheit"]

    # Dictionary of allowed units.
    _supported_units = {
        "KELVIN": _SireUnits.kelvin,
        "CELSIUS": _SireUnits.celsius,
        "FAHRENHEIT": _SireUnits.fahrenheit,
    }

    # Map unit abbreviations to the full name.
    _abbreviations = {"K": "KELVIN", "C": "CELSIUS", "F": "FAHRENHEIT"}

    # Print formatting.
    _print_format = {"KELVIN": "K", "CELSIUS": "C", "FAHRENHEIT": "F"}

    # Documentation strings.
    _doc_strings = {
        "KELVIN": "A temperature in Kelvin.",
        "CELSIUS": "A temperature in Celsius.",
        "FAHRENHEIT": "A temperature in Fahrenheit.",
    }

    # Null type unit for avoiding issue printing configargparse help.
    _default_unit = "KELVIN"

    # The dimension mask.
    _dimensions = tuple(list(_supported_units.values())[0].dimensions())

    def __init__(self, *args):
        """
        Constructor.

        ``*args`` can be a value and unit, or a string representation
        of the temperature, e.g. "298 K".

        Parameters
        ----------

        value : float
            The value.

        unit : str
            The unit.

        string : str
            A string representation of the temperature.

        Examples
        --------

        Create an object representing a temperature of 298 Kelvin then
        print the temperature in Celsius.

        >>> import BioSimSpace as BSS
        >>> temperature = BSS.Types.Temperature(298, "K")
        >>> print(temperature.celsius())

        The same as above, except passing a string representation of the
        temperature to the constructor.

        >>> import BioSimSpace as BSS
        >>> time = BSS.Types.Temperature("298 K")
        >>> print(temperature.celsius())

        The string matching is extremeley flexible, so all of the following
        would be valid arguments: "298 K", "298 kelvin", "2.98e2 k".
        """

        # Call the base class constructor.
        super().__init__(*args)

    def __add__(self, other):
        """Addition operator."""

        from ..Units import allow_offset

        # Addition of another object of the same type.
        if isinstance(other, self):
            # The temperatures have the same unit.
            if self._unit == other._unit:
                if self._unit != "KELVIN":
                    if not allow_offset:
                        raise ValueError(
                            "Ambiguous operation with offset unit: '%s'" % self._unit
                        )
                    else:
                        # Add the value in the original unit.
                        mag = self._value + other._value

                        # Return a new object of the same type with the original unit.
                        return Temperature(mag, self._unit)
                else:
                    return super().__add__(other)
            else:
                if not allow_offset:
                    raise ValueError(
                        "Ambiguous operation with offset unit: '%s'" % self._unit
                    )
                else:
                    # Left-hand operand takes precedence.
                    mag = self._value + other._convert_to(self._unit).value()

                    # Return a new object of the same type with the original unit.
                    return Temperature(mag, self._unit)

        # Addition of a string.
        elif isinstance(other, str):
            temp = self._from_string(other)
            return self + temp

        else:
            raise TypeError(
                "unsupported operand type(s) for +: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def __sub__(self, other):
        """Subtraction operator."""

        from ..Units import allow_offset

        # Subtraction of another object of the same type.
        if isinstance(other, self):
            # The temperatures have the same unit.
            if self._unit == other._unit:
                if self._unit != "KELVIN":
                    if not allow_offset:
                        raise ValueError(
                            "Ambiguous operation with offset unit: '%s'" % self._unit
                        )
                    else:
                        # Subtract the value in the original unit.
                        mag = self._value - other._value

                        # Return a new object of the same type with the original unit.
                        return Temperature(mag, self._unit)
                else:
                    return super().__sub__(other)
            else:
                if not allow_offset:
                    raise ValueError(
                        "Ambiguous operation with offset unit: '%s'" % self._unit
                    )
                else:
                    # Left-hand operand takes precedence.
                    mag = self._value - other._convert_to(self._unit).value()

                    # Return a new object of the same type with the original unit.
                    return Temperature(mag, self._unit)

        # Addition of a string.
        elif isinstance(other, str):
            temp = self._from_string(other)
            return self - temp

        else:
            raise TypeError(
                "unsupported operand type(s) for -: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def __mul__(self, other):
        """Multiplication operator."""

        if self._unit != "KELVIN":
            from ..Units import allow_offset

            if not allow_offset:
                raise ValueError(
                    "Ambiguous operation with offset unit: '%s'" % self._unit
                )

            else:
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
                    # Multiply value.
                    mag = self._value * other

                    # Return a new object of the same type with the original unit.
                    return Temperature(mag, self._unit)

                # Multiplication by another type.
                elif isinstance(other, _Type):
                    from ._general_unit import GeneralUnit as _GeneralUnit

                    return _GeneralUnit(self._to_sire_unit() * other._to_sire_unit())

                else:
                    raise TypeError(
                        "unsupported operand type(s) for *: '%s' and '%s'"
                        % (self.__class__.__qualname__, other.__class__.__qualname__)
                    )

        else:
            return super().__mul__(other)

    def __rmul__(self, other):
        """Multiplication operator."""

        # Multiplication is commutative: a*b = b*a
        return self.__mul__(other)

    def __truediv__(self, other):
        """Division operator."""

        if self._unit != "KELVIN":
            from ..Units import allow_offset

            if not allow_offset:
                raise ValueError(
                    "Ambiguous operation with offset unit: '%s'" % self._unit
                )

            else:
                # Convert int to float.
                if type(other) is int:
                    other = float(other)

                # Float division.
                if isinstance(other, float):
                    # Divide value.
                    mag = self._value / other

                    # Return a new object of the same type with the original unit.
                    return Temperature(mag, self._unit)

                # Division by another object of the same type.
                elif isinstance(other, self):
                    return self._value / other._convert_to(self._unit).value()

                # Division by another type.
                elif isinstance(other, _Type):
                    from ._general_unit import GeneralUnit as _GeneralUnit

                    return _GeneralUnit(self._to_sire_unit() / other._to_sire_unit())

                # Division by a string.
                elif isinstance(other, str):
                    obj = self._from_string(other)
                    return self / obj

                else:
                    raise TypeError(
                        "unsupported operand type(s) for /: '%s' and '%s'"
                        % (self.__class__.__qualname__, other.__class__.__qualname__)
                    )
        else:
            return super().__truediv__(other)

    def _kelvin(self):
        """Return the value of the temperature in Kelvin."""
        return (self._value * self._supported_units[self._unit]).value()

    def kelvin(self):
        """
        Return the temperature in Kelvin.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The temperature in Kelvin.
        """
        return Temperature(
            (self._value * self._supported_units[self._unit]).value(), "KELVIN"
        )

    def celsius(self):
        """
        Return the temperature in Celsius.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The temperature in Celsius.
        """
        return Temperature(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.celsius),
            "CELSIUS",
        )

    def fahrenheit(self):
        """
        Return the temperature in Fahrenheit.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The temperature in Fahrenheit.
        """
        return Temperature(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.fahrenheit),
            "FAHRENHEIT",
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

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The temperature in the default unit of Kelvin.
        """
        if mag is None:
            return self.kelvin()
        else:
            return Temperature(mag, "KELVIN")

    def _convert_to(self, unit):
        """
        Return the temperature in a different unit.

        Parameters
        ----------

        unit : str
            The unit to convert to.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The temperature in the specified unit.
        """
        if unit == "KELVIN":
            return self.kelvin()
        elif unit == "CELSIUS":
            return self.celsius()
        elif unit == "FAHRENHEIT":
            return self.fahrenheit()
        else:
            raise ValueError(
                "Supported units are: '%s'" % list(self._supported_units.keys())
            )

    @classmethod
    def _validate_unit(cls, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Strip all instances of "DEGREES", "DEGREE", "DEGS", & "DEG".
        unit = unit.replace("DEGREES", "")
        unit = unit.replace("DEGREE", "")
        unit = unit.replace("DEGS", "")
        unit = unit.replace("DEG", "")

        # Check that the unit is supported.
        if unit in cls._supported_units:
            return unit
        elif unit in cls._abbreviations:
            return cls._abbreviations[unit]
        elif len(unit) == 0:
            raise ValueError(f"Unit is not given. You must supply the unit.")
        else:
            raise ValueError(
                "Unsupported unit '%s'. Supported units are: '%s'"
                % (unit, list(cls._supported_units.keys()))
            )

    def _to_sire_unit(self):
        """
        Return the internal Sire Unit object to which this type corresponds.

        Returns
        -------

        sire_unit : sire.units.GeneralUnit
            The internal Sire Unit object that is being wrapped.
        """
        return self.kelvin().value() * _SireUnits.kelvin

    @classmethod
    def _from_sire_unit(cls, sire_unit):
        """
        Convert from a Sire Units object.

        Parameters
        ----------

        sire_unit : sire.units.GeneralUnit, sire.units.Celsius, sire.units.Fahrenheit
            The temperature as a Sire Units object.
        """

        if isinstance(sire_unit, _SireUnits.GeneralUnit):
            # Create a mask for the dimensions of the object.
            dimensions = (
                sire_unit.MASS(),
                sire_unit.LENGTH(),
                sire_unit.TIME(),
                sire_unit.CHARGE(),
                sire_unit.TEMPERATURE(),
                sire_unit.QUANTITY(),
                sire_unit.ANGLE(),
            )

            # Make sure the dimensions match.
            if dimensions != cls._dimensions:
                raise ValueError(
                    "The dimensions of the passed 'sire_unit' are incompatible with "
                    f"'{cls.__name__}'"
                )

            # Get the value in the default Sire unit for this type.
            value = sire_unit.to(cls._supported_units[cls._default_unit])

            # Return an object of this type using the value and unit.
            return cls(value, cls._default_unit)

        elif isinstance(sire_unit, (_SireUnits.Celsius, _SireUnits.Fahrenheit)):
            # Return an object of this type using the value and unit.
            return cls(sire_unit.value(), cls._default_unit)

        else:
            raise TypeError(
                "'sire_unit' must be of type 'sire.units.GeneralUnit', "
                "'sire.units.Celsius', or 'sire.units.Fahrenheit'"
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

        # Convert everything to Kelvin.
        unit = unit.replace("celsius", "274.15*kelvin")
        unit = unit.replace("fahrenheit", "255.9278*kelvin")

        # Convert plural to singular.
        unit = unit.replace("kelvins", "kelvin")

        # Convert powers. (Just 2nd and third for now.)
        unit = unit.replace("kelvin2", "(kelvin*kelvin)")
        unit = unit.replace("kelvin3", "(kelvin*kelvin*kelvin)")
        unit = unit.replace("kelvin-1", "(1/kelvin)")
        unit = unit.replace("kelvin-2", "(1/(kelvin*kelvin))")
        unit = unit.replace("kelvin-3", "(1/(kelvin*kelvin*kelvin))")

        return unit
