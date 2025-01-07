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

"""A base class for physical unit types."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Type"]

import math as _math
import re as _re

from sire.legacy import Units as _SireUnits

from .._Exceptions import IncompatibleError as _IncompatibleError


class Type:
    """A base class for custom types."""

    def __init__(self, *args):
        """
        Constructor.

        ``*args`` can be a value and unit, or a string representation
        of the type.

        Parameters
        ----------

        value : float
            The value.

        unit : str
            The unit.

        string : str
            A string representation of the type.
        """

        # Don't allow user to create an instance of this base class.
        if type(self) is Type:
            raise Exception("<Type> must be subclassed.")

        # The user has passed a value and a unit.
        if len(args) > 1:
            value = args[0]
            unit = args[1]

            # Check for a GeneralUnit value.
            if isinstance(value, _SireUnits.GeneralUnit):
                # Try to convert the Sire GeneralUnit to an object of this type.
                temp = self._from_sire_unit(args[0])

                # Validate the unit.
                self._unit = self._validate_unit(unit)

                # Convert to the desired unit.
                self._value = temp._convert_to(self._unit).value()

            else:
                if hasattr(value, "to_default"):
                    value = value.to_default()

                # Check that the value is valid.
                if type(value) is int:
                    self._value = float(value)
                elif isinstance(value, float):
                    self._value = value
                else:
                    raise TypeError("'value' must be of type 'int' or 'float'")

                # Check that the unit is supported.
                self._unit = self._validate_unit(unit)

        elif len(args) == 1:
            # The user has passed a Sire Unit object.
            if isinstance(
                args[0],
                (_SireUnits.GeneralUnit, _SireUnits.Celsius, _SireUnits.Fahrenheit),
            ):
                # Try to convert the Sire GeneralUnit to an object of this type.
                temp = self._from_sire_unit(args[0])
                self._value = temp._value
                self._unit = temp._unit

            # The user has passed a string representation of the type.
            elif isinstance(args[0], str):
                # Convert the string to an object of this type.
                obj = self._from_string(args[0])

                # Store the value and unit.
                self._value = obj._value
                self._unit = obj._unit

            else:
                raise TypeError(
                    "__init__() missing positional argument(s): "
                    "'string' or 'sire.units.GeneralUnit'"
                )

        # No arguments.
        else:
            raise TypeError(
                "__init__() missing positional argument(s): 'value' and 'unit', "
                "or 'string' or 'sire.units.GeneralUnit'"
            )

        # Set the documentation string.
        self.__doc__ = self._doc_strings[self._unit]

    def __str__(self):
        """Return a human readable string representation of the object."""
        if abs(self._value) > 1e4 or abs(self._value) < 1e-4:
            return "%.4e %s" % (self._value, self._print_format[self._unit])
        else:
            return "%5.4f %s" % (self._value, self._print_format[self._unit])

    def __repr__(self):
        """Return a human readable string representation of the object."""
        return self.__str__()

    def __pos__(self):
        """Unary + operator."""
        return type(self)(self.value(), self.unit())

    def __neg__(self):
        """Unary - operator."""
        return type(self)(-self.value(), self.unit())

    def __add__(self, other):
        """Addition operator."""

        # Addition of another object of the same type.
        if type(other) is type(self):
            # Add the values in a common unit.
            val = self._to_default_unit().value() + other._to_default_unit().value()

            # Return a new object of the same type with the original unit.
            return self._to_default_unit(val)._convert_to(self._unit)

        # Addition of a different type with the same dimensions.
        elif isinstance(other, Type) and self._dimensions == other.dimensions:
            # Invert the order of addition.
            return other + self

        # Addition of a string.
        elif isinstance(other, str):
            temp = self._from_string(other)
            return self + temp

        # Addition of a zero-valued integer or float.
        elif isinstance(other, (int, float)) and other == 0:
            return self

        else:
            raise TypeError(
                "unsupported operand type(s) for +: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def __radd__(self, other):
        """Addition operator."""

        # Addition is commutative: a+b = b+a
        return self.__add__(other)

    def __sub__(self, other):
        """Subtraction operator."""

        # Subtraction of another object of the same type.
        if type(other) is type(self):
            # Subtract the values in a common unit.
            val = self._to_default_unit().value() - other._to_default_unit().value()

            # Return a new object of the same type with the original unit.
            return self._to_default_unit(val)._convert_to(self._unit)

        # Subtraction of a different type with the same dimensions.
        elif isinstance(other, Type) and self._dimensions == other.dimensions:
            # Negate other and add.
            return -other + self

        # Subtraction of a string.
        elif isinstance(other, str):
            temp = self._from_string(other)
            return self - temp

        # Subtraction of a zero-valued integer or float.
        elif isinstance(other, (int, float)) and other == 0:
            return self

        else:
            raise TypeError(
                "unsupported operand type(s) for -: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def __rsub__(self, other):
        """Subtraction operator."""

        # Subtraction is not commutative: a-b != b-a
        return -self.__sub__(other)

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
            # Convert to default unit and multiply.
            val = self._to_default_unit().value() * other

            # Return a new object of the same type with the original unit.
            return self._to_default_unit(val)._convert_to(self._unit)

        # Multiplication by another type.
        elif isinstance(other, Type):
            from ._general_unit import GeneralUnit as _GeneralUnit

            return _GeneralUnit(self._to_sire_unit() * other._to_sire_unit())

        else:
            raise TypeError(
                "unsupported operand type(s) for *: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def __rmul__(self, other):
        """Multiplication operator."""

        # Multiplication is commutative: a*b = b*a
        return self.__mul__(other)

    def __pow__(self, other):
        """Power operator."""

        from ._general_unit import GeneralUnit as _GeneralUnit

        return _GeneralUnit(self._to_sire_unit(), no_cast=True) ** other

    def __truediv__(self, other):
        """Division operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Float division.
        if isinstance(other, float):
            # Convert to default unit and divide.
            val = self._to_default_unit().value() / other

            # Return a new object of the same type with the original unit.
            return self._to_default_unit(val)._convert_to(self._unit)

        # Division by another object of the same type.
        elif type(other) is type(self):
            return self._to_default_unit().value() / other._to_default_unit().value()

        # Division by another type.
        elif isinstance(other, Type):
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

    def __rtruediv__(self, other):
        """Reverse division operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Float division.
        if isinstance(other, float):
            from ._general_unit import GeneralUnit as _GeneralUnit

            return _GeneralUnit(other / self._to_sire_unit())

        # Division by another object of the same type.
        elif type(other) is type(self):
            return other._to_default_unit().value() / self._to_default_unit().value()

        # Division by another type.
        elif isinstance(other, Type):
            from ._general_unit import GeneralUnit as _GeneralUnit

            return _GeneralUnit(other._to_sire_unit() / self._to_sire_unit())

        # Division by a string.
        elif isinstance(other, str):
            obj = self._from_string(other)
            return obj / self

        else:
            raise TypeError(
                "unsupported operand type(s) for /: '%s' and '%s'"
                % (other.__class__.__qualname__, self.__class__.__qualname__)
            )

    def __lt__(self, other):
        """Less than operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._to_default_unit().value() < other._to_default_unit().value()

        # Compare with a string.
        elif isinstance(other, str):
            return (
                self._to_default_unit().value()
                < self._from_string(other)._to_default_unit().value()
            )

        else:
            raise TypeError(
                "unorderable types: '%s' < '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def __le__(self, other):
        """Less than or equal to operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._to_default_unit().value() <= other._to_default_unit().value()

        # Compare with a string.
        elif isinstance(other, str):
            return (
                self._to_default_unit().value()
                <= self._from_string(other)._to_default_unit().value()
            )

        else:
            raise TypeError(
                "unorderable types: '%s' <= '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def __eq__(self, other):
        """Equals to operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return _math.isclose(
                self._to_default_unit().value(),
                other._to_default_unit().value(),
            )

        # Compare with a string.
        elif isinstance(other, str):
            return _math.isclose(
                self._to_default_unit().value(),
                self._from_string(other)._to_default_unit().value(),
            )

        else:
            return False

    def __ge__(self, other):
        """Greater than or equal to operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._to_default_unit().value() >= other._to_default_unit().value()

        # Compare with a string.
        elif isinstance(other, str):
            return (
                self._to_default_unit().value()
                >= self._from_string(other)._to_default_unit().value()
            )

        else:
            raise TypeError(
                "unorderable types: '%s' >= '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def __gt__(self, other):
        """Greater than operator."""

        # Compare to another object of the same type.
        if type(other) is type(self):
            return self._to_default_unit().value() > other._to_default_unit().value()

        # Compare with a string.
        elif isinstance(other, str):
            return (
                self._to_default_unit().value()
                > self._from_string(other)._to_default_unit().value()
            )

        else:
            raise TypeError(
                "unorderable types: '%s' > '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def value(self):
        """
        Return the value.

        Returns
        -------

        value : float
            The value of the type.
        """
        return self._value

    def unit(self):
        """
        Return the unit.

        Returns
        -------

        unit : str
            The unit of the type.
        """
        return self._unit

    def to(self, unit):
        """
        Return the value of this general unit in the specified unit.

        Parameters
        ----------

        unit : str
            The unit to convert to.

        Returns
        -------

        value : float
            The value in the specified unit.
        """

        if not isinstance(unit, str):
            raise TypeError("'unit' must be of type 'str'.")

        try:
            return self._value * self._supported_units[self._unit].to(unit)
        except Exception as e:
            raise _IncompatibleError(
                str(e).split("incompatible_error:")[1].split("(call")[0].strip()
            )

    @classmethod
    def dimensions(cls):
        """
        Return the dimensions of this type. This is a tuple
        containing the power in each dimension.

        Returns : (int, int, int, int, int, int)
            The power in each dimension: 'mass', 'length', 'temperature',
            'charge', 'time', 'quantity', and 'angle'.
        """
        return cls._dimensions

    @classmethod
    def mass(cls):
        """
        Return the power in the 'mass' dimension.

        Returns
        -------

        mass : int
            The power in the 'mass' dimension.
        """
        return cls._dimensions[0]

    @classmethod
    def length(cls):
        """
        Return the power in the 'length' dimension.

        Returns
        -------

        length : int
            The power in the 'length' dimension.
        """
        return cls._dimensions[1]

    @classmethod
    def time(cls):
        """
        Return the power in the 'time' dimension.

        Returns
        -------

        time : int
            The power the 'time' dimension.
        """
        return cls._dimensions[2]

    @classmethod
    def charge(cls):
        """
        Return the power in the 'charge' dimension.

        Returns
        -------

        charge : int
            The power in the 'charge' dimension.
        """
        return cls._dimensions[3]

    @classmethod
    def temperature(cls):
        """
        Return the power in the 'temperature' dimension.

        Returns
        -------

        temperature : int
            The power in the 'temperature' dimension.
        """
        return cls._dimensions[4]

    @classmethod
    def quantity(cls):
        """
        Return the power in the 'quantity' dimension.

        Returns
        -------

        quantity : int
            The power in the 'quantity' dimension.
        """
        return cls._dimensions[5]

    @classmethod
    def angle(cls):
        """
        Return the power in the 'angle' dimension.

        Returns
        -------

        angle : int
            The power in the 'angle' dimension.
        """
        return cls._dimensions[6]

    @classmethod
    def _from_string(cls, string):
        """
        Convert a string to an object of the same type.

        Parameters
        ----------

        string : str
            The string to interpret.

        Returns
        -------

        type : :class:`Type <BioSimSpace.Types>`
            The type object.
        """

        if string == "==SUPPRESS==":
            return cls(0, cls._default_unit)

        string_copy = string

        if isinstance(string, str):
            # Strip white space from the string.
            string = string.replace(" ", "")

            # Try to match scientific format.
            match = _re.search(r"(\-?\d+\.?\d*e\-?\d+)(.*)", string, _re.IGNORECASE)

            # Try to match decimal format.
            if match is None:
                match = _re.search(r"(\-?\d+\.?\d*)(.*)", string, _re.IGNORECASE)

                # No matches, raise an error.
                if match is None:
                    raise ValueError(
                        "Could not interpret '%s' as '%s'" % (string_copy, cls.__name__)
                    )

            # Extract the value and unit.
            value, unit = match.groups()

            # Convert the value to a float.
            value = float(value)

            # Create and return a new object.
            return cls(value, unit)

        else:
            raise TypeError("'string' must be of type 'str'")

    def _to_sire_unit(self):
        """
        Return the internal Sire Unit object to which this type corresponds.

        Returns
        -------

        sire_unit : sire.units.GeneralUnit
            The internal Sire Unit object that is being wrapped.
        """
        return self._value * self._supported_units[self._unit]

    @classmethod
    def _from_sire_unit(cls, sire_unit):
        """
        Convert from a Sire GeneralUnit object.

        Parameters
        ----------

        sire_unit : sire.units.GeneralUnit
            A Sire GeneralUnit object.
        """

        if not isinstance(sire_unit, _SireUnits.GeneralUnit):
            raise TypeError("'sire_unit' must be of type 'sire.units.GeneralUnit'")

        # Create a mask for the dimensions of the object.
        dimensions = tuple(sire_unit.dimensions())

        # Make sure that this isn't zero.
        if hasattr(sire_unit, "is_zero"):
            if sire_unit.is_zero():
                dimensions = cls._dimensions

        # Make sure the dimensions match.
        if dimensions != cls._dimensions:
            raise ValueError(
                f"The dimensions of the passed 'sire_unit' {sire_unit} are incompatible with "
                f"'{cls.__name__}'"
            )

        # Get the value in the default Sire unit for this type.
        value = sire_unit.to(cls._supported_units[cls._default_unit])

        # Return an object of this type using the value and unit.
        return cls(value, cls._default_unit)
