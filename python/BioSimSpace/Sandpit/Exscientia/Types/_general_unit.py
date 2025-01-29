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

"""A general unit based type."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["GeneralUnit"]

import math as _math

from sire.legacy.Units import GeneralUnit as _GeneralUnit

from ._base_units import *
from ._type import Type as _Type


class GeneralUnit(_Type):
    """A general unit type."""

    _dimension_chars = [
        "M",  # Mass
        "L",  # Length
        "T",  # Time
        "C",  # Charge
        "t",  # Temperature
        "Q",  # Quantity
        "A",  # Angle
    ]

    def __new__(cls, *args, no_cast=False):
        """
        Constructor.

        ``*args`` can be a value and unit, or a string representation
        of the type, e.g. "30 kcal_per_mol / angstrom squared".

        Parameters
        ----------

        value : float
            The value.

        unit : str
            The unit.

        string : str
            A string representation of the unit type.

        no_cast: bool
            Whether to disable casting to a specific type.
        """

        # This operator may be called when unpickling an object. Catch empty
        # *args by calling __init__ immediately.
        if not args:
            return super(GeneralUnit, cls).__new__(cls)

        value = 1
        _args = list(args)

        # The user has passed a value and a unit.
        if len(args) > 1:
            value = args[0]
            unit = args[1]

            # Check that the value is valid.
            if type(value) is int:
                pass
            elif isinstance(value, float):
                pass
            else:
                raise TypeError("'value' must be of type 'int' or 'float'")

            # Delete the value so we can use the same logic below.
            del _args[0]

        if len(_args) == 1:
            # The user has passed a Sire GeneralUnit.
            if isinstance(_args[0], _GeneralUnit):
                general_unit = _args[0]

            # The user has passed a string representation of the type.
            elif isinstance(_args[0], str):
                # Extract the string.
                string = _args[0]

                # Try to parse and evalute the string.
                general_unit = cls._from_string(string)

                # We have converted to another type, return it immediately.
                if type(general_unit) is not GeneralUnit:
                    return value * general_unit
                else:
                    general_unit = general_unit._sire_unit

            else:
                raise TypeError(
                    "__new__() missing positional argument(s): "
                    "'string' or 'sire.units.GeneralUnit'. "
                    f"Passed argument {_args[0]} is a {type(_args[0])}"
                )

        # No arguments.
        else:
            raise TypeError(
                "__new__() missing positional argument(s): 'value' and 'unit', "
                "or 'string' or 'Sire.Units.GeneralUnit'"
            )

        # Scale the general unit by the value.
        general_unit = value * general_unit

        # Store the dimension mask.
        dimensions = tuple(general_unit.dimensions())

        # This is a dimensionless quantity, return the value as a float.
        if all(x == 0 for x in dimensions):
            return general_unit.value()

        # Check to see if the dimensions correspond to a supported type.
        # If so, return an object of that type.
        if not no_cast and dimensions in _base_dimensions:
            return _base_dimensions[dimensions](general_unit)
        # Otherwise, call __init__()
        else:
            return super(GeneralUnit, cls).__new__(cls)

    def __init__(self, *args, no_cast=False):
        """
        Constructor.

        ``*args`` can be a value and unit, or a string representation
        of the type, e.g. "30 kcal_per_mol / angstrom squared".

        Parameters
        ----------

        value : float
            The value.

        unit : str
            The unit.

        string : str
            A string representation of the unit type.

        no_cast: bool
            Whether to disable casting to a specific type.
        """

        value = 1
        _args = list(args)

        # The user has passed a value and a unit.
        if len(args) > 1:
            value = _args[0]
            unit = _args[1]

            # Check that the value is valid.
            if type(value) is int:
                self._value = float(value)
            elif isinstance(value, float):
                self._value = value
            else:
                raise TypeError("'value' must be of type 'int' or 'float'")

            # Delete the value so we can use the same logic below.
            del _args[0]

        if len(_args) == 1:
            # The user has passed a Sire GeneralUnit.
            if isinstance(_args[0], _GeneralUnit):
                general_unit = _args[0]

            # The user has passed a string representation of the type.
            elif isinstance(_args[0], str):
                # Extract the string.
                string = _args[0]

                # Try to parse and evalute the string.
                general_unit = self._from_string(string)._sire_unit

            else:
                raise TypeError(
                    "__init__() missing positional argument(s): "
                    "'string' or 'sire.units.GeneralUnit'. "
                    f"Passed value {_args[0]} is a {type(_args[0])}"
                )

        # No arguments.
        else:
            raise TypeError(
                "__init__() missing positional argument(s): 'value' and 'unit', "
                "or 'string' or 'sire.units.GeneralUnit'"
            )

        # Need to rescale the value since the units might have been
        # converted to the default for this GeneralUnit.
        self._sire_unit = value * general_unit
        self._value = self._sire_unit.value()

        # Store the dimension mask.
        self._dimensions = tuple(general_unit.dimensions())

        # Create the unit string.
        self._unit = ""
        for x, dim in enumerate(self._dimensions):
            if dim != 0:
                if len(self._unit) > 0:
                    self._unit += " "
                self._unit += self._dimension_chars[x]
                if dim != 1:
                    self._unit += f"{dim}"

    def __str__(self):
        """Return a human readable string representation of the object."""
        return str(self._sire_unit)

    def __repr__(self):
        """Return a human readable string representation of the object."""
        return self.__str__()

    def __pos__(self):
        """Unary + operator."""
        return type(self)(self._sire_unit)

    def __neg__(self):
        """Unary - operator."""
        return type(self)(-self._sire_unit)

    def __add__(self, other):
        """Addition operator."""

        # Addition of another object with the same dimensions.
        if isinstance(other, _Type) and other._dimensions == self._dimensions:
            temp = self._sire_unit + other._to_sire_unit()
            return GeneralUnit(temp)

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

        # Subtraction of another object with the same dimensions.
        if isinstance(other, _Type) and other._dimensions == self._dimensions:
            temp = self._sire_unit - other._to_sire_unit()
            return GeneralUnit(temp)

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

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Float.
        if isinstance(other, float):
            return GeneralUnit(other * self._sire_unit)

        # Another Type.
        elif isinstance(other, _Type):
            # Multipy the Sire unit objects.
            temp = self._sire_unit * other._to_sire_unit()

            # Get the dimension mask.
            dimensions = temp.dimensions()

            # Return as an existing type if the dimensions match.
            try:
                return _base_dimensions[dimensions](temp)
            except:
                return GeneralUnit(temp)

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

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Float division.
        if isinstance(other, float):
            return GeneralUnit(self._sire_unit / other)

        # Division by another Type.
        elif isinstance(other, _Type):
            # Divide the Sire unit objects.
            temp = self._sire_unit / other._to_sire_unit()

            # Create the dimension mask.
            dimensions = (
                temp.ANGLE(),
                temp.CHARGE(),
                temp.LENGTH(),
                temp.MASS(),
                temp.QUANTITY(),
                temp.TEMPERATURE(),
                temp.TIME(),
            )

            # Return as an existing type if the dimensions match.
            try:
                return _base_dimensions[dimensions](temp)
            except:
                return GeneralUnit(temp)

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
            return GeneralUnit(other / self._sire_unit)

        # Division by another Type.
        elif isinstance(other, _Type):
            # Divide the Sire unit objects.
            temp = other._to_sire_unit() / self._sire_unit

            # Create the dimension mask.
            dimensions = (
                temp.ANGLE(),
                temp.CHARGE(),
                temp.LENGTH(),
                temp.MASS(),
                temp.QUANTITY(),
                temp.TEMPERATURE(),
                temp.TIME(),
            )

            # Return as an existing type if the dimensions match.
            try:
                return _base_dimensions[dimensions](temp)
            except:
                return GeneralUnit(temp)

        # Division by a string.
        elif isinstance(other, str):
            obj = self._from_string(other)
            return obj / self

        else:
            raise TypeError(
                "unsupported operand type(s) for /: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def __pow__(self, other):
        """Power operator."""

        if not isinstance(other, (int, float)):
            raise TypeError(
                "unsupported operand type(s) for ^: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

        if other == 0:
            return GeneralUnit(self._sire_unit / self._sire_unit)

        # Convert to float.
        other = float(other)

        # Get the existing unit dimensions.
        dims = self.dimensions()

        # Compute the new dimensions, rounding floats to 16 decimal places.
        new_dims = [round(dim * other, 16) for dim in dims]

        # Make sure the new dimensions are integers.
        if not all(dim.is_integer() for dim in new_dims):
            raise ValueError(
                "The exponent must be a factor of all the unit dimensions."
            )

        # Convert to integers.
        new_dims = [int(dim) for dim in new_dims]

        # Compute the new value.
        value = self.value() ** other

        # Return a new GeneralUnit object.
        return GeneralUnit(_GeneralUnit(value, new_dims))

    def __lt__(self, other):
        """Less than operator."""

        # Compare to another object of the same type and dimensions.
        if isinstance(other, _Type) and other._dimensions == self._dimensions:
            return self._value < other._value

        # Compare with a string.
        elif isinstance(other, str):
            other = self._from_string(other)
            if other._dimensions == self._dimensions:
                return self._value < other._value
            else:
                raise TypeError(
                    "unorderable types: '%s' < '%s'"
                    % (self.__class__.__qualname__, other.__class__.__qualname__)
                )

        else:
            raise TypeError(
                "unorderable types: '%s' < '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def __le__(self, other):
        """Less than or equal to operator."""

        # Compare to another object of the same type and dimensions.
        if isinstance(other, _Type) and other._dimensions == self._dimensions:
            return self._value <= other._value

        # Compare with a string.
        elif isinstance(other, str):
            other = self._from_string(other)
            if other._dimensions == self._dimensions:
                return self._value <= other._value
            else:
                raise TypeError(
                    "unorderable types: '%s' <= '%s'"
                    % (self.__class__.__qualname__, other.__class__.__qualname__)
                )

        else:
            raise TypeError(
                "unorderable types: '%s' <= '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def __eq__(self, other):
        """Equals to operator."""

        # Compare to another object of the same type and dimensions.
        if isinstance(other, _Type) and other._dimensions == self._dimensions:
            return _math.isclose(self._value, other._value)

        # Compare with a string.
        elif isinstance(other, str):
            other = self._from_string(other)
            if other._dimensions == self._dimensions:
                return _math.isclose(self._value, other._value)
            else:
                return False

        else:
            return False

    def __ne__(self, other):
        """Not equals to operator."""

        # Compare to another object of the same type and dimensions.
        if isinstance(other, _Type) and other._dimensions == self._dimensions:
            return not _math.isclose(self._value, other._value)

        # Compare with a string.
        elif isinstance(other, str):
            other = self._from_string(other)
            if other._dimensions == self._dimensions:
                return not _math.isclose(self._value, other._value)
            else:
                return True

        else:
            return True

    def __ge__(self, other):
        """Greater than or equal to operator."""

        # Compare to another object of the same type and dimensions.
        if isinstance(other, _Type) and other._dimensions == self._dimensions:
            return self._value >= other._value

        # Compare with a string.
        elif isinstance(other, str):
            other = self._from_string(other)
            if other._dimensions == self._dimensions:
                return self._value >= other._value
            else:
                raise TypeError(
                    "unorderable types: '%s' >= '%s'"
                    % (self.__class__.__qualname__, other.__class__.__qualname__)
                )

        else:
            raise TypeError(
                "unorderable types: '%s' >= '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def __gt__(self, other):
        """Greater than operator."""

        # Compare to another object of the same type and dimensions.
        if isinstance(other, _Type) and other._dimensions == self._dimensions:
            return self._value > other._value

        # Compare with a string.
        elif isinstance(other, str):
            other = self._from_string(other)
            if other._dimensions == self._dimensions:
                return self._value > other._value
            else:
                raise TypeError(
                    "unorderable types: '%s' > '%s'"
                    % (self.__class__.__qualname__, other.__class__.__qualname__)
                )

        else:
            raise TypeError(
                "unorderable types: '%s' > '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def unit(self):
        """
        Return the powers of the unit in each dimension.

        Returns
        -------

        unit : str
            The powers of the unit in each dimension.
        """
        return self._unit

    def dimensions(self):
        """
        Return the dimensions of this type. This is a tuple
        containing the power in each dimension.

        Returns : (int, int, int, int, int, int)
            The power in each dimension: 'angle', 'charge', 'length',
            'mass', 'quantity', 'temperature', and 'time'.
        """
        return self._dimensions

    def mass(self):
        """
        Return the power of this general unit in the 'mass' dimension.

        Returns
        -------

        mass : int
            The power of the general unit in the 'mass' dimension.
        """
        return self._dimensions[0]

    def length(self):
        """
        Return the power of this general unit in the 'length' dimension.

        Returns
        -------

        length : int
            The power of the general unit in the 'length' dimension.
        """
        return self._dimensions[1]

    def time(self):
        """
        Return the power of this general unit in the 'time' dimension.

        Returns
        -------

        time : int
            The power of the general unit in the 'time' dimension.
        """
        return self._dimensions[2]

    def charge(self):
        """
        Return the power of this general unit in the 'charge' dimension.

        Returns
        -------

        charge : int
            The power of the general unit in the 'charge' dimension.
        """
        return self._dimensions[3]

    def temperature(self):
        """
        Return the power of this general unit in the 'temperature' dimension.

        Returns
        -------

        temperature : int
            The power of the general unit in the 'temperature' dimension.
        """
        return self._dimensions[4]

    def quantity(self):
        """
        Return the power of this general unit in the 'quantity' dimension.

        Returns
        -------

        quantity : int
            The power of the general unit in the 'quantity' dimension.
        """
        return self._dimensions[5]

    def angle(self):
        """
        Return the power of this general unit in the 'angle' dimension.

        Returns
        -------

        angle : int
            The power of the general unit in the 'angle' dimension.
        """
        return self._dimensions[6]

    def _to_sire_unit(self):
        """
        Return the internal Sire Unit object to which this type corresponds.

        Returns
        -------

        sire_unit : Sire.Units.GeneralUnit
            The internal Sire Unit object that is being wrapped.
        """
        return self._sire_unit

    @staticmethod
    def _from_sire_unit(sire_unit):
        """
        Convert from a Sire GeneralUnit object.

        Parameters
        ----------

        sire_unit : Sire.Units.GeneralUnit
            A Sire GeneralUnit object.
        """

        if not isinstance(sire_unit, _GeneralUnit):
            raise TypeError(
                "'sire_unit' must be of type 'sire.units.GeneralUnit'. "
                f"Value {sire_unit} is of type {type(sire_unit)}"
            )

        return GeneralUnit(sire_unit)

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

        string_copy = string

        if isinstance(string, str):
            # Try to parse as a GeneralUnit using Sire's inbuilt unit grammar.
            try:
                from sire import u as _unit

                return GeneralUnit(_unit(string))

            # Try manually parsing the string.
            except:
                # Convert to lower case and strip whitespace.
                string = string.lower().replace(" ", "")

                # Convert powers to common format.
                string = string.replace("squared", "2")
                string = string.replace("**2", "2")
                string = string.replace("^2", "2")
                string = string.replace("cubed", "3")
                string = string.replace("**3", "3")
                string = string.replace("^3", "3")
                string = string.replace("**-1", "-1")
                string = string.replace("^-1", "-1")
                string = string.replace("**-2", "-2")
                string = string.replace("^-2", "-2")
                string = string.replace("**-3", "-3")
                string = string.replace("^-3", "-3")

                for unit in _base_units:
                    string = unit._to_sire_format(string)

                try:
                    # Compile the eval expression to bytecode.
                    code = compile(string, "<string>", "eval")

                    # The bytecode must contain names.
                    if not code.co_names:
                        raise ValueError(
                            f"Could not infer GeneralUnit from string '{string}'"
                        ) from None

                    # Make sure the co_names only contains names within the allowed
                    # sire_units_local dictionary.
                    for name in code.co_names:
                        if name not in _sire_units_locals:
                            raise ValueError(
                                f"Could not infer GeneralUnit from string '{string}'"
                            ) from None

                    general_unit = eval(string, {}, _sire_units_locals)

                    # Create and return a new object.
                    return GeneralUnit(general_unit)

                except Exception as e:
                    raise ValueError(
                        f"Could not infer GeneralUnit from string '{string}'"
                    ) from None

        else:
            raise TypeError("'string' must be of type 'str'")

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
            return self
        else:
            return mag * self
