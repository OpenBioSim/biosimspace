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

"""A coordinate (position vector) type."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Coordinate"]

from ._length import Length as _Length
from ._vector import Vector as _Vector


class Coordinate:
    """A coordinate (position vector)."""

    def __init__(self, x, y, z):
        """
        Constructor.

        Parameters
        ----------

        x : :class: `Length <BioSimSpace.Types.Length>`
            The x position.

        y : :class: `Length <BioSimSpace.Types.Length>`
            The y position.

        z : :class: `Length <BioSimSpace.Types.Length>`
            The z position.
        """

        if not isinstance(x, _Length):
            raise TypeError("'x' must be of type 'BioSimSpace.Types.Length'")

        if not isinstance(y, _Length):
            raise TypeError("'y' must be of type 'BioSimSpace.Types.Length'")

        if not isinstance(z, _Length):
            raise TypeError("'z' must be of type 'BioSimSpace.Types.Length'")

        # Set the vector.
        self._vector = _Vector(
            x.angstroms().value(), y.angstroms().value(), z.angstroms().value()
        )

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "(%s, %s, %s)" % (self.x(), self.y(), self.z())

    def __repr__(self):
        """Return a human readable string representation of the object."""
        return "(%s, %s, %s)" % (self.x(), self.y(), self.z())

    def __pos__(self):
        """Unary + operator."""
        return Coordinate(self.x(), self.y(), self.z())

    def __neg__(self):
        """Unary - operator."""
        return Coordinate(-self.x(), -self.y(), -self.z())

    def __add__(self, other):
        """
        Addition operator.

        Parameters
        ----------

        other : :class: `Coordinate <BioSimSpace.Types.Coordinate>`, \
                :class: `Length <BioSimSpace.Types.Length>`

        Return
        ------

        result : :class: `Coordinate <BioSimSpace.Types.Coordinate>`
            The sum of the two coordinates.
        """
        if isinstance(other, Coordinate):
            return self.fromVector(self._vector + other._vector, _Length(1, "A"))

        elif isinstance(other, _Length):
            vector = self._vector + _Vector(
                other.angstroms().value(),
                other.angstroms().value(),
                other.angstroms().value(),
            )
            return self.fromVector(vector, _Length(1, "A"))

        else:
            raise TypeError(
                "unsupported operand type(s) for -: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def __sub__(self, other):
        """
        Subtraction operator.

        Parameters
        ----------

        other : :class: `Coordinate <BioSimSpace.Types.Coordinate>`, \
                :class: `Length <BioSimSpace.Types.Length>`
            Another coordinate or a length.

        Return
        ------

        result : :class: `Coordinate <BioSimSpace.Types.Coordinate>`
            The difference of the two coordinates.
        """
        if isinstance(other, Coordinate):
            return self.fromVector(self._vector - other._vector, _Length(1, "A"))

        elif isinstance(other, _Length):
            vector = self._vector - _Vector(
                other.angstroms().value(),
                other.angstroms().value(),
                other.angstroms().value(),
            )
            return self.fromVector(vector, _Length(1, "A"))

        else:
            raise TypeError(
                "unsupported operand type(s) for -: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def __mul__(self, other):
        """Multiplication operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Only support multiplication by float.
        if isinstance(other, float):
            return self.fromVector(other * self._vector, _Length(1, "A"))

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
            return self.fromVector(self._vector / other, _Length(1, "A"))

        else:
            raise TypeError(
                "unsupported operand type(s) for /: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def x(self):
        """
        Return the x component of the coordinate.

        Returns
        -------

        x : :class: `Length <BioSimSpace.Types.Length>`
            The x component of the coordinate.
        """
        return _Length(self._vector.x(), "A")

    def y(self):
        """
        Return the y component of the coordinate.

        Returns
        -------

        y : :class: `Length <BioSimSpace.Types.Length>`
            The y component of the coordinate.
        """
        return _Length(self._vector.y(), "A")

    def z(self):
        """
        Return the z component of the coordinate.

        Returns
        -------

        z : :class: `Length <BioSimSpace.Types.Length>`
            The z component of the coordinate.
        """
        return _Length(self._vector.z(), "A")

    def to_vector(self):
        """
        Convert to a unitless BioSimSpace.Types.Vector object.

        Returns
        -------

        vector : :class: `Vector <BioSimSpace.Types.Vector>`
            The unitless vector.
        """
        return self._vector

    @staticmethod
    def fromVector(vector, unit):
        """
        Convert from a unitless BioSimSpace.Types.Vector object.

        Returns
        -------

        vector : :class: `Vector <BioSimSpace.Types.Vector>`
            The unitless vector.

        unit : :class: `Length <BioSimSpace.Types.Length>`
            The coordinate unit.
        """
        if not isinstance(vector, _Vector):
            raise TypeError("'vector' must be of type 'BioSimSpace.Types.Vector'")

        if not isinstance(unit, _Length):
            raise TypeError("'unit' must be of type 'BioSimSpace.Types.Length'")

        return Coordinate(vector.x() * unit, vector.y() * unit, vector.z() * unit)

    @staticmethod
    def _from_sire_vector(vector):
        """
        Create a coordinate from a Sire.Maths.Vector object.

        Parameters
        ----------

        vector : Sire.Maths.Vector
            The Sire Vector object.

        Returns
        -------

        coordinate : :class: `Coordinate <BioSimSpace.Types.Coordinate>`
            A BioSimSpace Coordinate object.
        """
        # Create a new Coordinate using the x, y, z components.
        return Coordinate(
            _Length(vector.x().value(), "A"),
            _Length(vector.y().value(), "A"),
            _Length(vector.z().value(), "A"),
        )
