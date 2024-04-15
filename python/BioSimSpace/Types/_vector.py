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

"""A general three-vector type."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Vector"]

from sire.legacy.Maths import Vector as _Vector

from ._angle import Angle as _Angle


class Vector:
    """A three-vector."""

    def __init__(self, x, y, z):
        """
        Constructor.

        Parameters
        ----------

        x : float
            The x component of the vector.

        y : float
            The y component of the vector.

        z : float
            The z component of the vector.
        """

        try:
            x = float(x)
        except:
            raise TypeError("'x' must be of type 'float'")

        try:
            y = float(y)
        except:
            raise TypeError("'y' must be of type 'float'")

        try:
            z = float(z)
        except:
            raise TypeError("'z' must be of type 'float'")

        # Set the Sire Vector.
        self._sire_object = _Vector(x, y, z)

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "(%s, %s, %s)" % (self.x(), self.y(), self.z())

    def __repr__(self):
        """Return a human readable string representation of the object."""
        return "(%s, %s, %s)" % (self.x(), self.y(), self.z())

    def __pos__(self):
        """Unary + operator."""
        return Vector(self.x(), self.y(), self.z())

    def __neg__(self):
        """Unary - operator."""
        return Vector(-self.x(), -self.y(), -self.z())

    def __add__(self, other):
        """
        Addition operator.

        Parameters
        ----------

        other : :class: `Vector <BioSimSpace.Types.Vector>`
            Another vector.

        Return
        ------

        result : :class: `Vector <BioSimSpace.Types.Vector>`
            The sum of the two vectors.
        """
        if not isinstance(other, Vector):
            raise TypeError(
                "unsupported operand type(s) for +: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

        return self._from_sire_vector(self._sire_object + other._sire_object)

    def __sub__(self, other):
        """
        Subtraction operator.

        Parameters
        ----------

        other : :class: `Vector <BioSimSpace.Types.Vector>`
            Another vector.

        Return
        ------

        result : :class: `Vector <BioSimSpace.Types.Vector>`
            The difference of the two vectors.
        """
        if not isinstance(other, Vector):
            raise TypeError(
                "unsupported operand type(s) for -: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

        return self._from_sire_vector(self._sire_object - other._sire_object)

    def __mul__(self, other):
        """Multiplication operator."""

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Only support multiplication by float.
        if isinstance(other, float):
            # Return a new vector multiplied by other.
            return self._from_sire_vector(other * self._sire_object)

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
            # Return a new vector divided by other.
            return self._from_sire_vector(self._sire_object / other)

        else:
            raise TypeError(
                "unsupported operand type(s) for /: '%s' and '%s'"
                % (self.__class__.__qualname__, other.__class__.__qualname__)
            )

    def dot(self, other):
        """
        Return the dot (scalar) product with the other vector.

        Parameters
        ----------

        other : :class: `Vector <BioSimSpace.Types.Vector>`
            Another vector.

        Returns
        -------

        result : float
            The scalar product.
        """
        if not isinstance(other, Vector):
            raise TypeError("'other' must be of type 'BioSimSpace.Types.Vector'")

        return self._sire_object.dot(self._sire_object, other._sire_object)

    def cross(self, other):
        """
        Return the cross product with the other vector.

        Parameters
        ----------

        other : :class: `Vector <BioSimSpace.Types.Vector>`
            Another vector.

        Returns
        -------

        result : :class: `Vector <BioSimSpace.Types.Vector>`
            The cross product.
        """
        if not isinstance(other, Vector):
            raise TypeError("'other' must be of type 'BioSimSpace.Types.Vector'")

        x = self.y() * other.z() - self.z() * other.y()
        y = self.z() * other.x() - self.x() * other.z()
        z = self.x() * other.y() - self.y() * other.x()

        # Create a new Vector using the x, y, z components.
        return Vector(x, y, z)

    def angle(self, other):
        """
        Return the angle between this and the other vector.

        Parameters
        ----------

        other : :class: `Vector <BioSimSpace.Types.Vector>`
            Another vector.

        Returns
        -------

        angle : :class: `Angle <BioSimSpace.Types.Angle>`
            The angle between the two vectors.
        """
        if not isinstance(other, Vector):
            raise TypeError("'other' must be of type 'BioSimSpace.Types.Vector'")

        # Calculate the angle.
        angle = self._sire_object.angle(self._sire_object, other._sire_object)

        # Return as a BioSimSpace Angle type.
        return _Angle(angle.value(), "DEGREES").radians()

    def x(self):
        """
        Return the x component of the vector.

        Returns
        -------

        x : float
            The x component of the vector.
        """
        return self._sire_object.x().value()

    def y(self):
        """
        Return the y component of the vector.

        Returns
        -------

        y : float
            The y component of the vector.
        """
        return self._sire_object.y().value()

    def z(self):
        """
        Return the z component of the vector.

        Returns
        -------

        z : float
            The z component of the vector.
        """
        return self._sire_object.z().value()

    def magnitude(self):
        """
        Return the magnitude of the vector.

        Returns
        -------

        length : float
            The magnitude of the vector.
        """
        return self._sire_object.magnitude()

    def normalise(self):
        """
        Normalise the vector.

        Returns
        -------

        vector : :class: `Vector <BioSimSpace.Types.Vector>`
            The normalised vector.
        """
        return self._from_sire_vector(self._sire_object.normalise())

    @staticmethod
    def _from_sire_vector(vector):
        """
        Create a vector from a Sire.Maths.Vector object.

        Parameters
        ----------

        vector : Sire.Maths.Vector
            The Sire Vector object.

        Returns
        -------

        vector : :class: `Vector <BioSimSpace.Types.Vector>`
            A BioSimSpace Vector object.
        """
        return Vector(vector.x().value(), vector.y().value(), vector.z().value())
