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

"""Functionality for configuring restraints on collective variables."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Restraint"]


class Restraint:
    def __init__(self, value, force_constant=100.0, slope=0.0):
        """
        Constructor.

        Set a restraint on the value of a collective variable.

        The expression for the bias is:

        .. math::

            k/2 (x - a)^2 + m (x - a)

        The default restraint is purely harmonic.

        Parameters
        ----------

        value : int, float, :class:`Type <BioSimSpace.Types>`
            The value of the restraint. Use 'int' or 'float' for dimensionless
            collective variables.

        force_constant : float
            The force constant (k) for the harmonic term of the restraint.

        slope : float
            The slope (m) for the linar term of the restraint.
        """

        self.setValue(value)
        self.setForceConstant(force_constant)
        self.setSlope(slope)

    def __str__(self):
        """Return a human readable string representation of the object."""
        return (
            "<BioSimSpace.Metadynamics.Restraint: value=%s, force_constant=%s, slope=%s>"
            % (self._value, self._force_constant, self._slope)
        )

    def __repr__(self):
        """Return a human readable string representation of the object."""
        return self.__str__()

    def __eq__(self, other):
        """Equality operator."""
        return (
            self._value == other._value
            and self._force_constant == other._force_constant
            and self._slope == other._slope
        )

    def setValue(self, value):
        """
        Set the value of the bound.

        Parameters
        ----------

        value : int, float, :class:`Type <BioSimSpace.Types>`
            The value of the bound.
        """
        from ..Types._type import Type as _Type

        if not isinstance(value, (float, _Type)) and not type(value) is int:
            raise TypeError(
                "'value' must be of type 'int', 'float', or 'BioSimSpace.Types._type.Type'"
            )
        self._value = value

    def getValue(self):
        """
        Get the value of the bound.

        Returns
        -------

        value : int, float, :class:`Type <BioSimSpace.Types>`
            The value of the bound.
        """
        return self._value

    def setForceConstant(self, force_constant):
        """
        Set the force constant (k) for the harmonic term of the restraint.

        Parameters
        ----------

        force_constant : float
            The force constant for the harmonic term of the restraint.
        """
        try:
            self._force_constant = float(force_constant)
        except:
            raise TypeError("'force_constant' must be of type 'float'")

    def getForceConstant(self):
        """
        Get the force constant (k) for the harmonic term of the restraint.

        Returns
        -------

        force_constant : float
            The force constant for the harmonic term of the restraint.
        """
        return self._force_constant

    def setSlope(self, slope):
        """
        Set the slope (m) for the linear term of the restraint.

        Parameters
        ----------

        slope : float
            The slope for the linear term of the restraint.
        """
        try:
            self._slope = float(slope)
        except:
            raise TypeError("'slope' must be of type 'float'")

    def getSlope(self):
        """
        Get the slope (m) for the linear term of the restraint.

        Returns
        -------

        slope : float
            The slope for the linear term of the restraint.
        """
        return self._slope
