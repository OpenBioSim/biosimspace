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

"""Functionality for configuring bounds on collective variables."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Bound"]

from ..Types._type import Type as _Type


class Bound:
    def __init__(self, value, force_constant=100.0, exponent=2.0, epsilon=1.0):
        """
        Constructor.

        Set a bound on the value of a collective variable along with the
        parameters used to define the bias potential.

        The expression for the bias is:

        .. math::

            k ((x - a)/s)^e

        Parameters
        ----------

        value : int, float, :class:`Type <BioSimSpace.Types>`
            The value of the bound. Use 'int' or 'float' for dimensionless
            collective variables.

        force_constant : float
            The force constant (k) for the bias potential. Note that the
            units of the force constant aren't specified, i.e. the value
            takes the units of the collective variable to which it is later
            applied.

        exponent : float
            The exponent (e) for the bias potential.

        epsilon : float
            The rescaling factor (s) for the bias potential.
        """

        self.setValue(value)
        self.setForceConstant(force_constant)
        self.setExponent(exponent)
        self.setEpsilon(epsilon)

    def __str__(self):
        """Return a human readable string representation of the object."""
        return (
            "<BioSimSpace.Metadynamics.Bound: value=%s, force_constant=%s, exponent=%s, epsilon=%s>"
            % (self._value, self._force_constant, self._exponent, self._epsilon)
        )

    def __repr__(self):
        """Return a human readable string representation of the object."""
        return self.__str__()

    def __eq__(self, other):
        """Equality operator."""
        return (
            self._value == other._value
            and self._force_constant == other._force_constant
            and self._exponent == other._exponent
            and self._epsilon == other._epsilon
        )

    def setValue(self, value):
        """
        Set the value of the bound.

        Parameters
        ----------

        value : int, float, :class:`Type <BioSimSpace.Types>`
            The value of the bound.
        """
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
        Set the force constant (k) for the bias potential.

        Parameters
        ----------

        force_constant : float
            The force constant for the bias potential.
        """
        try:
            self._force_constant = float(force_constant)
        except:
            raise TypeError("'force_constant' must be of type 'float'")

    def getForceConstant(self):
        """
        Get the force constant (k) for the bias potential.

        Returns
        -------

        force_constant : float
            The force constant for the bias potential.
        """
        return self._force_constant

    def setExponent(self, exponent):
        """
        Set the exponent (e) for the bias potential.

        Parameters
        ----------

        exponent : float
            The exponent for the bias potential.
        """
        try:
            self._exponent = float(exponent)
        except:
            raise TypeError("'exponent' must be of type 'float'")

    def getExponent(self):
        """
        Get the exponent (e) for the bias potential.

        Returns
        -------

        exponent : float
            The exponent for the bias potential.
        """
        return self._exponent

    def setEpsilon(self, epsilon):
        """
        Set the rescaling factor (s) for the bias potential.

        Parameters
        ----------

        epsilon : float
            The rescaling factor for the bias potential.
        """
        try:
            self._epsilon = float(epsilon)
        except:
            raise TypeError("'epsilon' must be of type 'float'")

    def getEpsilon(self):
        """
        Get the rescaling factor (s) for the bias potential.

        Returns
        -------

        epsilon : float
            The rescaling factor for the bias potential.
        """
        return self._epsilon
