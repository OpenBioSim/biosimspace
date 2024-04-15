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

"""Functionality for configuring grids for metadynamics simulation."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Grid"]

from ..Types._type import Type as _Type


class Grid:
    def __init__(self, minimum, maximum, num_bins=None):
        """
        Constructor.

        Define a grid for use in metadynamics simulation.

        Parameters
        ----------

        mimumum : int, float, :class:`Type <BioSimSpace.Types>`
            The minimum value of the grid. Use 'int' or 'float' for
            dimensionless collective variables.

        maximum : int, float, :class:`Type <BioSimSpace.Types>`
            The maximum value of the grid. Use 'int' or 'float' for
            dimensionless collective variables.

        num_bins : int
            The number of bins in the grid. If None, then the number will
            be automatically generated from the metadynamics hill width.
        """

        if num_bins is not None:
            try:
                num_bins = int(num_bins)
            except:
                raise TypeError("'num_bins' must be of type 'int'")

        self.setMinimum(minimum)
        self.setMaximum(maximum)
        if num_bins is not None:
            self.setBins(num_bins)
        else:
            self._num_bins = None

    def __str__(self):
        """Return a human readable string representation of the object."""
        return (
            "<BioSimSpace.Metadynamics.Grid: minimum=%s, maximum=%s, num_bins=%s>"
            % (self._minimum, self._maximum, self._num_bins)
        )

    def __repr__(self):
        """Return a human readable string representation of the object."""
        return self.__str__()

    def __eq__(self, other):
        """Equality operator."""
        return (
            self._minimum == other._minimum
            and self._maximum == other._maximum
            and self._num_bins == other._num_bins
        )

    def setMinimum(self, minimum):
        """
        Set the minimum value of the grid.

        Parameters
        ----------

        minimum : int, float, :class:`Type <BioSimSpace.Types>`
            The minimum value of the grid.
        """
        if not isinstance(minimum, (float, _Type)) and not type(minimum) is int:
            raise TypeError(
                "'minimum' must be of type 'int', 'float', or 'BioSimSpace.Types._type.Type'"
            )
        self._minimum = minimum

        if hasattr(self, "_maximum"):
            if type(self._minimum) is not type(self._maximum):
                raise TypeError("'minimum' and 'maximum' must be of the same type.")
            if self._minimum > self._maximum:
                raise ValueError("'minimum' must be less than 'maximum'")

    def getMinimum(self):
        """
        Get the minimum value of the grid.

        Returns
        -------

        minimum : int, float, :class:`Type <BioSimSpace.Types>`
            The minimum value of the grid.
        """
        return self._minimum

    def setMaximum(self, maximum):
        """
        Set the maximum value of the grid.

        Parameters
        ----------

        maximum : int, float, :class:`Type <BioSimSpace.Types>`
            The maximum value of the grid.
        """
        if not isinstance(maximum, (float, _Type)) and not type(maximum) is int:
            raise TypeError(
                "'maximum' must be of type 'int', 'float', or 'BioSimSpace.Types._type.Type'"
            )
        self._maximum = maximum

        if hasattr(self, "_minimum"):
            if type(self._minimum) is not type(self._maximum):
                raise TypeError("'minimum' and 'maximum' must be of the same type.")
            if self._minimum > self._maximum:
                raise ValueError("'minimum' must be less than 'maximum'")

    def getMaximum(self):
        """
        Get the maximum value of the grid.

        Returns
        -------

        maximum : int, float, :class:`Type <BioSimSpace.Types>`
            The maximum value of the grid.
        """
        return self._maximum

    def setBins(self, num_bins):
        """
        Set the number of bins in the grid.

        Parameters
        ----------

        num_bins : int
            The number of bins in the grid.
        """
        try:
            num_bins = int(num_bins)
        except:
            raise TypeError("'num_bins' must be of type 'int'")

        if num_bins < 1:
            raise ValueError("'num_bins' must be a positive integer.")

        self._num_bins = num_bins

    def getBins(self):
        """
        Get the number of bins in the grid.

        Returns
        -------

        num_bins : int
            The number of bins in the grid.
        """
        return self._num_bins
