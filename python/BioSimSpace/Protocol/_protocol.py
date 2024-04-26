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

"""Functionality for handling simulation protocols."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Protocol"]


class Protocol:
    """A base class for holding simulation protocols."""

    def __init__(self):
        """Constructor."""

        # Don't allow user to create an instance of this base class.
        if type(self) is Protocol:
            raise Exception("<Protocol> must be subclassed.")

        # Flag that the protocol hasn't been customised.
        self._is_customised = False

    def getRestraint(self):
        """
        Return the type of restraint.

        Returns
        -------

        restraint : str, [int]
            The type of restraint, either a keyword or a list of atom indices.
        """
        from ._position_restraint_mixin import _PositionRestraintMixin

        if isinstance(self, _PositionRestraintMixin):
            return self._restraint
        else:
            return None

    def _setCustomised(self, is_customised):
        """
        Internal function to flag whether a protocol has been customised.

        Parameters
        ----------

        is_customised : bool
            Whether the protocol has been customised.
        """
        if not isinstance(is_customised, bool):
            raise TypeError("'is_customised' must be of type 'bool'.")

        self._is_customised = is_customised
