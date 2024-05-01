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

"""
Functionality for storing custom simulation protocols.
Author: Lester Hedges <lester.hedges@gmail.com>.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Custom"]

import os as _os

from ._protocol import Protocol as _Protocol


class Custom(_Protocol):
    """A class for storing custom protocols."""

    def __init__(self, config):
        """
        Constructor.

        Parameters
        ----------

        config : str, [ str ]
            The custom protocol configuration.
        """

        # Call the base class constructor.
        super().__init__()

        # Set the protocol configuration.
        self.setConfig(config)

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.Protocol.Custom>"

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.Protocol.Custom>"

    def __eq__(self, other):
        """Equality operator."""
        return self._config == other._config

    def getConfig(self):
        """
        Return the custom configuration.

        Returns
        -------

        config : [ str ]
            Return the list of configuration strings.
        """
        return self._config.copy()

    def setConfig(self, config):
        """
        Set the custom configuration.

        Parameters
        ----------

        config : str, [ str ]
            A config file, or list of configuration strings.
        """
        # Check that the passed configuration is a list of strings.
        if _is_list_of_strings(config):
            self._config = config

        # The user has passed a path to a file.
        elif _os.path.isfile(config):
            # Clear the existing config.
            self._config = []

            # Read the contents of the file.
            with open(config, "r") as file:
                for line in file:
                    self._config.append(line.rstrip())
        else:
            raise ValueError("'config' must be a list of strings, or a file path.")


def _is_list_of_strings(lst):
    """Check whether the passed argument is a list of strings."""
    if lst and isinstance(lst, list):
        return all(isinstance(elem, str) for elem in lst)
    else:
        return False
