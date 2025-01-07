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

"""Functionality for finding and managing hardware resources."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["ResourceManager"]

import argparse as _argparse


class ResourceManager:
    """A class for finding and managing hardware resources."""

    def __init__(self):
        """Constructor."""

        # Set default values.
        self._nodes = None
        self._cpus = None
        self._gpus = None

        # Create the argument parser.
        self._parser = _argparse.ArgumentParser(
            prog="BioSimSpace",
            description="Command-line parser for hardware resources",
            add_help=False,
            allow_abbrev=False,
        )

        # Add the arguments.
        self._parser.add_argument(
            "--nodes", type=int, help="The number of hardware nodes."
        )
        self._parser.add_argument(
            "--cpus", type=int, help="The number of hardware central processing units."
        )
        self._parser.add_argument(
            "--gpus", type=int, help="The number of hardware graphics processors."
        )

    def _initialise(self):
        """Initialise the resource manager."""

        # Parse the arguments into a dictionary.
        args = vars(self._parser.parse_known_args()[0])

        # Loop over all of the arguments and set the values.
        for key, value in args.items():
            if key == "nodes":
                if value is not None:
                    self.setNodes(value)
            elif key == "cpus":
                if value is not None:
                    self.setCPUs(value)
            elif key == "gpus":
                if value is not None:
                    self.setGPUs(value)

    def getNodes(self):
        """
        Return the number of nodes.

        Returns
        -------

        nodes : int
            The number of nodes.
        """
        return self._nodes

    def setNodes(self, nodes):
        """
        Set the number of nodes.

        Parameters
        ----------

        nodes : int
            The number of nodes.
        """

        if not type(nodes) is int:
            raise TypeError("'nodes' must be of type 'int'.")

        if nodes < 0:
            raise ValueError("'nodes' cannot be negative!")

        self._nodes = nodes

    def getCPUs(self):
        """
        Return the number of cpus.

        Returns
        -------

        cpus : int
            The number of CPUs.
        """
        return self._cpus

    def setCPUs(self, cpus):
        """
        Set the number of CPUs.

        Parameters
        ----------

        cpus : int
            The number of cpus.
        """

        if not type(cpus) is int:
            raise TypeError("'cpus' must be of type 'int'.")

        if cpus < 0:
            raise ValueError("'cpus' cannot be negative!")

        self._cpus = cpus

    def getGPUs(self):
        """
        Return the number of GPUs.

        Returns
        -------

        gpus : int
            The number of GPUs.
        """
        return self._gpus

    def setGPUs(self, gpus):
        """
        Set the number of GPUs.

        Parameters
        ----------

        gpus : int
            The number of GPUs.
        """

        if not type(gpus) is int:
            raise TypeError("'gpus' must be of type 'int'.")

        if gpus < 0:
            raise ValueError("'gpus' cannot be negative!")

        self._gpus = gpus
