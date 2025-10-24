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

"""Functionality for generating configuration files for molecular dynamics engines."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Config"]

import math as _math
import warnings as _warnings

from ..Protocol._protocol import Protocol as _ProtocolBase
from .._SireWrappers import System as _System

from .. import Protocol as _Protocol


class Config:
    """Base class for generating configuration files for molecular dynamics engines."""

    def __init__(self, system, protocol, property_map={}):
        """
        Constructor.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.

        protocol : :class:`Protocol <BioSimSpace.Protocol>`
            The protocol for the process.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Don't allow user to create an instance of this base class.
        if type(self) is Config:
            raise Exception("<Config> must be subclassed.")

        # Validate the input.

        if not isinstance(system, _System):
            raise TypeError(
                "'system' must be of type 'BioSimSpace._SireWrappers.System'"
            )

        if not isinstance(protocol, _ProtocolBase):
            raise TypeError("'protocol' must be of type 'BioSimSpace.Protocol'")

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        # Store the attributes.
        self._system = system
        self._protocol = protocol
        self._property_map = property_map

    @staticmethod
    def hasBox(system, property_map={}):
        """
        Whether the system has a box.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        has_box : bool
            Whether the system has a simulation box.
        """

        if not isinstance(system, _System):
            raise TypeError(
                "'system' must be of type 'BioSimSpace._SireWrappers.System'"
            )

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")

        space_prop = property_map.get("space", "space")
        if space_prop in system._sire_object.property_keys():
            try:
                # Make sure that we have a periodic box. The system will now have
                # a default cartesian space.
                box = system._sire_object.property(space_prop)
                has_box = box.is_periodic()
            except:
                has_box = False
        else:
            has_box = False
            _warnings.warn("No simulation box found. Assuming gas phase simulation.")

        return has_box

    @staticmethod
    def hasWater(system):
        """
        Whether the system is contains water molecules.

        Parameters

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.

        Returns
        -------

        has_water : bool
            Whether the system contains water molecules.
        """
        if not isinstance(system, _System):
            raise TypeError(
                "'system' must be of type 'BioSimSpace._SireWrappers.System'"
            )

        return system.nWaterMolecules() > 0

    def reportInterval(self):
        """
        Return the report interval based on the protocol value.

        Returns
        -------

        report_interval : int
            The report interval in integration steps.
        """
        if isinstance(self._protocol, _Protocol.Minimisation):
            return 100
        else:
            report_interval = self._protocol.getReportInterval()
            if report_interval > self.steps():
                report_interval = self.steps()
        return report_interval

    def isRestart(self):
        """
        Return whether this is a restart simulation.

        Returns
        -------

        is_restart : bool
            Whether this is a restart simulation.
        """
        try:
            return self._protocol.isRestart()
        except:
            return False

    def restartInterval(self):
        """
        Return the restart interval based on the protocol value.

        Returns
        -------

        restart_interval : int
            The restart interval in integration steps.
        """
        if isinstance(self._protocol, _Protocol.Minimisation):
            return None
        else:
            restart_interval = self._protocol.getRestartInterval()
            if restart_interval > self.steps():
                restart_interval = self.steps()
        return restart_interval

    def steps(self):
        """
        Return the number of integration steps based on the protocol value.

        Returns
        -------

        steps : int
            The number of integration steps.
        """
        if isinstance(self._protocol, _Protocol.Minimisation):
            steps = self._protocol.getSteps()
        else:
            steps = _math.ceil(
                round(self._protocol.getRunTime() / self._protocol.getTimeStep(), 6)
            )
        return steps
