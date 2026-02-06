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

"""Utility functions."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["protocols", "createProtocol"]

from ._custom import *
from ._equilibration import *
from ._free_energy_minimisation import *
from ._free_energy_equilibration import *
from ._free_energy_production import *
from ._metadynamics import *
from ._minimisation import *
from ._production import *
from ._steering import *

_protocols = []  # List of protocols (actual name).
_protocols_lower = []  # List of lower case names.
_protocol_dict = {}  # Mapping between upper case name and class.
import sys as _sys

_namespace = _sys.modules[__name__]
for _var in dir():
    if _var[0] != "_" and _var != "Custom":
        _protocols.append(_var)
        _protocols_lower.append(_var.lower())
        _protocol_dict[_var.lower()] = getattr(_namespace, _var)
del _namespace
del _sys
del _var


def protocols():
    """
    Return a list of the supported Molecular Dynamics protocols.

    Returns
    -------

    protocols : [str]
       A list of the supported Molecular Dynamics protocols.
    """
    return _protocols


def createProtocol(protocol):
    """
    Create a default simulation protocol.

    Parameters
    ----------

    protocol : str
        The name of the simulation protocol.

    Returns
    -------

    protocol : :class:`Protocol <BioSimSpace.Protocol>`
        The chosen simulation protocol.
    """

    # Strip whitespace and convert to lower case.
    _protocol = protocol.replace(" ", "").lower()

    if _protocol not in _protocols_lower:
        raise KeyError(
            "Unsupported protocol '%s', supported protocols are %s"
            % (protocol, _protocols)
        )

    return _protocol_dict[_protocol]()
