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

__all__ = ["engines", "createProcess"]

from ._amber import *
from ._gromacs import *
from ._namd import *
from ._openmm import *
from ._process_runner import *
from ._somd import *

_engines = []  # List of supported engines (actual name).
_engines_lower = []  # List of lower case engine names.
_engine_dict = {}  # Mapping between lower case name and class.
import sys as _sys

_namespace = _sys.modules[__name__]
for _var in dir():
    if _var[0] != "_" and _var != "ProcessRunner":
        _engines.append(_var)
        _engines_lower.append(_var.lower())
        _engine_dict[_var.lower()] = getattr(_namespace, _var)
del _namespace
del _sys
del _var


def engines():
    """
    Return a list of the supported Molecular Dynamics engines.

    Returns
    -------

    engines : [str]
        The list of supported Molecular Dynamics engines.
    """
    return _engines


def createProcess(system, protocol, engine, **kwargs):
    """
    Create a simulation process.

    Parameters
    ----------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The molecular system.

    protocol : :class:`Protocol <BioSimSpace.Protocol>`
        The protocol for the process.

    engine : str
        The name of the simulation engine.

    kwargs : dict
        A dictionary of optional keyword arguments neeeded by the engine.

    Returns
    -------

    process : :class:`Process <BioSimSpace.Process>`
        The process object for the specific simulation engine.
    """

    # Strip whitespace and convert to lower case.
    _engine = engine.replace(" ", "").lower()

    if _engine not in _engines_lower:
        raise KeyError(
            "Unsupported engine '%s', supported engines are %s" % (engine, _engines)
        )

    return _engine_dict[_engine](system, protocol, **kwargs)
