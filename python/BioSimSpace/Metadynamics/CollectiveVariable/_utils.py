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

__all__ = ["list"]

from ._distance import *
from ._funnel import *
from ._torsion import *

# Create a list of the supported collective variables.
_colvars = []
import sys as _sys

_namespace = _sys.modules[__name__]
for _var in dir():
    if _var[0] != "_" and _var[0] != _var[0].lower():
        _colvars.append(_var)
del _namespace
del _sys
del _var


def list():
    """
    Return a list of the supported collective variables.

    Returns
    -------

    colvars : [str]
       A list of the supported collective variables.
    """
    return _colvars
