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

# The set of supported base units. General unit based types can be created
# by combining these, e.g. via multiplication or division.

__all__ = ["_base_units", "_base_dimensions", "_sire_units_locals"]

import sys as _sys

from ._angle import *
from ._area import *
from ._charge import *
from ._energy import *
from ._length import *
from ._pressure import *
from ._temperature import *
from ._time import *
from ._volume import *

_namespace = _sys.modules[__name__]

# Create the list of base unit types.
_base_units = [
    getattr(_namespace, var) for var in dir() if var[0] != "_" and isinstance(var, type)
]

_base_dimensions = {}
for unit in _base_units:
    _base_dimensions[unit._dimensions] = unit

# Create a local namespace dictionary for the supported Sire units. This
# maps between the unit name and the Sire unit object, allowing us to
# use eval to parse arbitrary expressions based in these units from the
# command-line.
import sire.legacy.Units as _SireUnits

_sire_units_locals = {}

for unit in _base_units:
    for sire_unit in unit._sire_units:
        _sire_units_locals[sire_unit] = getattr(_SireUnits, sire_unit)
