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

"""
.. currentmodule:: BioSimSpace.Types

Classes
=======

.. autosummary::
    :toctree: generated/

    Angle
    Area
    Charge
    Coordinate
    Energy
    Length
    Mass
    Pressure
    Temperature
    Time
    Vector
    Volume
"""

import sire as _sr

_sr.use_new_api()
del _sr

from ._angle import *
from ._area import *
from ._charge import *
from ._coordinate import *
from ._energy import *
from ._length import *
from ._mass import *
from ._pressure import *
from ._temperature import *
from ._time import *
from ._vector import *
from ._volume import *

# Hide GeneralUnit since it will be automatically created from combinations
# of the unit based types above.
from ._general_unit import GeneralUnit as _GeneralUnit
