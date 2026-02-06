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
.. currentmodule:: BioSimSpace.Units

Length units
============

.. autosummary::
    :toctree: generated/

    Length.meter
    Length.angstrom
    Length.nanometer
    Length.picometer

Area units
==========

.. autosummary::
    :toctree: generated/

    Area.meter2
    Area.angstrom2
    Area.nanometer2
    Area.picometer2

Angle units
===========

.. autosummary::
    :toctree: generated/

    Angle.radian
    Angle.degree

Volume units
============

.. autosummary::
    :toctree: generated/

    Volume.meter3
    Volume.angstrom3
    Volume.nanometer3
    Volume.picometer3

Charge units
============

.. autosummary::
    :toctree: generated/

    Charge.electron_charge
    Charge.coulomb

Mass units
==========

.. autosummary::
    :toctree: generated/

    Mass.kilogram
    Mass.gram
    Mass.milligram
    Mass.microgram
    Mass.nanogram
    Mass.picogram
    Mass.femtogram

Energy units
============

.. autosummary::
    :toctree: generated/

    Energy.kcal_per_mol
    Energy.kj_per_mol
    Energy.kt

Pressure units
==============

.. autosummary::
    :toctree: generated/

    Pressure.atm
    Pressure.bar

Temperature units
=================

.. autosummary::
    :toctree: generated/

    Temperature.kelvin
    Temperature.celsius
    Temperature.fahrenheit

Time units
==========

.. autosummary::
    :toctree: generated/

    Time.day
    Time.hour
    Time.minute
    Time.second
    Time.millisecond
    Time.nanosecond
    Time.picosecond
    Time.femtosecond
"""

import sire as _sr

_sr.use_new_api()
del _sr

from . import Angle
from . import Area
from . import Charge
from . import Energy
from . import Length
from . import Mass
from . import Pressure
from . import Temperature
from . import Time
from . import Volume

# Whether to allow operations between offset units, see here for details:
# http://pint.readthedocs.io/en/latest/nonmult.html
allow_offset = True
