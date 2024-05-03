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

"""Time units."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = [
    "day",
    "hour",
    "minute",
    "second",
    "millisecond",
    "nanosecond",
    "picosecond",
    "femtosecond",
]

from ...Types import Time as _Time

day = _Time(1, "day")
hour = _Time(1, "hour")
minute = _Time(1, "minute")
second = _Time(1, "second")
millisecond = _Time(1, "millisecond")
nanosecond = _Time(1, "nanosecond")
picosecond = _Time(1, "picosecond")
femtosecond = _Time(1, "femtosecond")
