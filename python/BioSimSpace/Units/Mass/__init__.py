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

"""Mass units."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = [
    "kilogram",
    "gram",
    "milligram",
    "microgram",
    "nanogram",
    "picogram",
    "femtogram",
]


from ...Types import Mass as _Mass

kilogram = _Mass(1, "kilogram")
gram = _Mass(1, "gram")
milligram = _Mass(1, "milligram")
microgram = _Mass(1, "microgram")
nanogram = _Mass(1, "nanogram")
picogram = _Mass(1, "picogram")
femtogram = _Mass(1, "femtogram")
