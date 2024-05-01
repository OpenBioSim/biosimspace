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

"""Custom exceptions for error handling."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = [
    "AlignmentError",
    "AnalysisError",
    "ConversionError",
    "IncompatibleError",
    "MissingSoftwareError",
    "ParameterisationError",
    "StreamError",
    "ThirdPartyError",
]


class AlignmentError(Exception):
    """Exception thrown when molecular alignment fails."""

    pass


class AnalysisError(Exception):
    """Exception thrown when analysis on existing simulation data fails."""


class ConversionError(Exception):
    """Exception thrown when molecular conversion fails."""


class IncompatibleError(Exception):
    """Exception thrown when objects are incompatible with each other."""

    pass


class MissingSoftwareError(Exception):
    """Exception thrown when external software dependencies are missing."""

    pass


class ParameterisationError(Exception):
    """Exception thrown when molecular parameterisation fails."""

    pass


class StreamError(Exception):
    """Exception thrown when streaming fails."""

    pass


class ThirdPartyError(Exception):
    """Exception thrown by a third party package."""
