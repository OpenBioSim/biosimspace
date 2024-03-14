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

"""
Functionality for handling parameterisation protocols.
Author: Lester Hedges <lester.hedges@gmail.com>.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Protocol"]


class Protocol:
    """A base class for parameterisation protocols."""

    def __init__(self, forcefield, ensure_compatible=True, property_map={}):
        """
        Constructor.

        Parameters
        ----------

        forcefield : str
            The name of the force field.

        ensure_compatible : bool
            Whether to ensure that the topology of the parameterised molecule is
            compatible with that of the original molecule. An exception will be
            raised if this isn't the case, e.g. if atoms have been added. Set
            this to False is parameterising lone oxygen atoms corresponding to
            structural (crystal) water molecules.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Don't allow user to create an instance of this base class.
        if type(self) is Protocol:
            raise Exception("<Protocol> must be subclassed.")

        # Validate and set the force field name.
        if not isinstance(forcefield, str):
            raise TypeError("'forcefield' must be of type 'str'")
        else:
            self._forcefield = forcefield

        # Validate and set the topology compatibilty flag.
        if not isinstance(ensure_compatible, bool):
            raise TypeError("'ensure_compatible' must be of type 'bool'")
        else:
            self._ensure_compatible = ensure_compatible

        # Validate and set the property map.
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")
        else:
            self._property_map = property_map.copy()
