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

"""Functionality for distance based collective variables."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Distance"]

from math import ceil as _ceil

from ._collective_variable import CollectiveVariable as _CollectiveVariable
from .._bound import Bound as _Bound
from .._grid import Grid as _Grid
from ...Types import Coordinate as _Coordinate
from ...Types import Length as _Length


class Distance(_CollectiveVariable):
    """A class for distance based collective variables."""

    def __init__(
        self,
        atom0,
        atom1,
        hill_width=_Length(0.1, "nanometer"),
        weights0=None,
        weights1=None,
        is_com0=None,
        is_com1=None,
        lower_bound=None,
        upper_bound=None,
        grid=None,
        component=None,
        pbc=True,
    ):
        """
        Constructor.

        Parameters
        ----------

        atom0 : int, [int, int, ...], :class:`Coordinate <BioSimSpace.Types.Coordinate>`
            The atom, group of atoms, or coordinate, that the distance
            will be measured from.

        atom1 : int, [int, int, ...], :class:`Coordinate <BioSimSpace.Types.Coordinate>`
            The atom, group of atoms, or coordinate, that the distance
            will be measured to.

        hill_width : :class:`Length <BioSimSpace.Types.Length>`
            The width of the Gaussian hill used to sample this variable.

        weights0 : [float]
            A list of weights to be used when computing the center of the
            first atom group. This is ignored when a single index is passed
            for 'atom0'.

        weights1 : [float]
            A list of weights to be used when computing the center of the
            second atom group. This is ignored when a single index is passed
            for 'atom1'.

        is_com0 : bool
            Whether to compute the center of mass of the first atom group.
            If True, this option will take precedence over any weights
            passed in via 'weights0'.

        is_com1 : bool
            Whether to compute the center of mass of the second atom group.
            If True, this option will take precedence over any weights
            passed in via 'weights1'.

        lower_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
            A lower bound on the value of the collective variable.

        upper_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
            An upper bound on the value of the collective variable.

        grid : :class:`Grid <BioSimSpace.Metadynamics.Grid>`
            The grid on which the collective variable will be sampled.
            This can help speed up long metadynamics simulations where
            the number of Gaussian kernels can become prohibitive.

        component : str
            Whether to use the 'x', 'y', or 'z' component of the distance
            as the collective variable. If None, then the full Euclidean
            distance is used.

        pbc : bool
            Whether to use periodic boundary conditions when computing the
            collective variable.
        """

        # Call the base class constructor.
        super().__init__()

        # Set the types associated with this collective variable.
        self._types = [_Length]

        # Initialise member data.
        self._atom0 = None
        self._atom1 = None
        self._weights0 = None
        self._weights1 = None
        self._is_com0 = None
        self._is_com1 = None
        self._lower_bound = None
        self._upper_bound = None
        self._grid = None
        self._component = None

        # Set the required parameters.
        self.setAtom0(atom0)
        self.setAtom1(atom1)
        self.setHillWidth(hill_width)
        self.setPeriodicBoundaries(pbc)

        # Set the optional parameters.
        if weights0 is not None:
            self.setWeights0(weights0)
        if weights1 is not None:
            self.setWeights1(weights1)
        if is_com0 is not None:
            self.setCoM0(is_com0)
        if is_com1 is not None:
            self.setCoM1(is_com1)
        if lower_bound is not None:
            self.setLowerBound(lower_bound)
        if upper_bound is not None:
            self.setUpperBound(upper_bound)
        if grid is not None:
            self.setGrid(grid)
        if component is not None:
            self.setComponent(component)

        # Validate that the state is self-consistent.
        self._validate()

        # Flag that the object has been instantiated, i.e. it is no longer "new".
        self._is_new_object = False

    def __str__(self):
        """Return a human readable string representation of the object."""
        string = "<BioSimSpace.Metadynamics.CollectiveVariable.Distance: "
        string += "atom0=%s" % self._atom0
        string += ", atom1=%s" % self._atom1
        string += ", hill_width=%s" % self._hill_width
        if self._weights0 is not None:
            string += ", weights0=%s" % self._weights0
        if self._weights1 is not None:
            string += ", weights1=%s" % self._weights1
        if self._is_com0 is not None:
            string += ", is_com0=%s" % self._is_com0
        if self._is_com1 is not None:
            string += ", is_com1=%s" % self._is_com1
        if self._lower_bound is not None:
            string += ", lower_bound=%s" % self._lower_bound
        if self._upper_bound is not None:
            string += ", upper_bound=%s" % self._upper_bound
        if self._grid is not None:
            string += ", grid=%s" % self._grid
        if self._component is not None:
            string += ", component=%r" % self._component
        string += ", pbc=%s" % self._pbc
        string += ">"
        return string

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return self.__str__()

    def __eq__(self, other):
        """Equality operator."""
        return (
            self._atom0 == other._atom0
            and self._atom1 == other._atom1
            and self._weights0 == other._weights0
            and self._weights1 == other._weights1
            and self._is_com0 == other._is_com0
            and self._is_com1 == other._is_com1
            and self._lower_bound == other._lower_bound
            and self._upper_bound == other._upper_bound
            and self._grid == other._grid
            and self._component == other._component
            and self._pbc == other._pbc
        )

    def setAtom0(self, atom0):
        """
        Set the atom, atoms, or coordinate from which the distance will be
        measured.

        Parameters
        ----------

        atom0 : int, [int, int, ...], [:class:`Length <BioSimSpace.Types.Length>`]
            The atom, group of atoms, or coordinate, that the distance
            will be measured from.
        """

        # Convert tuples to a list.
        if isinstance(atom0, tuple):
            atom0 = list(atom0)

        # Single atom index.
        if type(atom0) is int:
            pass

        # List of atom indices.
        elif isinstance(atom0, list) and all(type(x) is int for x in atom0):
            pass

        # A coordinate.
        elif isinstance(atom0, _Coordinate):
            pass

        # Invalid type.
        else:
            raise TypeError(
                "'atom0' must be of type 'int', a list of 'int' types, "
                "or a 'BioSimSpace.Types.Coordinate' type."
            )

        # Store the existing value.
        old_value = self._atom0

        # All okay, set the value.
        self._atom0 = atom0

        # If we are modifying an existing object, then check for consistency.
        if not self._is_new_object:
            try:
                self._validate()
            except:
                self._atom0 = old_value
                raise

    def getAtom0(self):
        """
        Return the atom, atoms, or coordinate, that the distance will be
        measured from.

        Returns
        -------

        atom0 : int, [int, int, ...], [:class:`Length <BioSimSpace.Types.Length>`]
            The atom, group of atoms, or coordinate, that the distance
            will be measured from.
        """
        return self._atom0

    def setAtom1(self, atom1):
        """
        Set the atom, atoms, or coordinate to which the distance will be
        measured.

        Parameters
        ----------

        atom1 : int, [int, int, ...], [:class:`Length <BioSimSpace.Types.Length>`]
            The atom, group of atoms, or coordinate, that the distance
            will be measured to.
        """

        # Convert tuples to a list.
        if isinstance(atom1, tuple):
            atom1 = list(atom1)

        # Single atom index.
        if type(atom1) is int:
            pass

        # List of atom indices.
        elif isinstance(atom1, list) and all(type(x) is int for x in atom1):
            pass

        # A coordinate.
        elif isinstance(atom1, _Coordinate):
            pass

        # Invalid type.
        else:
            raise TypeError(
                "'atom1' must be of type 'int', a list of 'int' types, "
                "or a 'BioSimSpace.Types.Coordinate' type."
            )

        # Store the existing value.
        old_value = self._atom1

        # All okay, set the value.
        self._atom1 = atom1

        # If we are modifying an existing object, then check for consistency.
        if not self._is_new_object:
            try:
                self._validate()
            except:
                self._atom1 = old_value
                raise

    def getAtom1(self):
        """
        Return the atom, atoms, or coordinate, that the distance will be
        measured to.

        Returns
        -------

        atom1 : int, [int, int, ...], [:class:`Length <BioSimSpace.Types.Length>`]
            The atom, group of atoms, or coordinate, that the distance
            will be measured to.
        """
        return self._atom1

    def setHillWidth(self, hill_width):
        """
        Set the width of the Gaussian hills used to bias this collective
        variable.

        hill_width : :class:`Length <BioSimSpace.Types.Length>`
            The width of the Gaussian hill.
        """
        if not isinstance(hill_width, _Length):
            raise TypeError("'hill_width' must be of type 'BioSimSpace.Types.Length'")

        if hill_width.value() < 0:
            raise ValueError("'hill_width' must have a value of > 0")

        # Convert to the internal unit.
        self._hill_width = hill_width.nanometers()

    def getHillWidth(self):
        """
        Return the width of the Gaussian hill used to bias this collective
        variable.

        Returns
        -------

        hill_width : :class:`Length <BioSimSpace.Types.Length>`
            The width of the Gaussian hill.
        """
        return self._hill_width

    def setWeights0(self, weights0=None):
        """
        Set the weights to be used when computing the center of the first
        atom group. Can be called with no arguments to clear the weights.

        Parameters
        ----------

        weights0 : [float]
            A list of weights to be used when computing the center of the
            first atom group.
        """

        if weights0 is None:
            self._weights0 = None
            return

        if isinstance(weights0, (list, tuple)):
            weights = []

            # Try converting the weights to floats.
            for w in weights0:
                try:
                    weights.append(float(w))
                except:
                    raise TypeError("'weights0' should be a list of 'float' types.")
        else:
            raise TypeError("'weights0' should be a list of 'float' types.")

        # Store the existing value.
        old_value = self._weights0

        # All okay, set the value.
        self._weights0 = list(weights)

        # If we are modifying an existing object, then check for consistency.
        if not self._is_new_object:
            try:
                self._validate()
            except:
                self._weights0 = old_value
                raise

    def getWeights0(self):
        """
        Get the weights to be used when computing the center of the first
        atom group.

        Returns
        -------

        weights0 : [float]
            A list of weights to be used when computing the center of the
            first atom group.
        """
        if self._weights0 is None:
            return None
        else:
            return self._weights0.copy()

    def setWeights1(self, weights1=None):
        """
        Set the weights to be used when computing the center of the second
        atom group. Can be called with no arguments to clear the weights.

        Parameters
        ----------

        weights1 : [float]
            A list of weights to be used when computing the center of the
            second atom group.
        """

        if weights1 is None:
            self._weights1 = None
            return

        if isinstance(weights1, (list, tuple)):
            weights = []

            # Try converting the weights to floats.
            for w in weights1:
                try:
                    weights.append(float(w))
                except:
                    raise TypeError("'weights1' should be a list of 'float' types.")
        else:
            raise TypeError("'weights1' should be a list of 'float' types.")

        # Store the existing value.
        old_value = self._weights1

        # All okay, set the value.
        self._weights1 = list(weights)

        # If we are modifying an existing object, then check for consistency.
        if not self._is_new_object:
            try:
                self._validate()
            except:
                self._weights1 = old_value
                raise

    def getWeights1(self):
        """
        Get the weights to be used when computing the center of the second
        atom group.

        Returns
        -------

        weights1 : [float]
            A list of weights to be used when computing the center of the
            second atom group.
        """
        if self._weights1 is None:
            return None
        else:
            return self._weights1.copy()

    def setCoM0(self, is_com=None):
        """
        Set whether to compute the center of mass of the first atom group.
        If True, this option will take precedence over any weights that may
        have been set. Can be called with no arguments to clear the data.

        Parameters
        ----------

        is_com : bool
            Whether to compute the center of mass of each atom group.
        """
        if is_com is None:
            self._is_com0 = None
            return

        if not isinstance(is_com, bool):
            raise TypeError("'is_com' must be of type 'bool'")

        # Store the existing value.
        old_value = self._is_com0

        self._is_com0 = is_com

        # If we are modifying an existing object, then check for consistency.
        if not self._is_new_object:
            try:
                self._validate()
            except:
                self._is_com0 = old_value
                raise

    def getCoM0(self):
        """
        Whether to compute the center of mass of the first atom group. If
        True, this option will take precedence over any weights that may
        have been set.

        Returns
        -------

        is_com0 : bool
            Whether to compute the center of mass of the first atom group.
        """
        return self._is_com0

    def setCoM1(self, is_com=None):
        """
        Set whether to compute the center of mass of the second atom group.
        If True, this option will take precedence over any weights that may
        have been set. Can be called with no arguments to clear the data.

        Parameters
        ----------

        is_com : bool
            Whether to compute the center of mass of each atom group.
        """
        if is_com is None:
            self._is_com0 = None
            return

        if not isinstance(is_com, bool):
            raise TypeError("'is_com' must be of type 'bool'")

        # Store the existing value.
        old_value = self._is_com1

        self._is_com1 = is_com

        # If we are modifying an existing object, then check for consistency.
        if not self._is_new_object:
            try:
                self._validate()
            except:
                self._is_com1 = old_value
                raise

    def getCoM1(self):
        """
        Whether to compute the center of mass of the first second group. If
        True, this option will take precedence over any weights that may
        have been set.

        Returns
        -------

        is_com1 : bool
            Whether to compute the center of mass of the first atom group.
        """
        return self._is_com1

    def setComponent(self, component=None):
        """
        Whether to use the 'x', 'y', or 'z' component of the distance
        as the collective variable. If unset, then the full Euclidean
        distance is used. Can be called with no argument to clear the
        data.

        Parameters
        ----------

        component : str
            'x', 'y', or 'z'
        """
        if component is None:
            self._component = None
            return

        if not isinstance(component, str):
            raise TypeError("'component' must be of type 'str'")

        allowed = ["x", "y", "z"]

        # Strip whitespace and convert to lower case.
        component = component.replace(" ", "").lower()

        if component not in allowed:
            raise ValueError("'component' should either be 'x', 'y', or 'z'")

        self._component = component

    def getComponent(self):
        """
        Whether to use the 'x', 'y', or 'z' component of the distance
        as the collective variable. If unset, then the full Euclidean
        distance is used.

        Returns
        -------

        component : str
            'x', 'y', or 'z'
        """
        return self._component

    def setPeriodicBoundaries(self, pbc):
        """
        Set whether to use periodic_boundaries when calculating the
        collective variable.

        Parameters
        ----------

        pbc : bool
            Whether to use periodic boundaries conditions.
        """
        if not isinstance(pbc, bool):
            raise TypeError("'pbc' must be of type 'bool'")
        self._pbc = pbc

    def getPeriodicBoundaries(self):
        """
        Return whether to take account of periodic boundary conditions
        when computing the collective variable.

        Returns
        -------

        pbc : bool
            Whether to use periodic boundaries conditions.
        """
        return self._pbc

    def _validate(self):
        """Internal function to check that the object is in a consistent state."""

        if self._weights0 is not None:
            if not isinstance(self._atom0, list):
                raise ValueError(
                    "'weights0' only valid when 'atom0' is a " "list of atom indices."
                )
            elif len(self._weights0) != len(self._atom0):
                raise ValueError(
                    "'weights0' not consistent with 'atom0': "
                    "len(weights0) = %d, len(atom0) = %d"
                    % (len(self._weights0), len(self._atom0))
                )

        if self._weights1 is not None:
            if not isinstance(self._atom1, list):
                raise ValueError(
                    "'weights1' only valid when 'atom1' is a " "list of atom indices."
                )
            elif len(self._weights1) != len(self._atom1):
                raise ValueError(
                    "'weights1' not consistent with 'atom1': "
                    "len(weights1) = %d, len(atom1) = %d"
                    % (len(self._weights1), len(self._atom1))
                )

        if self._is_com0 == True:
            if not isinstance(self._atom0, list):
                raise ValueError(
                    "'is_com0=True but atom0 is not a list of indices. "
                    "Cannot compute center without atom group!"
                )
        if self._is_com1 == True:
            if not isinstance(self._atom1, list):
                raise ValueError(
                    "'is_com1=True but atom1 is not a list of indices. "
                    "Cannot compute center without atom group!"
                )

        if self._lower_bound is not None:
            if not isinstance(self._lower_bound.getValue(), _Length):
                raise TypeError(
                    "'lower_bound' must be of type 'BioSimSpace.Types.Length'"
                )
            # Convert to default unit.
            self._lower_bound.setValue(self._lower_bound.getValue().nanometers())
        if self._upper_bound is not None:
            if not isinstance(self._upper_bound.getValue(), _Length):
                raise TypeError(
                    "'upper_bound' must be of type 'BioSimSpace.Types.Length'"
                )
            # Convert to default unit.
            self._upper_bound.setValue(self._upper_bound.getValue().nanometers())
        if self._lower_bound is not None and self._upper_bound is not None:
            if self._lower_bound.getValue() >= self._upper_bound.getValue():
                raise TypeError("'lower_bound' must less than 'upper_bound'")

        if self._grid is not None:
            if not isinstance(self._grid.getMinimum(), _Length):
                raise TypeError(
                    "'grid' minimum must be of type 'BioSimSpace.Types.Length'"
                )
            # Convert to default unit.
            self._grid.setMinimum(self._grid.getMinimum().nanometers())
            if not isinstance(self._grid.getMaximum(), _Length):
                raise TypeError(
                    "Grid 'maximum' must be of type 'BioSimSpace.Types.Length'"
                )
            # Convert to default unit.
            self._grid.setMaximum(self._grid.getMaximum().nanometers())
            if (
                self._lower_bound is not None
                and self._grid.getMinimum() > self._lower_bound.getValue()
            ):
                raise ValueError("'lower_bound' is less than 'grid' minimum.")
            if (
                self._upper_bound is not None
                and self._grid.getMaximum() < self._upper_bound.getValue()
            ):
                raise ValueError("'upper_bound' is greater than 'grid' maximum.")

            # If the number of bins isn't specified, estimate it out from the hill width.
            if self._grid.getBins() is None:
                grid_range = (self._grid.getMaximum() - self._grid.getMinimum()).value()
                num_bins = _ceil(5.0 * (grid_range / self._hill_width.value()))
                self._grid.setBins(num_bins)
