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

__all__ = ["Funnel", "makeFunnel", "viewFunnel"]


from ._collective_variable import CollectiveVariable as _CollectiveVariable
from .._bound import Bound as _Bound
from .._grid import Grid as _Grid
from ...Types import Length as _Length


class Funnel(_CollectiveVariable):
    """A class for a funnel collective variable."""

    def __init__(
        self,
        atoms0,
        atoms1,
        hill_width=(_Length(0.025, "nanometer"), _Length(0.05, "nanometer")),
        width=_Length(0.6, "nanometers"),
        buffer=_Length(0.15, "nanometers"),
        steepness=1.5,
        inflection=_Length(2.0, "nanometers"),
        lower_bound=_Bound(_Length(0.5, "nanometers"), force_constant=2000),
        upper_bound=_Bound(_Length(4.0, "nanometers"), force_constant=2000),
        grid=(
            _Grid(_Length(0.0, "nanometers"), _Length(4.5, "nanometers"), num_bins=400),
            _Grid(_Length(0.0, "nanometers"), _Length(0.9, "nanometers"), num_bins=400),
        ),
    ):
        """
        Constructor.

        Parameters
        ----------

        atoms0 : [int, int, ...]
            A list of atom indices that define the origin of the funnel.

        atoms1 : [int, int, ...]
            A list of atom indices that define the inflection point of the funnel.

        width : :class:`Length <BioSimSpace.Types.Length>`
           The funnel "wall width".

        buffer : :class:`Length <BioSimSpace.Types.Length>`
           The funnel "wall buffer".

        steepness : float
            The steepness of the funnel at the inflection point.

        inflection : :class:`Length <BioSimSpace.Types.Length>`
            The inflection point as a value of the projection along the
            funnel axis.

        hill_width : (:class:`Length <BioSimSpace.Types.Length>`, "class:`Length <BioSimSpace.Types.Length>`)
            The width of the Gaussian hill used to sample each component
            of the collective variable.

        lower_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
            A lower bound on the value of the collective variable. This is
            used to constrain the "projection on axis" component of the
            collective variable, i.e. the distance from the funnel axis to
            the origin, along the axis.

        upper_bound : :class:`Bound <BioSimSpace.Metadynamics.Bound>`
            An upper bound on the value of the collective variable. This is
            used to constrain the "projection on axis" component of the
            collective variable, i.e. the distance from the funnel axis to
            the origin, along the axis.

        grid : (:class:`Grid <BioSimSpace.Metadynamics.Grid>`, :class:`Grid <BioSimSpace.Metadynamics.Grid>`)
            The grid on which the collective variable will be sampled.
            This can help speed up long metadynamics simulations where
            the number of Gaussian kernels can become prohibitive. The
            collective variable has two components: the "projection on
            the funnel axis" (the distance from the axis to the origin,
            along the axis) and the orthogonal distance between the ligand
            and the axis. The grid should be passed as a tuple with
            the parameters for each component.
        """

        # Call the base class constructor.
        super().__init__()

        # Set the types associated with this collective variable.
        self._types = [_Length]

        # This collective variable has two components: the projection along
        # the funnel axis and the extension orthogonal to it.
        self._num_components = 2

        # Initialise member data.
        self._atoms0 = None
        self._atoms1 = None
        self._lower_bound = None
        self._upper_bound = None
        self._grid = None

        # Set the required parameters.
        self.setAtoms0(atoms0)
        self.setAtoms1(atoms1)
        self.setWidth(width)
        self.setBuffer(buffer)
        self.setSteepness(steepness)
        self.setInflection(inflection)
        self.setHillWidth(hill_width)

        # Set the optional parameters.
        if lower_bound is not None:
            self.setLowerBound(lower_bound)
        if upper_bound is not None:
            self.setUpperBound(upper_bound)
        if grid is not None:
            self.setGrid(grid)

        # Validate that the state is self-consistent.
        self._validate()

        # Flag that the object has been instantiated, i.e. it is no longer "new".
        self._is_new_object = False

    def __str__(self):
        """Return a human readable string representation of the object."""
        string = "<BioSimSpace.Metadynamics.CollectiveVariable.Funnel: "
        string += "atoms0=%s" % self._atoms0
        string += ", atoms1=%s" % self._atoms1
        string += ", width=%s" % self._width
        string += ", buffer=%s" % self._buffer
        string += ", steepness=%s" % self._steepness
        string += ", inflection=%s" % self._inflection
        string += ", hill_width=%s" % (self._hill_width,)
        if self._lower_bound is not None:
            string += ", lower_bound=%s" % self._lower_bound
        if self._upper_bound is not None:
            string += ", upper_bound=%s" % self._upper_bound
        if self._grid is not None:
            string += ", grid=(%s, %s)" % self._grid
        string += ">"
        return string

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return self.__str__()

    def __eq__(self, other):
        """Equality operator."""
        return (
            self._atoms0 == other._atoms0
            and self._atoms1 == other._atoms1
            and self._width == other._width
            and self._buffer == other._buffer
            and self._steepness == other._steepness
            and self._inflection == other._inflection
            and self._hill_width == other._hill_width
            and self._lower_bound == other._lower_bound
            and self._upper_bound == other._upper_bound
            and self._grid == other._grid
        )

    def setAtoms0(self, atoms0):
        """
        Set the list of atom indices whose center-of-mass define the origin
        of the funnel.

        Parameters
        ----------

        atoms0 : [int, int, ...]
            The list of atom indices whose center-of-mass define the origin
            of the funnel.
        """

        # List/tuple of atom indices.
        if isinstance(atoms0, (list, tuple)) and all(type(x) is int for x in atoms0):
            pass

        # Invalid type.
        else:
            raise TypeError("'atoms0' must be of list of 'int' types.")

        # All okay, set the value.
        self._atoms0 = list(atoms0)

    def getAtoms0(self):
        """
        Return indices of the atoms whose center-of-mass defines the origin
        of the funnel.

        Returns
        -------

        atoms0 : [int]
            The list of atom indices whose center-of-mass define the origin
            of the funnel.
        """
        return self._atoms0

    def setAtoms1(self, atoms1):
        """
        Set the list of atom indices whose center-of-mass define the
        inflection point of the funnel.

        Parameters
        ----------

        atoms1 : [int, int, ...]
            The list of atom indices whose center-of-mass define the
            inflection point of the funnel.
        """

        # List/tuple of atom indices.
        if isinstance(atoms1, (list, tuple)) and all(type(x) is int for x in atoms1):
            pass

        # Invalid type.
        else:
            raise TypeError("'atoms1' must be of list of 'int' types.")

        # All okay, set the value.
        self._atoms1 = list(atoms1)

    def getAtoms1(self):
        """
        Return indices of the atoms whose center-of-mass defines the
        inflection point of the funnel.

        Returns
        -------

        atoms1 : [int]
            The list of atom indices whose center-of-mass define the
            inflection point of the funnel.
        """
        return self._atoms1

    def setWidth(self, width):
        """
        Set the funnel "wall width".

        Parameters
        ----------

        width : :class:`Length <BioSimSpace.Types.Length>`
            The width of the funnel wall.
        """
        if not isinstance(width, _Length):
            raise TypeError("'width' must be of type 'BioSimSpace.Types.Length'")

        if width.value() < 0:
            raise ValueError("'width' must have a value of > 0")

        # Convert to the internal unit.
        self._width = width.nanometers()

    def getWidth(self):
        """
        Return the funnel "wall width".

        Returns
        -------

        width : :class:`Length <BioSimSpace.Types.Length>`
            The funnel "wall width".
        """
        return self._width

    def setBuffer(self, buffer):
        """
        Set the funnel "wall buffer".

        Parameters
        ----------

        buffer : :class:`Length <BioSimSpace.Types.Length>`
            The width of the funnel wall buffer.
        """
        if not isinstance(buffer, _Length):
            raise TypeError("'buffer' must be of type 'BioSimSpace.Types.Length'")

        if buffer.value() < 0:
            raise ValueError("'buffer' must have a value of > 0")

        # Convert to the internal unit.
        self._buffer = buffer.nanometers()

    def getBuffer(self):
        """
        Return the funnel "wall buffer".

        Returns
        -------

        buffer : :class:`Length <BioSimSpace.Types.Length>`
            The funnel "wall buffer".
        """
        return self._buffer

    def setSteepness(self, steepness):
        """
        Set the steepness of the funnel at the inflection point.

        Parameters
        ----------

        steepness : float
            The steepness of the funnel at the inflection point.
        """
        # Convert int to float.
        if type(steepness) is int:
            steepness = float(steepness)

        if not isinstance(steepness, float):
            raise TypeError("'steepness' must be of type 'float'")

        if steepness < 0:
            raise ValueError("'steepness' must be > 0")

        self._steepness = steepness

    def getSteepness(self):
        """
        Return the steepness of the funnel at the inflection point.

        Returns
        -------

        steepness : float
            The steepness of the funnel at the inflection point.
        """
        return self._steepness

    def setInflection(self, inflection):
        """
        Set the inflection point as a value of the projection along the
        funnel axis.

        Parameters
        ----------

        inflection : :class:`Length <BioSimSpace.Types.Length>`
            The inflection point as avalue of the projection along the
            funnel axis.
        """
        if not isinstance(inflection, _Length):
            raise TypeError("'inflection' must be of type 'BioSimSpace.Types.Length'")

        if inflection.value() < 0:
            raise ValueError("'inflection' must have a value of > 0")

        # Convert to the internal unit.
        self._inflection = inflection.nanometers()

    def getInflection(self):
        """
        Return the inflection point as a value of the projection along the
        funnel axis.

        Returns
        -------

        inflection : :class:`Length <BioSimSpace.Types.Length>`
            The inflection point as avalue of the projection along the
            funnel axis.
        """
        return self._inflection

    def setHillWidth(self, hill_width):
        """
        Set the width of the Gaussian hills used to bias this collective
        variable.

        Parameters
        ----------

        hill_width : :class:`Length <BioSimSpace.Types.Length>`
            The width of the Gaussian hill for the two components of the
            collective variable: the distance along the funnel projection
            axis, and the orthogonal extent from the axis.
        """

        # Convert list to tuple.
        if isinstance(hill_width, list):
            hill_width = tuple(hill_width)

        if isinstance(hill_width, tuple):
            if len(hill_width) != 2 or not all(
                isinstance(x, _Length) for x in hill_width
            ):
                raise ValueError(
                    "'hill_width' must be a two-component tuple of of type 'BioSimSpace.Metadynamics.Length'"
                )

        for width in hill_width:
            if width.value() < 0:
                raise ValueError("'hill_width' must have a value of > 0")

        # Convert to the internal unit.
        self._hill_width = tuple(x.nanometers() for x in hill_width)

    def getHillWidth(self):
        """
        Return the width of the Gaussian hill used to bias this collective
        variable.

        Returns
        -------

        hill_width : (:class:`Length <BioSimSpace.Types.Length>`, "class:`Length <BioSimSpace.Types.Length>`)
            The width of the Gaussian hill for each component of the
            collective variable.
        """
        return self._hill_width

    def setGrid(self, grid=None):
        """
        Set a grid on which the collective variable will be sampled.
        Call with no arguments to clear the grid.

        Parameters
        ----------

        grid : (:class:`Grid <BioSimSpace.Metadynamics.Grid>`,: class:`Grid <BioSimSpace.Metadynamics.Grid>`)
            A grid for the two components of the collective variable:
            the distance along the funnel projection axis, and the
            orthogonal extent from the axis.
        """

        if grid is None:
            self._grid = None
            return

        # Convert list to tuple.
        if type(grid) is list:
            grid = tuple(grid)

        if isinstance(grid, tuple):
            if len(grid) != 2 or not all(isinstance(x, _Grid) for x in grid):
                raise ValueError(
                    "'grid' must be a two-component tuple of of type 'BioSimSpace.Metadynamics.Grid'"
                )

        # Store the existing value.
        old_value = self._grid

        # Set the new value.
        self._grid = grid

        # If we are modifying an existing object, then check for consistency.
        if not self._is_new_object:
            try:
                self._validate()
            except:
                self._grid = old_value
                raise

    def getGrid(self):
        """
        Get the grid on which the collective variable is sampled.

        Returns
        -------

        grid : (:class:`Grid <BioSimSpace.Metadynamics.Grid>`,: class:`Grid <BioSimSpace.Metadynamics.Grid>`)
            The grid on which the two-components of the collective variable are sampled.
            A grid for the two components of the collective variable:
            the distance along the funnel projection axis, and the
        """
        return self._grid

    def getCorrection(
        self, proj_min=None, proj_max=None, delta=_Length(0.01, "nanometers")
    ):
        """
        Get the funnel correction term. This is the correction factor for
        free-energy estimates that takes into account the reduction in
        entropy caused by the funnel restraint.

        Parameters
        ----------

        proj_min : :class:`Length <BioSimSpace.Types.Length>`
            The minimum coordinate along the projection axis of the funnel.

        proj_max : :class:`Length <BioSimSpace.Types.Length>`
            The maximum coordinate along the projection axis of the funnel.

        delta : :class:`Length <BioSimSpace.Types.Length>`
            The delta for integrating the volume.

        Returns
        -------

        correction : :class:`Energy <BioSimSpace.Types.Energy>`
            The funnel correction.
        """
        import math as _math
        from ...Types import Energy as _Energy
        from ...Types import Volume as _Volume

        if proj_min is None:
            if proj_max is None:
                proj_max = self._upper_bound.getValue()
            else:
                if not isinstance(proj_max, _Length):
                    raise TypeError(
                        "'proj_max' must be of type 'BioSimSpace.Types.Length'."
                    )

            # Default to the last 5 nanometers.
            proj_min = proj_max - _Length(5, "nanometers")

        else:
            if not isinstance(proj_min, _Length):
                raise TypeError(
                    "'proj_min' must be of type 'BioSimSpace.Types.Length'."
                )

        if proj_max is None:
            # Default to the upper bound.
            proj_max = self._upper_bound.getValue()

        else:
            if not isinstance(proj_max, _Length):
                raise TypeError(
                    "'proj_max' must be of type 'BioSimSpace.Types.Length'."
                )

        if proj_min >= proj_max:
            raise ValueError("'proj_min' must be less than 'proj_max'.")

        if not isinstance(delta, _Length):
            raise TypeError("'delta' must be of type 'BioSimSpace.Types.Length'.")

        if delta >= (proj_max - proj_min):
            raise ValueError("'delta' must be larger than 'proj_max' - 'proj_min'.")

        import numpy as _np

        # Create an array of 0.01 nm spaced points between the lower and upper bounds.
        coords = _np.arange(
            proj_min.nanometers().value(),
            proj_max.nanometers().value(),
            delta.nanometers().value(),
        ).tolist()

        # Get the extent values.
        funnel = [
            self.getExtent(_Length(x, "nanometers")).nanometers().value()
            for x in coords
        ]

        # Now integrate the data.
        result = 0
        delta = delta.nanometers().value()
        for x, y in zip(funnel[:-1], funnel[1:]):
            result += (x**2 + y**2) * delta * 0.5

        # Work out the volume of the unbound area.
        volume = _Volume(_math.pi * result, "nanometers cubed")

        # Estimate the average area of the restraint (in Angstrom squared).
        area = (volume / (proj_max - proj_min)).angstroms2()

        # Compute the correction. (1/1660 A-3 is the standard concentration.)
        correction = _Energy(_math.log((area / 1660).value()), "kt")

        return correction

    def getExtent(self, projection):
        """
        Return a value of the funnel extent for a given distance along
        the projection axis.

        Parameters
        ----------

        proj : :class:`Length <BioSimSpace.Types.Length>`
            The distance along the projection axis.

        Returns
        -------

        extent : :class:`Length <BioSimSpace.Types.Length>`
            The distance along the extent axis.
        """
        import math as _math

        if not isinstance(projection, _Length):
            raise TypeError("'projection' must be of type 'BioSimSpace.Types.Length'.")

        extent = (
            self.getWidth().nanometers().value()
            / (
                1
                + _math.exp(
                    self.getSteepness()
                    * (projection - self.getInflection()).nanometers().value()
                )
            )
            + self.getBuffer().nanometers().value()
        )

        return _Length(extent, "nanometers")

    def _validate(self):
        """Internal function to check that the object is in a consistent state."""
        import math as _math

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
            # Check the two components of the grid.
            if not isinstance(self._grid[0].getMinimum(), _Length):
                raise TypeError(
                    "'grid' minimum must be of type 'BioSimSpace.Types.Length'"
                )
            if not isinstance(self._grid[1].getMinimum(), _Length):
                raise TypeError(
                    "'grid' minimum must be of type 'BioSimSpace.Types.Length'"
                )
            # Convert to default unit.
            self._grid[0].setMinimum(self._grid[0].getMinimum().nanometers())
            self._grid[1].setMinimum(self._grid[1].getMinimum().nanometers())
            if not isinstance(self._grid[0].getMaximum(), _Length):
                raise TypeError(
                    "Grid 'maximum' must be of type 'BioSimSpace.Types.Length'"
                )
            if not isinstance(self._grid[1].getMaximum(), _Length):
                raise TypeError(
                    "Grid 'maximum' must be of type 'BioSimSpace.Types.Length'"
                )
            # Convert to default unit.
            self._grid[0].setMaximum(self._grid[0].getMaximum().nanometers())
            self._grid[1].setMaximum(self._grid[1].getMaximum().nanometers())
            # Lower and upper bounds only apply to the grid parameters for the
            # first component of the collective variable, the projection on the
            # funnel axis.
            if (
                self._lower_bound is not None
                and self._grid[0].getMinimum() > self._lower_bound.getValue()
            ):
                raise ValueError("'lower_bound' is less than 'grid' minimum.")
            if (
                self._upper_bound is not None
                and self._grid[0].getMaximum() < self._upper_bound.getValue()
            ):
                raise ValueError("'upper_bound' is greater than 'grid' maximum.")

            # If the number of bins isn't specified, estimate it out from the hill width.
            if self._grid[0].getBins() is None:
                grid_range = (
                    self._grid[0].getMaximum() - self._grid[0].getMinimum()
                ).value()
                num_bins = _math.ceil(5.0 * (grid_range / self._hill_width.value()))
                self._grid[0].setBins(num_bins)
            if self._grid[1].getBins() is None:
                grid_range = (
                    self._grid[1].getMaximum() - self._grid[1].getMinimum()
                ).value()
                num_bins = _math.ceil(5.0 * (grid_range / self._hill_width.value()))
                self._grid[1].setBins(num_bins)


def makeFunnel(
    system, protein=None, ligand=None, alpha_carbon_name="CA", property_map={}
):
    """
    Calculate the two sets of atom indices, atoms0 and atoms1, that are used
    to define the funnel collective variable.

    Parameters
    ----------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The system of interest. This is assumed to be a solvated
        protein-ligand system.

    protein : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, int
        The protein molecule. This can either be a Molecule object, or
        the index of the protein within the passed system. If None is
        passed then we assume that the protein is the largest molecule
        in the system.

    ligand : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, int, \
            :class:`Coordinate <BioSimSpace.Types.Coordinate>`
        The ligand molecule. This can either be a Molecule object, the
        index of the ligand within the passed system, or the position of
        the binding site within the protein. If None is passed then we
        assume that the ligand is the second largest molecule in the system.

    alpha_carbon_name : str
        The name of alhpa carbon atoms in the system topology.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    atoms0 : [int]
        A list of atom indices that define the origin of the funnel.

    atoms1 : [int]
        A list of atom indices that define the inflection point of the funnel.
    """
    from ..._Exceptions import IncompatibleError as _IncompatibleError
    from ..._SireWrappers import Molecule as _Molecule
    from sire.legacy.Maths import Vector as _SireVector
    from ..._SireWrappers import System as _System
    from ...Types import Coordinate as _Coordinate

    # Validate the input.

    # System.
    if not isinstance(system, _System):
        raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'.")

    # Protein.
    if isinstance(protein, _Molecule):
        # Make sure the molecule exists in the system.
        try:
            system.getIndex(protein)
        except KeyError:
            raise (ValueError("The passed 'protein' was not found in the 'system'."))

        protein = protein.copy()

    elif type(protein) is int:
        try:
            protein = system[protein]
        except:
            ValueError(f"Couldn't find 'protein' index {protein} in the 'system'.")

    elif isinstance(protein, _Coordinate):
        pass

    elif protein is None:
        # Store the number of atoms in each molecule.
        atom_nums = []
        for mol in system:
            atom_nums.append(mol.nAtoms())

        # Sort the molecule indices by the number of atoms they contain.
        sorted_nums = sorted((value, idx) for idx, value in enumerate(atom_nums))

        # Set the protein to the largest molecule.
        protein = system[sorted_nums[-1][1]]

    else:
        raise TypeError(
            "'protein' must be of type 'BioSimSpace._SireWrappers.Molecule' or 'None'."
        )

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'.")

    # Get the "coordinates" property from the user mapping.
    coordinates = property_map.get("coordinates", "coordinates")
    if not protein._sire_object.has_property(coordinates):
        raise ValueError(
            f"The 'protein' molecule doesn't have a {coordinates} property!"
        )

    # Ligand.
    if isinstance(ligand, _Molecule):
        # Make sure the molecule exists in the system.
        try:
            system.getIndex(ligand)
        except KeyError:
            raise (ValueError("The passed 'ligand' was not found in the 'system'."))

    elif type(ligand) is int:
        try:
            ligand = system[ligand]
        except:
            ValueError(f"Couldn't find 'ligand' index {ligand} in the 'system'.")
        binding_site = ligand

    elif isinstance(ligand, _Coordinate):
        binding_site = ligand
        pass

    elif ligand is None:
        # Store the number of atoms in each molecule.
        atom_nums = []
        for mol in system:
            atom_nums.append(mol.nAtoms())

        # Sort the molecule indices by the number of atoms they contain.
        sorted_nums = sorted((value, idx) for idx, value in enumerate(atom_nums))

        # Set the ligand to the second largest molecule.
        ligand = system[sorted_nums[-2][1]]

    else:
        raise TypeError(
            "'ligand' must be of type 'BioSimSpace._SireWrappers.Molecule' or 'None'."
        )

    # Alpha carbon name.
    if not isinstance(alpha_carbon_name, str):
        raise TypeError("'alpha_carbon_name' must be of type 'str'.")

    # Get the "space" property from the user map.
    space_prop = property_map.get("space", "space")
    if space_prop not in system._sire_object.property_keys():
        raise _IncompatibleError("The system contains no simulation box property!")

    # Store the space.
    space = system._sire_object.property(space_prop)

    # The following is adapted from funnel_maker.py by Dominykas Lukauskis.

    # To compute the funnel vector we project a grid on the binding site and
    # determine the average x,y,z coordinate that doesn't contain any protein
    # atoms within a specified search radius. This is more robust than searching
    # for solvent molecules directly, since it allows for dry binding sites.

    # If the ligand is a Molecule, then assume the binding site is the ligand
    # center of mass.
    if isinstance(ligand, _Molecule):
        com = ligand._getCenterOfMass(space=space, property_map=property_map)
        binding_site = _Coordinate(*com)

    # Set up the grid.

    # 20 Angstrom grid edge.
    grid_length = _Length(20, "A")

    # Number of points along each edge.
    num_edge = 5

    # Search radius.
    search_radius = (grid_length / num_edge) / 2

    # Work out the grid extent.
    grid_min = binding_site - 0.5 * grid_length
    grid_max = binding_site + 0.5 * grid_length

    # Initialise a vector to store the average non-protein coordinate.
    non_protein = _SireVector()
    num_non_protein = 0

    # Create a property map to use for the search.
    search_map = {"space": space}

    import numpy as _np

    # Loop over x grid points.
    for x in _np.linspace(
        grid_min.x().angstroms().value(), grid_max.x().angstroms().value(), num_edge
    ):
        # Loop over y grid points.
        for y in _np.linspace(
            grid_min.y().angstroms().value(), grid_max.y().angstroms().value(), num_edge
        ):
            # Loop over z grid points.
            for z in _np.linspace(
                grid_min.z().angstroms().value(),
                grid_max.z().angstroms().value(),
                num_edge,
            ):
                # Generate the search string.
                string = (
                    f"atoms within {search_radius.angstroms().value()} of {x},{y},{z}"
                )

                # Search the protein for atoms with the search radius of the
                # point x,y,z.
                try:
                    search = protein.search(string, property_map=search_map)

                # If there are no protein atoms then add the grid coordinate
                # to our running total.
                except:
                    non_protein += _SireVector(x, y, z)
                    num_non_protein += 1

    # Work out the average non-protein coordinate.
    non_protein /= num_non_protein
    non_protein = _Coordinate._from_sire_vector(non_protein)

    # Now select all alpha carbon atoms within 10 Angstrom of the ligand or grid.

    # Generate the search string.
    x = binding_site.x().angstroms().value()
    y = binding_site.y().angstroms().value()
    z = binding_site.z().angstroms().value()
    string = f"(atoms within 10 of {x},{y},{z}) and atomname {alpha_carbon_name}"

    # Perform the search.
    try:
        search = system.search(string, property_map=search_map)

    # Raise exception if no atoms were found.
    except:
        raise ValueError(
            "No alpha carbon atoms found within 10 Angstrom of "
            "the binding pocket center. Try explicitly setting "
            "the center using a 'BioSimSpace.Types.Coordinate' "
            "or using a different option for 'alpha_carbon_name'."
        )

    # Work out the center of mass of the alpha carbon atoms. (In Angstrom.)
    com = _Coordinate(_Length(0, "A"), _Length(0, "A"), _Length(0, "A"))
    atoms1 = []
    for atom in search:
        com += atom.coordinates(property_map=property_map)
        atoms1.append(system.getIndex(atom))
    com /= search.nResults()

    # Compute the normal vector for the funnel.
    initial_funnel_normal_vector = (non_protein - com).to_vector().normalise()

    # Compute the location of a point 10 Angstom in the direction of the funnel
    # normal vector in the direction of the protein.
    into_the_protein = com.to_vector() - 10 * initial_funnel_normal_vector

    # Search for all alpha carbons within 7 Angstrom of the point.

    # Generate the search string.
    x = into_the_protein.x()
    y = into_the_protein.y()
    z = into_the_protein.z()
    string = f"(atoms within 7 of {x},{y},{z}) and atomname {alpha_carbon_name}"

    # Perform the search.
    try:
        search = system.search(string)

    # Raise exception if no atoms were found.
    except:
        raise ValueError(
            "No alpha carbon atoms found within 10 Angstrom of "
            "the binding pocket center. Try explicitly setting "
            "the center using a 'BioSimSpace.Types.Coordinate' "
            "or using a different option for 'alpha_carbon_name'."
        )

    # Append the indices of these atoms to the atoms0 vector.
    atoms0 = []
    for atom in search:
        atoms0.append(system.getIndex(atom))

    return atoms0, atoms1


def viewFunnel(system, collective_variable, property_map={}):
    """
    Visualise the shape of the funnel defined by the collective variable.
    This is useful for checking that the funnel configuration is sane
    prior to running any metadynamics simulations. The function returns a
    BioSimSpace.Notebook.View object, allowing the user to customise the
    visulisation.

    Parameters
    ----------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The system of interest. This is assumed to be a solvated
        protein-ligand system.

    collective_variable : :class:`Funnel <BioSimSpace.Metadynamics.CollectiveVariable.Funnel>`
        The funnel collective variable.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    view : :class:`View <BioSimSpace.Notebook.View>`
        A view object showing the system and funnel.
    """
    from ..._SireWrappers import Molecule as _Molecule
    from ... import _is_notebook
    from sire.legacy.Maths import Vector as _SireVector
    from ..._SireWrappers import System as _System
    import sire.legacy.Mol as _SireMol

    # The following is adapted from funnel_maker.py by Dominykas Lukauskis.

    # Don't do anything if this isn't called from within a notebook.
    if not _is_notebook:
        return None

    # Validate the input.

    if not isinstance(system, _System):
        raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'.")

    if not isinstance(collective_variable, Funnel):
        raise TypeError(
            "'collective_variable' must be of type "
            "'BioSimSpace.Metadynamics.CollectiveVariable.Funnel'"
        )

    # Store the number of atoms in each molecule.
    atom_nums = []
    for mol in system:
        atom_nums.append(mol.nAtoms())

    # Sort the molecule indices by the number of atoms they contain.
    sorted_nums = sorted((value, idx) for idx, value in enumerate(atom_nums))

    # Set the protein to the largest molecule.
    protein = system[sorted_nums[-1][1]]

    # Get the "coordinates" property from the user mapping.
    coordinates = property_map.get("coordinates", "coordinates")
    if not protein._sire_object.has_property(coordinates):
        raise ValueError(
            f"The 'protein' molecule doesn't have a {coordinates} property!"
        )

    # Get the atoms that define the origin and inflection point of the funnel.
    atoms0 = collective_variable.getAtoms0()
    atoms1 = collective_variable.getAtoms1()

    # Store the protein atoms.
    atoms = protein.getAtoms()

    # Work out the center-of-mass of each set of atoms.
    com0 = _SireVector()
    com1 = _SireVector()

    for atom in atoms0:
        try:
            com0 += atoms[atom]._sire_object.property(coordinates)
        except:
            raise ValueError(f"Could not obtain coordinates for atom index '{atom}'!")
    com0 /= len(atoms0)
    com0 = [x.value() for x in com0]

    for atom in atoms1:
        try:
            com1 += atoms[atom]._sire_object.property(coordinates)
        except:
            raise ValueError(f"Could not obtain coordinates for atom index '{atom}'!")
    com1 /= len(atoms1)
    com1 = [x.value() for x in com1]

    # Create a new molecule to hold the funnel.
    funnel_mol = _SireMol.Molecule("Funnel")

    # Add a single residue called FUN.
    res = funnel_mol.edit().add(_SireMol.ResNum(1))
    res.rename(_SireMol.ResName("FUN"))

    # Create a single cut-group.
    cg = res.molecule().add(_SireMol.CGName("1"))

    # Counter for the number of atoms.
    num = 1

    # Funnel plot variables.
    vec_step = 2
    n_angle_samples = 8

    # Extract collective variable data.
    lower_wall = collective_variable.getLowerBound().getValue().angstroms().value()
    upper_wall = collective_variable.getUpperBound().getValue().angstroms().value()
    wall_width = collective_variable.getWidth().angstroms().value()
    beta_cent = collective_variable.getSteepness()
    s_cent = collective_variable.getInflection().angstroms().value()
    wall_buffer = collective_variable.getBuffer().angstroms().value()

    # Get the element property from the map.
    element = property_map.get("element", "element")

    import numpy as _np

    # Calculate the vector defined by points p0 and p1.
    vec = _np.array(com1, dtype=float) - _np.array(com0, dtype=float)

    # BEWARE: inconsistency with linalg, if vec is a list and not an array!
    # Make it a unit vector.
    unit_vec = vec / _np.linalg.norm(vec)

    # Now to get orthogonal vectors:
    # https://math.stackexchange.com/questions/133177/finding-a-unit-vector-perpendicular-to-another-vector

    # Determine 1st orthogonal vector.
    a0 = _np.random.randint(1, 10)
    a1 = _np.random.randint(1, 10)
    a2 = -(a0 * vec[0] + a1 * vec[1]) / vec[2]
    a = _np.asarray([a0, a1, a2])
    unit_a = a / _np.linalg.norm(a)

    # Determine 2nd orthogonal vector.

    unit_b = _np.cross(unit_a, unit_vec)
    # Iterate along the selected vector
    funnel_coords = []
    for step in _np.arange(lower_wall, upper_wall + vec_step, vec_step):
        # iterate around a circle with its radius defined by the sigmoid function
        radius = (wall_width / (1 + _np.exp(beta_cent * (step - s_cent)))) + wall_buffer
        for angle in _np.arange(-_np.pi, _np.pi, 2 * _np.pi / n_angle_samples):
            # Calculate parametric functions for this specific case
            # https://math.stackexchange.com/questions/73237/parametric-equation-of-a-circle-in-3d-space
            # Generate pseudoatoms along the axis.
            coord = (
                com0
                + unit_vec * step
                + radius * (_np.cos(angle) * unit_a + _np.sin(angle) * unit_b)
            )
            funnel_coords.append(_SireVector(coord))

            # Add the pseudoatom.
            added = cg.add(_SireMol.AtomName("HE"))
            added.renumber(_SireMol.AtomNum(num))
            added.reparent(_SireMol.ResIdx(0))
            num += 1

    # Recreate the funnel molecule and make it editable.
    funnel_mol = cg.molecule().commit().edit()

    # Create a Helium element.
    helium = _SireMol.Element("He")

    # Add the coordinaates and element property.
    for x in range(0, funnel_mol.nAtoms()):
        idx = _SireMol.AtomIdx(x)
        funnel_mol = (
            funnel_mol.atom(idx).set_property(coordinates, funnel_coords[x]).molecule()
        )
        funnel_mol = funnel_mol.atom(idx).set_property(element, helium).molecule()

    # Add the funnel pseudoatoms to the system.
    new_system = system + _Molecule(funnel_mol.commit())

    from ...Notebook import View as _View

    # Create a BioSimSpace notebook View.
    view = _View(new_system)

    return view
