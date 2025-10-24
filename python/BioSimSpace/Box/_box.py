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

"""Functionality for generating box parameters."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = [
    "boxTypes",
    "generateBoxParameters",
    "cubic",
    "rhombicDodecahedronSquare",
    "rhombicDodecahedronHexagon",
    "truncatedOctahedron",
]

from sire.legacy.Maths import Vector as _Vector
from sire.legacy.Vol import TriclinicBox as _TriclinicBox

from ..Types import Angle as _Angle
from ..Types import Length as _Length


def generateBoxParameters(box_type, image_distance):
    """
    Generate parameters for the named box type with specified image distance.

    Parameters
    ----------

    box_type : str
        The name of the box type. Run BioSimSpace.Box.boxTypes() to get a
        list of the supported boxes.

    image_distance : :class:`Length <BioSimSpace.Types.Length>`
        The image distance.

    Returns
    -------

    box : [:class:`Length <BioSimSpace.Types.Length>`]
        The box vector magnitudes.

    angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
        The box vector angles: yz, xz, and xy.
    """

    if not isinstance(box_type, str):
        raise TypeError("'box_type' must be of type 'str'")
    else:
        # Strip whitespace and convert to lower case.
        box_type = box_type.replace(" ", "").lower()

        if box_type not in _box_types_lower:
            raise ValueError("Supported box types are: %s" % boxTypes())

    return _box_types_dict[box_type](image_distance)


def cubic(image_distance):
    """
    Generate parameters for a cubic box.

    Parameters
    ----------

    image_distance : :class:`Length <BioSimSpace.Types.Length>`
        The image distance.

    Returns
    -------

    box : [:class:`Length <BioSimSpace.Types.Length>`]
        The box vector magnitudes.

    angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
        The box vector angles: yz, xz, and xy.
    """

    # Validate arguments.

    if not isinstance(image_distance, _Length):
        raise TypeError("'image_distance' must be of type 'BioSimSpace.Types.Length'.")

    if image_distance.value() <= 0:
        raise ValueError("'image_distance' must be greater than zero.")

    box = 3 * [image_distance]
    angles = 3 * [_Angle(90, "degrees")]

    return box, angles


def rhombicDodecahedronSquare(image_distance):
    """
    Generate parameters for a square rhombic dodecahedron.

    Parameters
    ----------

    image_distance : :class:`Length <BioSimSpace.Types.Length>`
        The image distance.

    Returns
    -------

    box : [:class:`Length <BioSimSpace.Types.Length>`]
        The box vector magnitudes.

    angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
        The box vector angles: yz, xz, and xy.
    """

    # Validate arguments.

    if not isinstance(image_distance, _Length):
        raise TypeError("'image_distance' must be of type 'BioSimSpace.Types.Length'.")

    if image_distance.value() <= 0:
        raise ValueError("'image_distance' must be greater than zero.")

    # Create the triclinic box.

    triclinic_box = _TriclinicBox.rhombic_dodecahedron_square(
        image_distance.angstroms().value()
    )

    return _get_box_parameters(triclinic_box)


def rhombicDodecahedronHexagon(image_distance):
    """
    Generate parameters for a hexagonal rhombic dodecahedron.

    Parameters
    ----------

    image_distance : :class:`Length <BioSimSpace.Types.Length>`
        The image distance.

    Returns
    -------

    box : [:class:`Length <BioSimSpace.Types.Length>`]
        The box vector magnitudes.

    angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
        The box vector angles: yz, xz, and xy.
    """

    # Validate arguments.

    if not isinstance(image_distance, _Length):
        raise TypeError("'image_distance' must be of type 'BioSimSpace.Types.Length'.")

    if image_distance.value() <= 0:
        raise ValueError("'image_distance' must be greater than zero.")

    # Create the triclinic box.

    triclinic_box = _TriclinicBox.rhombic_dodecahedron_hexagon(
        image_distance.angstroms().value()
    )

    return _get_box_parameters(triclinic_box)


def truncatedOctahedron(image_distance):
    """
    Generate parameters for a truncated octahedron.

    Parameters
    ----------

    image_distance : :class:`Length <BioSimSpace.Types.Length>`
        The image distance.

    Returns
    -------

    box : [:class:`Length <BioSimSpace.Types.Length>`]
        The box vector magnitudes.

    angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
        The box vector angles: yz, xz, and xy.
    """

    # Validate arguments.

    if not isinstance(image_distance, _Length):
        raise TypeError("'image_distance' must be of type 'BioSimSpace.Types.Length'.")

    if image_distance.value() <= 0:
        raise ValueError("'image_distance' must be greater than zero.")

    # Create the triclinic box.

    triclinic_box = _TriclinicBox.truncated_octahedron(
        image_distance.angstroms().value()
    )

    return _get_box_parameters(triclinic_box)


def _get_box_parameters(triclinic_box):
    """
    Internal helper function to get parameters for the passed triclinic box.

    Parameters
    ----------

    triclinic_box : :class `TriclinicBox <Sire.Vol.TriclinicBox>`

    Returns
    -------

    box : [:class:`Length <BioSimSpace.Types.Length>`]
        The box vector magnitudes.
    """

    box = [
        _Length(triclinic_box.vector0().magnitude(), "angstrom"),
        _Length(triclinic_box.vector1().magnitude(), "angstrom"),
        _Length(triclinic_box.vector2().magnitude(), "angstrom"),
    ]

    angles = [
        _Angle(
            _Vector.angle(triclinic_box.vector1(), triclinic_box.vector2()).value(),
            "radians",
        ).degrees(),
        _Angle(
            _Vector.angle(triclinic_box.vector0(), triclinic_box.vector2()).value(),
            "radians",
        ).degrees(),
        _Angle(
            _Vector.angle(triclinic_box.vector0(), triclinic_box.vector1()).value(),
            "radians",
        ).degrees(),
    ]

    return box, angles


# Create a list of the box type names.
# This needs to come after all of the box functions.
_box_types = []  # List of box types (actual names).
_box_types_lower = []  # List of lower case names.
_box_types_dict = {}  # Mapping between lower case names and functions.
import sys as _sys

_namespace = _sys.modules[__name__]
for _var in dir():
    if _var[0] != "_" and _var[0].upper() != "G":
        _box_types.append(_var)
        _box_types_lower.append(_var.lower())
        _box_types_dict[_var.lower()] = getattr(_namespace, _var)
del _namespace
del _sys
del _var


def boxTypes():
    """
    Return a list of the supported box types.

    Returns
    -------

    box_types_fields : [str]
        A list of the supported box types.
    """
    return _box_types
