######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2023
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

"""Functionality for solvating molecular systems."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["solvate", "spc", "spce", "tip3p", "tip4p", "tip5p", "waterModels"]

import os as _os
import re as _re
import subprocess as _subprocess
import shlex as _shlex
import sys as _sys
import warnings as _warnings

from sire.legacy import Base as _SireBase
from sire.legacy import IO as _SireIO
from sire.legacy import Mol as _SireMol
from sire.legacy.Maths import Vector as _Vector
from sire.legacy.Vol import TriclinicBox as _TriclinicBox

from sire.legacy.Units import degree as _degree

from .. import _gmx_exe, _gmx_path
from .. import _isVerbose

from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import System as _System
from .._SireWrappers import Molecule as _Molecule
from .._SireWrappers import Molecules as _Molecules
from ..Types import Coordinate as _Coordinate
from ..Types import Angle as _Angle
from ..Types import Length as _Length

from .. import IO as _IO
from .. import _Utils


def solvate(
    model,
    molecule=None,
    box=None,
    angles=3 * [_Angle(90, "degrees")],
    shell=None,
    ion_conc=0,
    is_neutral=True,
    is_aligned=False,
    match_water=True,
    work_dir=None,
    property_map={},
):
    """
    Solvate with the specified water model.

    Parameters
    ----------

    model : str
        The name of the water model.

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
               :class:`Molecule <BioSimSpace._SireWrappers.Molecules>`, \
               :class:`System <BioSimSpace._SireWrappers.System>`
        A molecule, or container/system of molecules.

    box : [:class:`Length <BioSimSpace.Types.Length>`]
        A list containing the box size in each dimension: x, y, and z.

    angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
        A list containing the angles between the box vectors: yz, xz, and xy.

    shell : :class:`Length` <BioSimSpace.Types.Length>`
        Thickness of the water shell around the solute. Note that the
        base length of the resulting box must be at least twice as large
        as the cutoff used by the chosen molecular dynamics engine. As such,
        the shell option is often unsuitable for small molecules.

    ion_conc : float
        The ion concentration in (mol per litre).

    is_neutral : bool
        Whether to neutralise the system.

    is_aligned : bool
        Whether to align the principal axes of the molecule to those of the
        solvent box.

    match_water : bool
        Whether to update the naming of existing water molecules to match the
        expected convention for GROMACS, which is used as the solvation engine.

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The solvated molecular system.
    """

    if not isinstance(model, str):
        raise TypeError("'model' must be of type 'str'")
    else:
        # Strip whitespace and convert to lower case.
        model = model.replace(" ", "").lower()

        # Check that this water model is supported.
        if model not in _models_lower:
            raise ValueError("Supported water models are: %s" % waterModels())

    return _model_dict[model](
        molecule,
        box,
        angles,
        shell,
        ion_conc,
        is_neutral,
        is_aligned,
        match_water,
        work_dir,
        property_map,
    )


def spc(
    molecule=None,
    box=None,
    angles=3 * [_Angle(90, "degrees")],
    shell=None,
    ion_conc=0,
    is_neutral=True,
    is_aligned=False,
    match_water=True,
    work_dir=None,
    property_map={},
):
    """
    Add SPC solvent.

    Parameters
    ----------

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
               :class:`Molecule <BioSimSpace._SireWrappers.Molecules>`, \
               :class:`System <BioSimSpace._SireWrappers.System>`
        A molecule, or container/system of molecules.

    box : [:class:`Length <BioSimSpace.Types.Length>`]
        A list containing the box size in each dimension.

    angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
        A list containing the angles between the box vectors: yz, xz, and xy.

    shell : :class:`Length` <BioSimSpace.Types.Length>`
        Thickness of the water shell around the solute. Note that the
        base length of the resulting box must be at least twice as large
        as the cutoff used by the chosen molecular dynamics engine. As such,
        the shell option is often unsuitable for small molecules.

    ion_conc : float
        The ion concentration in (mol per litre).

    is_neutral : bool
        Whether to neutralise the system.

    is_aligned : bool
        Whether to align the principal axes of the molecule to those of the
        solvent box.

    match_water : bool
        Whether to update the naming of existing water molecules to match the
        expected convention for GROMACS, which is used as the solvation engine.

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The solvated molecular system.
    """

    if _gmx_exe is None or _gmx_path is None:
        raise _MissingSoftwareError(
            "'BioSimSpace.Solvent.spc' is not supported. "
            "Please install GROMACS (http://www.gromacs.org)."
        )

    # Validate arguments.
    molecule, box, angles, shell, work_dir, property_map = _validate_input(
        "spc",
        molecule,
        box,
        angles,
        shell,
        ion_conc,
        is_neutral,
        is_aligned,
        match_water,
        work_dir,
        property_map,
    )

    # Create the solvated system.
    return _solvate(
        molecule,
        box,
        angles,
        shell,
        "spc",
        3,
        ion_conc,
        is_neutral,
        is_aligned,
        match_water,
        work_dir=work_dir,
        property_map=property_map,
    )


def spce(
    molecule=None,
    box=None,
    angles=3 * [_Angle(90, "degrees")],
    shell=None,
    ion_conc=0,
    is_neutral=True,
    is_aligned=False,
    match_water=True,
    work_dir=None,
    property_map={},
):
    """
    Add SPC/E solvent.

    Parameters
    ----------

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
               :class:`Molecule <BioSimSpace._SireWrappers.Molecules>`, \
               :class:`System <BioSimSpace._SireWrappers.System>`
        A molecule, or container/system of molecules.

    box : [:class:`Length <BioSimSpace.Types.Length>`]
        A list containing the box size in each dimension.

    angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
        A list containing the angles between the box vectors: yz, xz, and xy.

    shell : :class:`Length` <BioSimSpace.Types.Length>`
        Thickness of the water shell around the solute. Note that the
        base length of the resulting box must be at least twice as large
        as the cutoff used by the chosen molecular dynamics engine. As such,
        the shell option is often unsuitable for small molecules.

    ion_conc : float
        The ion concentration in (mol per litre).

    is_neutral : bool
        Whether to neutralise the system.

    is_aligned : bool
        Whether to align the principal axes of the molecule to those of the
        solvent box.

    match_water : bool
        Whether to update the naming of existing water molecules to match the
        expected convention for GROMACS, which is used as the solvation engine.

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The solvated molecular system.
    """

    if _gmx_exe is None:
        raise _MissingSoftwareError(
            "'BioSimSpace.Solvent.spce' is not supported. "
            "Please install GROMACS (http://www.gromacs.org)."
        )

    # Validate arguments.
    molecule, box, angles, shell, work_dir, property_map = _validate_input(
        "spce",
        molecule,
        box,
        angles,
        shell,
        ion_conc,
        is_neutral,
        is_aligned,
        match_water,
        work_dir,
        property_map,
    )

    # Create the solvated system.
    return _solvate(
        molecule,
        box,
        angles,
        shell,
        "spce",
        3,
        ion_conc,
        is_neutral,
        is_aligned,
        match_water,
        work_dir=work_dir,
        property_map=property_map,
    )


def tip3p(
    molecule=None,
    box=None,
    angles=3 * [_Angle(90, "degrees")],
    shell=None,
    ion_conc=0,
    is_neutral=True,
    is_aligned=False,
    match_water=True,
    work_dir=None,
    property_map={},
):
    """
    Add TIP3P solvent.

    Parameters
    ----------

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
               :class:`Molecule <BioSimSpace._SireWrappers.Molecules>`, \
               :class:`System <BioSimSpace._SireWrappers.System>`
        A molecule, or container/system of molecules.

    box : [:class:`Length <BioSimSpace.Types.Length>`]
        A list containing the box size in each dimension.

    angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
        A list containing the angles between the box vectors: yz, xz, and xy.

    shell : :class:`Length` <BioSimSpace.Types.Length>`
        Thickness of the water shell around the solute. Note that the
        base length of the resulting box must be at least twice as large
        as the cutoff used by the chosen molecular dynamics engine. As such,
        the shell option is often unsuitable for small molecules.

    ion_conc : float
        The ion concentration in (mol per litre).

    is_neutral : bool
        Whether to neutralise the system.

    is_aligned : bool
        Whether to align the principal axes of the molecule to those of the
        solvent box.

    match_water : bool
        Whether to update the naming of existing water molecules to match the
        expected convention for GROMACS, which is used as the solvation engine.

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The solvated molecular system.
    """

    if _gmx_exe is None:
        raise _MissingSoftwareError(
            "'BioSimSpace.Solvent.tip3p' is not supported. "
            "Please install GROMACS (http://www.gromacs.org)."
        )

    # Validate arguments.
    molecule, box, angles, shell, work_dir, property_map = _validate_input(
        "tip3p",
        molecule,
        box,
        angles,
        shell,
        ion_conc,
        is_neutral,
        is_aligned,
        match_water,
        work_dir,
        property_map,
    )

    # Create the solvated system.
    return _solvate(
        molecule,
        box,
        angles,
        shell,
        "tip3p",
        3,
        ion_conc,
        is_neutral,
        is_aligned,
        match_water,
        work_dir=work_dir,
        property_map=property_map,
    )


def tip4p(
    molecule=None,
    box=None,
    angles=3 * [_Angle(90, "degrees")],
    shell=None,
    ion_conc=0,
    is_neutral=True,
    is_aligned=False,
    match_water=True,
    work_dir=None,
    property_map={},
):
    """
    Add TIP4P solvent.

    Parameters
    ----------

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
               :class:`Molecule <BioSimSpace._SireWrappers.Molecules>`, \
               :class:`System <BioSimSpace._SireWrappers.System>`
        A molecule, or container/system of molecules.

    box : [:class:`Length <BioSimSpace.Types.Length>`]
        A list containing the box size in each dimension.

    angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
        A list containing the angles between the box vectors: yz, xz, and xy.

    shell : :class:`Length` <BioSimSpace.Types.Length>`
        Thickness of the water shell around the solute. Note that the
        base length of the resulting box must be at least twice as large
        as the cutoff used by the chosen molecular dynamics engine. As such,
        the shell option is often unsuitable for small molecules.

    ion_conc : float
        The ion concentration in (mol per litre).

    is_neutral : bool
        Whether to neutralise the system.

    is_aligned : bool
        Whether to align the principal axes of the molecule to those of the
        solvent box.

    match_water : bool
        Whether to update the naming of existing water molecules to match the
        expected convention for GROMACS, which is used as the solvation engine.

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The solvated molecular system.
    """

    if _gmx_exe is None:
        raise _MissingSoftwareError(
            "'BioSimSpace.Solvent.tip4p' is not supported. "
            "Please install GROMACS (http://www.gromacs.org)."
        )

    # Validate arguments.
    molecule, box, angles, shell, work_dir, property_map = _validate_input(
        "tip4p",
        molecule,
        box,
        angles,
        shell,
        ion_conc,
        is_neutral,
        is_aligned,
        match_water,
        work_dir,
        property_map,
    )

    # Return the solvated system.
    return _solvate(
        molecule,
        box,
        angles,
        shell,
        "tip4p",
        4,
        ion_conc,
        is_neutral,
        is_aligned,
        match_water,
        work_dir=work_dir,
        property_map=property_map,
    )


def tip5p(
    molecule=None,
    box=None,
    angles=3 * [_Angle(90, "degrees")],
    shell=None,
    ion_conc=0,
    is_neutral=True,
    is_aligned=False,
    match_water=True,
    work_dir=None,
    property_map={},
):
    """
    Add TIP5P solvent.

    Parameters
    ----------

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
               :class:`Molecule <BioSimSpace._SireWrappers.Molecules>`, \
               :class:`System <BioSimSpace._SireWrappers.System>`
        A molecule, or container/system of molecules.

    box : [:class:`Length <BioSimSpace.Types.Length>`]
        A list containing the box size in each dimension.

    angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
        A list containing the angles between the box vectors: yz, xz, and xy.

    shell : :class:`Length` <BioSimSpace.Types.Length>`
        Thickness of the water shell around the solute. Note that the
        base length of the resulting box must be at least twice as large
        as the cutoff used by the chosen molecular dynamics engine. As such,
        the shell option is often unsuitable for small molecules.

    ion_conc : float
        The ion concentration in (mol per litre).

    is_neutral : bool
        Whether to neutralise the system.

    is_aligned : bool
        Whether to align the principal axes of the molecule to those of the
        solvent box.

    match_water : bool
        Whether to update the naming of existing water molecules to match the
        expected convention for GROMACS, which is used as the solvation engine.

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The solvated molecular system.
    """

    if _gmx_exe is None:
        raise _MissingSoftwareError(
            "'BioSimSpace.Solvent.tip5p' is not supported. "
            "Please install GROMACS (http://www.gromacs.org)."
        )

    # Validate arguments.
    molecule, box, angles, shell, work_dir, property_map = _validate_input(
        "tip5p",
        molecule,
        box,
        angles,
        shell,
        ion_conc,
        is_neutral,
        is_aligned,
        match_water,
        work_dir,
        property_map,
    )

    # Return the solvated system.
    return _solvate(
        molecule,
        box,
        angles,
        shell,
        "tip5p",
        5,
        ion_conc,
        is_neutral,
        is_aligned,
        match_water,
        work_dir=work_dir,
        property_map=property_map,
    )


def _validate_input(
    model,
    molecule,
    box,
    angles,
    shell,
    ion_conc,
    is_neutral,
    is_aligned,
    match_water,
    work_dir,
    property_map,
):
    """
    Internal function to validate function arguments.

    Parameters
    ----------

    model : str
        The name of the water model.

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
                :class:`Molecule <BioSimSpace._SireWrappers.Molecules>`, \
                :class:`System <BioSimSpace._SireWrappers.System>`
        A molecule, or container/system of molecules.

    box : [:class:`Length <BioSimSpace.Types.Length>`]
        A list containing the box size in each dimension.

    angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
        A list containing the angles between the box vectors: yz, xz, and xy.

    shell : :class:`Length` <BioSimSpace.Types.Length>`
        Thickness of the water shell around the solute. Note that the
        base length of the resulting box must be at least twice as large
        as the cutoff used by the chosen molecular dynamics engine. As such,
        the shell option is often unsuitable for small molecules.

    ion_conc : float
        The ion concentration in (mol per litre).

    is_neutral : bool
        Whether to neutralise the system.

    is_aligned : bool
        Whether to align the principal axes of the molecule to those of the
        solvent box.

    match_water : bool
        Whether to update the naming of existing water molecules to match the
        expected convention for GROMACS, which is used as the solvation engine.

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    (molecule, box, angles, shell, work_dir, property_map) : tuple
        The validated input arguments.
    """

    # Whether to check the box size.
    check_box = True

    # Validate the molecule and create a local copy called _molecule to ensure
    # that the passed molecule is preserved.
    if molecule is not None:
        if isinstance(molecule, _Molecule):
            _molecule = _Molecule(molecule)
        elif isinstance(molecule, _Molecules):
            _molecule = molecule.toSystem()
        elif isinstance(molecule, _System):
            _molecule = _System(molecule)
        else:
            raise TypeError(
                "'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule' "
                "'BioSimSpace._SireWrappers.Molecules', or 'BioSimSpace._SireWrappers.System'"
            )

        # Try to extract the box dimensions from the system.
        if isinstance(_molecule, _System) and box is None and shell is None:
            try:
                check_box = False
                prop = property_map.get("space", "space")
                box = molecule._sire_object.property(prop).dimensions()
                # Convert to a list of Length objects.
                box = [_Length(box[0], "A"), _Length(box[1], "A"), _Length(box[2], "A")]
            except:
                raise ValueError(
                    "The system has no box information. Please use "
                    "the 'box' keyword argument."
                )
        else:
            if box is None and shell is None:
                raise ValueError("Missing 'box' keyword argument!")

        # Warn the user if any of the molecules contain structural ions
        # parameterised for a different water model.
        if isinstance(_molecule, _System):
            for mol in _molecule:
                ion_water_model = mol._ion_water_model
                if ion_water_model is not None and ion_water_model != model:
                    _warnings.warn(
                        "Mismatch with water model used to parameterise "
                        f"structural ions: '{ion_water_model}'"
                    )
                    break
        else:
            ion_water_model = _molecule._ion_water_model
            if ion_water_model is not None and ion_water_model != model:
                _warnings.warn(
                    "Mismatch with water model used to parameterise "
                    f"structural ions: '{ion_water_model}'"
                )

    else:
        _molecule = None

        if box is None:
            raise ValueError("Missing 'box' keyword argument!")

        if shell is not None:
            _warnings.warn("Ignoring 'shell' keyword argument as solute is missing.")
            shell = None

    if box is not None:
        # Convert tuple to list.
        if isinstance(box, tuple):
            box = list(box)

        # Convert Coordinate to list.
        if isinstance(box, _Coordinate):
            box = [box.x(), box.y(), box.z()]

        # Validate.
        if len(box) != 3:
            raise ValueError("The 'box' must have x, y, and z size information.")
        else:
            if not all(isinstance(x, _Length) for x in box):
                raise ValueError(
                    "The box dimensions must be of type 'BioSimSpace.Types.Length'"
                )
            if not all(x.value() >= 0 for x in box):
                raise ValueError("All box dimensions must be greater than zero.")

    if angles is not None:
        # Convert tuple to list.
        if isinstance(angles, tuple):
            angles = list(angles)

        # Validate.
        if len(angles) != 3:
            raise ValueError("'angles' must have three components: yz, xz, and xy.")
        else:
            if not all(isinstance(x, _Angle) for x in angles):
                raise ValueError(
                    "The angle between box vectors must be of type 'BioSimSpace.Types.Angle'"
                )
    else:
        # Default to periodic box.
        angles = 3 * [_Angle(90, "degrees")]

    if shell is not None:
        if not isinstance(shell, _Length):
            raise ValueError("'shell' must must be of type 'BioSimSpace.Types.Length'")

        if box is not None:
            _warnings.warn(
                "Ignoring 'box' keyword argument as 'shell' takes precedence."
            )

        # Work out the box size based on axis-aligned bounding box.
        # We take the maximum dimension as the base length of our box.
        base_length = max(2 * molecule._getAABox().halfExtents())

        # Now add the shell thickness.
        base_length = _Length(base_length, "A") + shell

        # If we need to add ions, make sure the box is at least 2.54 nanometers
        # wide, i.e. twice the rlist cutoff used by GROMACS protocols.
        if base_length < _Length(2.54, "nm"):
            base_length = _Length(2.54, "nm")

        # Create the dimensions for a cubic box.
        box = 3 * [base_length]

    # Check that the ion concentration is valid.
    if not isinstance(ion_conc, float) and not type(ion_conc) is int:
        raise TypeError("'ion_conc' must be of type 'int' or 'float'.")
    elif ion_conc < 0:
        raise ValueError("'ion_conc' cannot be negative!")

    if not isinstance(is_neutral, bool):
        raise TypeError("'is_neutral' must be of type 'bool'.")

    if not isinstance(is_aligned, bool):
        raise TypeError("'is_aligned' must be of type 'bool'.")

    if not isinstance(match_water, bool):
        raise TypeError("'match_water' must be of type 'bool'.")

    # Check that the working directory is valid.
    if work_dir is not None and not isinstance(work_dir, str):
        raise TypeError("'work_dir' must be of type 'str'")

    # Check that the property map is valid.
    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    # Check that the box is large enough to hold the molecule.
    if check_box:
        if (
            molecule is not None
            and shell is None
            and not _check_box_size(molecule, box, property_map)
        ):
            raise ValueError("The 'box' is not large enough to hold the 'molecule'")

    return (_molecule, box, angles, shell, work_dir, property_map)


def _solvate(
    molecule,
    box,
    angles,
    shell,
    model,
    num_point,
    ion_conc,
    is_neutral,
    is_aligned,
    match_water,
    work_dir=None,
    property_map={},
):
    """
    Internal function to add solvent using 'gmx solvate'.

    Parameters
    ----------

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
               :class:`System <BioSimSpace._SireWrappers.System>`
        A molecule, or system of molecules.

    box : [:class:`Length <BioSimSpace.Types.Length>`]
        A list containing the box size in each dimension.

    angles : [:class:`Angle <BioSimSpace.Types.Angle>`]
        A list containing the angles between the box vectors: yz, xz, and xy.

    shell : :class:`Length` <BioSimSpace.Types.Length>`
        Thickness of the water shell around the solute.

    model : str
        The name of the water model.

    num_point : int
        The number of atoms in the water model.

    ion_conc : float
        The ion concentration in (mol per litre).

    is_neutral : bool
        Whether to neutralise the system.

    is_aligned : bool
        Whether to align the principal axes of the molecule to those of the
        solvent box.

    match_water : bool
        Whether to update the naming of existing water molecules to match the
        expected convention for GROMACS, which is used as the solvation engine.

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The solvated system.
    """

    if molecule is not None:
        # Get the axis aligned bounding box.
        aabox_min, aabox_max = molecule.getAxisAlignedBoundingBox()

        # Work out the aabox center.
        center = [
            0.5 * (aabox_max[x] + aabox_min[x]).angstroms().value() for x in range(0, 3)
        ]

        # Generate a TriclinicBox based on the box magnitudes and angles.
        triclinic_box = _TriclinicBox(
            box[0].angstroms().value(),
            box[1].angstroms().value(),
            box[2].angstroms().value(),
            angles[0].degrees().value() * _degree,
            angles[1].degrees().value() * _degree,
            angles[2].degrees().value() * _degree,
        )

        # Work out the center of the triclinic cell.
        box_center = triclinic_box.cellMatrix() * _Vector(0.5, 0.5, 0.5)

        # Work out the offset between the molecule and box centers.
        shift = [
            _Length(box_center[x].value() - center[x], "Angstrom") for x in range(0, 3)
        ]

        # Center the solute in the box.
        molecule.translate(shift)

        # Intitialise the original waters.
        original_waters = None

        if isinstance(molecule, _System):
            # Store the existing water molecules.
            original_waters = molecule.getWaterMolecules()

            # Reformat all of the water molecules so that they match the
            # expected GROMACS topology template, but flagged as crystal.
            molecule._set_water_topology(
                "GROMACS", is_crystal=True, property_map=property_map
            )

            # Extract the crystal waters.
            crystal_waters = molecule.getWaterMolecules()

            # Make sure the water molecules are at the end of the topology
            # since gmx genion requires that they are contiguous.
            molecule.removeWaterMolecules()
            molecule = molecule + crystal_waters

    # Create the working directory.
    work_dir = _Utils.WorkDir(work_dir)

    # Write to 6dp, unless precision is specified by use.
    _property_map = property_map.copy()
    if "precision" not in _property_map:
        _property_map["precision"] = _SireBase.wrap(6)

    # Run the solvation in the working directory.
    with _Utils.cd(work_dir):
        # First, generate a box file corresponding to the requested geometry.
        if molecule is not None:
            # Write the molecule/system to a GRO files.
            _property_map["crystal_water"] = "XTL"
            _IO.saveMolecules(
                "input",
                molecule,
                "gro87",
                match_water=False,
                property_map=_property_map,
            )

        # We need to create a dummy input file with no molecule in it.
        else:
            with open("input.gro", "w") as file:
                file.write("BioSimSpace_System\n")
                file.write("    0\n")
                file.write("   0.00000  0.00000  0.00000\n")

        # Create the editconf command.
        command = (
            "%s editconf -f input.gro -bt triclinic" % _gmx_exe
            + " -box %f %f %f"
            % (
                box[0].nanometers().value(),
                box[1].nanometers().value(),
                box[2].nanometers().value(),
            )
            + " -angles %f %f %f"
            % (
                angles[0].degrees().value(),
                angles[1].degrees().value(),
                angles[2].degrees().value(),
            )
            + f" -noc -{'no' if not is_aligned else ''}princ -o box.gro"
        )

        with open("README.txt", "w") as file:
            # Write the command to file.
            file.write("# gmx editconf was run with the following command:\n")
            file.write("%s\n" % command)

        # Create files for stdout/stderr.
        stdout = open("editconf.out", "w")
        stderr = open("editconf.err", "w")

        # Select the "System" if choosing to align.
        proc_echo = _subprocess.Popen(
            ["echo", "0"], shell=False, stdout=_subprocess.PIPE
        )
        # Run the command as a subprocess.
        proc = _subprocess.Popen(
            _Utils.command_split(command),
            shell=False,
            stdin=proc_echo.stdout,
            stdout=stdout,
            stderr=stderr,
        )
        proc.wait()
        proc_echo.stdout.close()

        stdout.close()
        stderr.close()

        # gmx doesn't return sensible error codes, so we need to check that
        # the expected output was generated.
        if not _os.path.isfile("box.gro"):
            raise RuntimeError(
                "'gmx editconf failed to generate required box! "
                + "Check your lattice vectors and angles."
            )

        # Create the gmx command.
        if num_point == 3:
            mod = "spc216"
        else:
            mod = model
        command = "%s solvate -cs %s" % (_gmx_exe, mod)

        # Add the shell information.
        if molecule is not None and shell is not None:
            command += " -shell %f" % shell.nanometers().value()

        command += " -cp box.gro -o output.gro"

        with open("README.txt", "a") as file:
            # Write the command to file.
            file.write("\n# gmx solvate was run with the following command:\n")
            file.write("%s\n" % command)

        # Create files for stdout/stderr.
        stdout = open("solvate.out", "w")
        stderr = open("solvate.err", "w")

        # Run gmx solvate as a subprocess.
        proc = _subprocess.run(
            _Utils.command_split(command), shell=False, stdout=stdout, stderr=stderr
        )
        stdout.close()
        stderr.close()

        # gmx doesn't return sensible error codes, so we need to check that
        # the expected output was generated.
        if not _os.path.isfile("output.gro"):
            raise RuntimeError("'gmx solvate failed to generate output!")

        # Extract the water lines from the GRO file.
        water_lines = []
        with open("output.gro", "r") as file:
            # Only search lines that weren't part of the existing molecule.
            if molecule is None:
                num_atoms = 0
            else:
                num_atoms = molecule.nAtoms()
            for line in file.readlines()[num_atoms + 2 :]:
                if _re.search("SOL", line):
                    # Store the SOL atom record.
                    water_lines.append(line)

            # Add any box information. This is the last line in the GRO file.
            water_lines.append(line)

        # Write a GRO file that contains only the water atoms.
        if len(water_lines) - 1 > 0:
            with open("water.gro", "w") as file:
                file.write("BioSimSpace %s water box\n" % model.upper())
                file.write("%d\n" % (len(water_lines) - 1))

                for line in water_lines:
                    file.write("%s" % line)
        else:
            raise ValueError(
                "No water molecules were generated. Try increasing "
                "the 'box' size or 'shell' thickness."
            )

        # Create a TOP file for the water model. By default we use the Amber03
        # force field to generate a dummy topology for the water model.
        with open("water_ions.top", "w") as file:
            file.write("#define FLEXIBLE 1\n\n")
            file.write("; Include AmberO3 force field\n")
            file.write('#include "amber03.ff/forcefield.itp"\n\n')
            file.write("; Include %s water topology\n" % model.upper())
            file.write('#include "amber03.ff/%s.itp"\n\n' % model)
            file.write("; Include ions\n")
            file.write('#include "amber03.ff/ions.itp"\n\n')
            file.write("[ system ] \n")
            file.write("BioSimSpace %s water box\n\n" % model.upper())
            file.write("[ molecules ] \n")
            file.write(";molecule name    nr.\n")
            file.write("SOL               %d\n" % ((len(water_lines) - 1) / num_point))

        # Load the water box.
        water = _IO.readMolecules(
            ["water.gro", "water_ions.top"], property_map=_property_map
        )

        # Create a new system by adding the water to the original molecule.
        if molecule is not None:
            if isinstance(molecule, _System):
                system = molecule + water
            else:
                system = molecule.toSystem() + water

            # Add all of the water box properties to the new system.
            for prop in water._sire_object.propertyKeys():
                prop = _property_map.get(prop, prop)

                # Add the space property from the water system.
                system._sire_object.setProperty(prop, water._sire_object.property(prop))
        else:
            system = water

        # Now we add ions to the system and neutralise the charge.
        if ion_conc > 0 or is_neutral:
            try:
                # Write the molecule + water system to file.
                _IO.saveMolecules(
                    "solvated",
                    system,
                    "gro87",
                    match_water=False,
                    property_map=_property_map,
                )
                _IO.saveMolecules(
                    "solvated",
                    system,
                    "grotop",
                    match_water=False,
                    property_map=_property_map,
                )
            except Exception as e:
                msg = (
                    "Failed to write GROMACS topology file. "
                    "Is your molecule parameterised?"
                )
                if _isVerbose():
                    raise IOError(msg) from e
                else:
                    raise IOError(msg) from None

            # First write an mdp file.
            with open("ions.mdp", "w") as file:
                file.write("; Neighbour searching\n")
                file.write("cutoff-scheme           = Verlet\n")
                file.write("rlist                   = 1.1\n")
                file.write("pbc                     = xyz\n")
                file.write("verlet-buffer-tolerance = -1\n")
                file.write("\n; Electrostatics\n")
                file.write("coulombtype             = cut-off\n")
                file.write("\n; VdW\n")
                file.write("rvdw                    = 1.0\n")

            # Create the grompp command.
            command = (
                "%s grompp -f ions.mdp -po ions.out.mdp -c solvated.gro -p solvated.top -o ions.tpr"
                % _gmx_exe
            )

            with open("README.txt", "a") as file:
                # Write the command to file.
                file.write("\n# gmx grompp was run with the following command:\n")
                file.write("%s\n" % command)

            # Create files for stdout/stderr.
            stdout = open("grommp.out", "w")
            stderr = open("grommp.err", "w")

            # Run grompp as a subprocess.
            proc = _subprocess.run(
                _Utils.command_split(command), shell=False, stdout=stdout, stderr=stderr
            )
            stdout.close()
            stderr.close()

            # Flag whether to break out of the ion adding stage.
            is_break = False

            # Check for the tpr output file.
            if not _os.path.isfile("ions.tpr"):
                if shell is None:
                    raise RuntimeError(
                        "'gmx grommp' failed to generate the required output for "
                        "'gmx genion'. Perhaps your box is too small?"
                    )
                else:
                    is_break = True
                    _warnings.warn(
                        "Unable to achieve target ion concentration, try using "
                        "'box' option instead of 'shell'."
                    )

            # Only continue if grommp was successful. This allows us to skip the remainder
            # of the code if the ion addition failed when the 'shell' option was chosen, i.e.
            # because the estimated simulation box was too small.
            if not is_break:
                is_break = False

                # The ion concentration is unset.
                if ion_conc == 0:
                    # Get the current molecular charge.
                    charge = system.charge()

                    # Round to the nearest integer value.
                    charge = round(charge.value())

                    # Create the genion command.
                    command = (
                        "%s genion -s ions.tpr -o solvated_ions.gro -p solvated.top -neutral"
                        % _gmx_exe
                    )

                    # Add enough counter ions to neutralise the charge.
                    if charge > 0:
                        command += " -nn %d" % abs(charge)
                    else:
                        command += " -np %d" % abs(charge)
                else:
                    # Create the genion command.
                    command = (
                        "%s genion -s ions.tpr -o solvated_ions.gro -p solvated.top -%s -conc %f"
                        % (_gmx_exe, "neutral" if is_neutral else "noneutral", ion_conc)
                    )

                with open("README.txt", "a") as file:
                    # Write the command to file.
                    file.write("\n# gmx genion was run with the following command:\n")
                    file.write("%s\n" % command)

                # Create files for stdout/stderr.
                stdout = open("genion.out", "w")
                stderr = open("genion.err", "w")

                # Run genion as a subprocess.
                proc_echo = _subprocess.Popen(
                    ["echo", "SOL"], shell=False, stdout=_subprocess.PIPE
                )
                proc = _subprocess.Popen(
                    _Utils.command_split(command),
                    shell=False,
                    stdin=proc_echo.stdout,
                    stdout=stdout,
                    stderr=stderr,
                )
                proc.wait()
                proc_echo.stdout.close()
                stdout.close()
                stderr.close()

                # Check for the output GRO file.
                if not _os.path.isfile("solvated_ions.gro"):
                    if shell is None:
                        raise RuntimeError(
                            "'gmx genion' failed to add ions! Perhaps your box is too small?"
                        )
                    else:
                        is_break = True
                        _warnings.warn(
                            "Unable to achieve target ion concentration, try using "
                            "'box' option instead of 'shell'."
                        )

                if not is_break:
                    # Counters for the number of SOL, NA, and CL atoms.
                    num_sol = 0
                    num_na = 0
                    num_cl = 0

                    # We now need to loop through the GRO file to extract
                    # the lines corresponding to water or ion atoms.
                    water_ion_lines = []

                    with open("solvated_ions.gro", "r") as file:
                        if molecule is None:
                            num_atoms = 0
                        # Make sure we don't don't search for ions from the
                        # original system.
                        else:
                            num_atoms = molecule.nAtoms()
                        for line in file.readlines()[num_atoms + 2 :]:
                            # This is a Sodium atom.
                            if _re.search("NA", line):
                                water_ion_lines.append(line)
                                num_na += 1

                            # This is a Chlorine atom.
                            elif _re.search("CL", line):
                                water_ion_lines.append(line)
                                num_cl += 1

                            # This is a water atom.
                            elif _re.search("SOL", line):
                                water_ion_lines.append(line)
                                num_sol += 1

                    # Add any box information. This is the last line in the GRO file.
                    water_ion_lines.append(line)

                    # Write a GRO file that contains only the water and ion atoms.
                    if len(water_ion_lines) - 1 > 0:
                        with open("water_ions.gro", "w") as file:
                            file.write("BioSimSpace %s water box\n" % model.upper())
                            file.write("%d\n" % (len(water_ion_lines) - 1))

                            for line in water_ion_lines:
                                file.write("%s" % line)

                    # Ions have been added. Update the TOP file for the water model
                    # with the new atom counts.
                    if num_na > 0 or num_cl > 0:
                        with open("water_ions.top", "w") as file:
                            file.write("#define FLEXIBLE 1\n\n")
                            file.write("; Include AmberO3 force field\n")
                            file.write('#include "amber03.ff/forcefield.itp"\n\n')
                            file.write("; Include %s water topology\n" % model.upper())
                            file.write('#include "amber03.ff/%s.itp"\n\n' % model)
                            file.write("; Include ions\n")
                            file.write('#include "amber03.ff/ions.itp"\n\n')
                            file.write("[ system ] \n")
                            file.write("BioSimSpace %s water box\n\n" % model.upper())
                            file.write("[ molecules ] \n")
                            file.write(";molecule name    nr.\n")
                            file.write("SOL               %d\n" % (num_sol / num_point))
                            if num_na > 0:
                                file.write("NA                %d\n" % num_na)
                            if num_cl > 0:
                                file.write("CL                %d\n" % num_cl)

                    # Load the water/ion box.
                    water_ions = _IO.readMolecules(["water_ions.gro", "water_ions.top"])

                    # Create a new system by adding the water and ions to the original molecule.
                    if molecule is not None:
                        if isinstance(molecule, _System):
                            # Preserve the original water naming.
                            if match_water is False:
                                # Flag each original water as non searchable.
                                flagged_waters = []
                                for water in original_waters:
                                    water._sire_object = (
                                        water._sire_object.edit()
                                        .setProperty(
                                            "is_non_searchable_water",
                                            _SireBase.wrap(True),
                                        )
                                        .molecule()
                                        .commit()
                                    )
                                    flagged_waters.append(water)
                                molecule.removeWaterMolecules()
                                system = molecule + flagged_waters + water_ions
                            # Update crystal waters to match the standard water topology.
                            else:
                                molecule._set_water_topology(
                                    "GROMACS", property_map=property_map
                                )
                                system = molecule + water_ions
                        else:
                            system = molecule.toSystem() + water_ions

                        # Add all of the system properties from the water molecules
                        # to the new system.
                        for prop in water_ions._sire_object.propertyKeys():
                            prop = _property_map.get(prop, prop)
                            system._sire_object.setProperty(
                                prop, water_ions._sire_object.property(prop)
                            )

                    else:
                        system = water_ions

        # Store the name of the water model as a system property.
        system._sire_object.setProperty("water_model", _SireBase.wrap(model))

    return system


def _check_box_size(molecule, box, property_map={}):
    """
    Internal function to check that box is big enough for the molecule.

    Parameters
    ----------

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, \
               :class:`System <BioSimSpace._SireWrappers.System>`
        A molecule, or system of molecules.

    box : [:class:`Length <BioSimSpace.Types.Length>`]
        A list containing the box size in each dimension.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    is_okay : True
        Whether the box is large enough.
    """

    # Get the axis-aligned bounding box of the molecule/system.
    aabox = molecule._getAABox(property_map)

    # Calculate the box size in each dimension, storing each component as a
    # length in Angstroms.
    mol_box = [_Length(2 * x, " A") for x in aabox.halfExtents()]

    # Make sure the box is big enough in each dimension.
    for len1, len2 in zip(box, mol_box):
        if len1 < len2:
            return False

    # We made it this far, all dimensions are large enough.
    return True


def _rename_water_molecule(molecule):
    """
    Internal function to rename residues/atoms in a water molecule to match
    the naming conventions used by GROMACS.

    Parameters
    ----------

    molecule : Sire.Mol.Molecule
        A Sire Molecule object.

    Returns
    -------

    molecule : Sire.Mol.Molecule
        The updated Sire Molecule object.
    """

    # Make the molecule editable.
    molecule = molecule.edit()

    # In GROMACS, all water molecules must be given the residue label "SOL".
    # We extract all of the waters from the system and relabel the
    # residues as appropriate.
    #
    # We need to work out what to do if existing water molecules don't contain
    # all of the required atoms, e.g. if we have crystal water oxygen sites
    # from a PDB file, or if the system is to be solvated with a 4- or 5-point
    # water model (where we need to add virtual sites).

    # Update the molecule with the new residue name.
    molecule = (
        molecule.residue(_SireMol.ResIdx(0)).rename(_SireMol.ResName("SOL")).molecule()
    )

    # Index for the hydrogen atoms.
    hydrogen_idx = 1

    # Gromacs water models use HW1/HW2 for hydrogen atoms at OW for water.
    for atom in molecule.atoms():
        try:
            # Hydrogen.
            if atom.property("element") == _SireMol.Element("H"):
                molecule = (
                    molecule.atom(atom.number())
                    .rename(_SireMol.AtomName("HW%d" % hydrogen_idx))
                    .molecule()
                )
                hydrogen_idx += 1
            # Oxygen.
            elif atom.property("element") == _SireMol.Element("O"):
                molecule = (
                    molecule.atom(atom.number())
                    .rename(_SireMol.AtomName("OW"))
                    .molecule()
                )

        # Otherwise, try to infer the element from the atom name.
        except:
            # Strip all digits from the name.
            name = "".join([x for x in atom.name().value() if not x.isdigit()])

            # Remove any whitespace.
            name = name.replace(" ", "")

            # Try to infer the element.
            element = _SireMol.Element.biologicalElement(name)

            # Hydrogen.
            if element == _SireMol.Element("H"):
                molecule = (
                    molecule.atom(atom.number())
                    .rename(_SireMol.AtomName("HW%d" % hydrogen_idx))
                    .molecule()
                )
                hydrogen_idx += 1
            # Oxygen.
            elif element == _SireMol.Element("O"):
                molecule = (
                    molecule.atom(atom.number())
                    .rename(_SireMol.AtomName("OW"))
                    .molecule()
                )

    # Commit and return the updated molecule.
    return molecule.commit()


# Create a list of the water models names.
# This needs to come after all of the solvation functions.
_models = []  # List of water models (actual names).
_models_lower = []  # List of lower case names.
_model_dict = {}  # Mapping between lower case names and functions.
import sys as _sys

_namespace = _sys.modules[__name__]
for _var in dir():
    if _var[0] != "_" and _var != "solvate" and _var[0] != "M":
        _models.append(_var)
        _models_lower.append(_var.lower())
        _model_dict[_var.lower()] = getattr(_namespace, _var)
del _namespace
del _sys
del _var


def waterModels():
    """
    Return a list of the supported water models.

    Returns
    -------

    models : [str]
       A list of the supported water models.
    """
    return _models
