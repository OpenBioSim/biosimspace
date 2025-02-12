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

"""Functionality for converting between molecular representations."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = [
    "smiles",
    "supportedFormats",
    "to",
    "toBioSimSpace",
    "toOpenMM",
    "toRDKit",
    "toSire",
]

from .._Utils import _try_import, _have_imported

_openmm = _try_import("openmm")

import os as _os

from rdkit.Chem.rdchem import Mol as _RDMol

import rdkit.Chem as _Chem

from sire import convert as _sire_convert
from sire import smiles as _sire_smiles

import sire.legacy.Base as _SireBase
import sire.legacy.Mol as _SireMol
import sire.legacy.System as _SireSystem
import sire.legacy.Vol as _SireVol

import sire.system as _NewSireSystem

from .._Exceptions import ConversionError as _ConversionError
from .. import IO as _IO
from .. import _SireWrappers


def smiles(
    smiles_string, add_hydrogens=True, generate_coordinates=True, property_map={}
):
    """
    Generate a BioSimSpace Molecule froma SMILES string.

    Parameters
    ----------

    smiles_string : str
        The molecule in SMILES string format.

    add_hydrogens : bool
        Whether or not to automatically add hydrogens.

    generate_coordinates : bool
       Whether or not to automatically generate coordinates. Note that
       generating the coordinates will automatically switch on addition
       of hydrogens.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
        A BioSimSpace molecule.
    """

    if not isinstance(smiles_string, str):
        raise TypeError("'smiles_string' must be of type 'str'.")

    if not isinstance(add_hydrogens, bool):
        raise TypeError("'add_hydrogens' must be of type 'bool'.")

    if not isinstance(generate_coordinates, bool):
        raise TypeError("'generate_coordinates' must be of type 'bool'.")

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'.")

    # Strip whitespace.
    smiles_string = smiles_string.replace(" ", "")

    try:
        molecule = _SireWrappers.Molecule(
            _sire_smiles(
                smiles_string,
                add_hydrogens=add_hydrogens,
                generate_coordinates=generate_coordinates,
                map=property_map,
            )
        )

        # Rename the molecule with the original SMILES string.
        # Since the name is written to topology file formats, we
        # need to ensure that it doesn't start with an [ character,
        # which would break GROMACS.
        if smiles_string.startswith("["):
            name = f"smiles:{smiles_string}"
            edit_mol = molecule._sire_object.edit()
            edit_mol = edit_mol.rename(name).molecule()
            molecule._sire_object = edit_mol.commit()

        return molecule

    except:
        raise _ConversionError(
            f"Unable to create a BioSimSpace Molecule from SMILES string: {smiles_string}"
        )


def supportedFormats():
    """
    Return a list of the supported formats.

    Returns
    -------

    formats : [str]
        The supported formats.
    """
    return _sire_convert.supported_formats()


def to(obj, format="biosimspace", property_map={}, **kwargs):
    """
    Convert an object to a specified format.

    Parameters
    ----------

    obj :
        The input object to convert.

    format : str
        The format to convert to.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    converted_obj :
       The object in the converted format.
    """

    # Validate the input.

    if not isinstance(format, str):
        raise TypeError("'format' must be of type 'str'.")

    # Convert to lower case and strip whitespace.
    format = format.lower().replace(" ", "")

    if not format in supportedFormats():
        raise ValueError(
            f"Unsupported format '{format}', options are: {', '.join(supportedFormats())}."
        )

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'.")

    # Check for force_stereo_inference in kwargs.
    if "force_stereo_inference" in kwargs:
        force_stereo_inference = kwargs["force_stereo_inference"]
        if not isinstance(force_stereo_inference, bool):
            raise TypeError("'force_stereo_inference' must be of type 'bool'.")
        property_map["force_stereo_inference"] = _SireBase.wrap(force_stereo_inference)

    # Special handling for OpenMM conversion. Currently this is a one-way (toOpenMM)
    # conversion only and is only supported for specific Sire and BioSimSpace types.
    if format == "openmm":
        if not _have_imported(_openmm):
            raise ConversionError(
                "Conversion to OpenMM format is currently not supported on this platform."
            )

        # If this is already an OpenMM context, then simply return it.
        if isinstance(obj, _openmm.openmm.Context):
            return obj

        # BioSimSpace objects.
        if isinstance(obj, _SireWrappers._sire_wrapper.SireWrapper):
            if not isinstance(obj, _SireWrappers.System):
                # Convert to a system where possible.
                try:
                    obj = obj.toSystem()
                except:
                    # Otherwise, convert residues/atoms to a molecule, then a system.
                    try:
                        obj = obj.toMolecule().toSystem()
                    except:
                        raise _ConversionError(
                            "Unable to convert object to OpenMM format!"
                        )

            # Get the user-defined space property name.
            prop = property_map.get("space", "space")

            # Try to get the space property. Use an infinite cartesian space
            # if none is present.
            try:
                space = obj._sire_object.property(prop)
            except:
                space = _SireVol.Cartesian()

            # Set a shared space property.
            obj._sire_object.addSharedProperty(prop, space)

            # Now try to convert the object to OpenMM format.
            try:
                return _sire_convert.to(
                    obj,
                    "openmm",
                    map=property_map,
                )
            except:
                raise _ConversionError("Unable to convert object to OpenMM format!")

        # Sire objects.
        elif isinstance(
            obj,
            (
                _NewSireSystem.System,
                _SireSystem.System,
                _SireMol.MoleculeGroup,
                _SireMol.Molecules,
                _SireMol.Molecule,
                _SireMol.Residue,
            ),
        ):
            try:
                # First convert the object to BioSimSpace format.
                obj = to(obj, format="biosimspace", property_map=property_map)

                # Now try converting to OpenMM format.
                return to(obj, format="openmm", property_map=property_map)

            except:
                raise _ConversionError("Unable to convert object to OpenMM format!")

        else:
            raise _ConversionError(
                f"Currently unable to convert object of type {type(obj)} to OpenMM format!"
            )

    # Now do some work to handle different types of input. Where possible, we
    # try to convert to the "expected" format, i.e. atom to atom, residue
    # to residue, etc. This means that we need to use sire.legacy.Mol.PartialMolecule
    # to extract the sub-components we require.

    # BioSimSpace objects, i.e. wrapped Sire objects.
    if isinstance(obj, _SireWrappers._sire_wrapper.SireWrapper):
        # Return the object itself.
        if format == "biosimspace":
            return obj

        # Return the underlying object.
        elif format == "sire":
            if isinstance(obj, _SireWrappers.Molecules):
                return _SireMol.SelectorMol(obj._sire_object)
            else:
                return obj._sire_object

        elif format == "rdkit":
            # Handle System and Molecules objects.
            if isinstance(obj, (_SireWrappers.System, _SireWrappers.Molecules)):
                try:
                    return _sire_convert.to(
                        obj._sire_object.molecules(), "rdkit", map=property_map
                    )
                except:
                    raise _ConversionError("Unable to convert object to RDKit format!")
            # Handle Residues and Atoms by extracting the PartialMolecule.
            else:
                try:
                    return _sire_convert.to(
                        _SireMol.PartialMolecule(obj._sire_object).extract(),
                        "rdkit",
                        map=property_map,
                    )
                except:
                    raise _ConversionError("Unable to convert object to RDKit format!")

    # Sire objects.
    elif isinstance(
        obj,
        (
            _NewSireSystem.System,
            _SireSystem.System,
            _SireMol.SelectorMol,
            _SireMol.MoleculeGroup,
            _SireMol.Molecules,
            _SireMol.Molecule,
            _SireMol.Residue,
            _SireMol.Atom,
        ),
    ):
        if format == "biosimspace":
            if isinstance(obj, _NewSireSystem.System):
                try:
                    return _SireWrappers.System(obj._system)
                except:
                    raise _ConversionError(
                        "Unable to convert object to BioSimSpace format!"
                    )

            if isinstance(obj, _SireSystem.System):
                try:
                    return _SireWrappers.System(obj)
                except:
                    raise _ConversionError(
                        "Unable to convert object to BioSimSpace format!"
                    )

            elif isinstance(obj, _SireMol.SelectorMol):
                try:
                    return _SireWrappers.Molecules(obj.toMoleculeGroup())
                except:
                    raise _ConversionError(
                        "Unable to convert object to BioSimSpace format!"
                    )

            elif isinstance(obj, _SireMol.Molecules):
                try:
                    molgrp = _SireMol.MoleculeGroup("all")
                    molgrp.add(obj)
                    return _SireWrappers.Molecules(molgrp)
                except:
                    raise _ConversionError(
                        "Unable to convert object to BioSimSpace format!"
                    )

            elif isinstance(obj, _SireMol.MoleculeGroup):
                try:
                    return _SireWrappers.Molecules(obj)
                except:
                    raise _ConversionError(
                        "Unable to convert object to BioSimSpace format!"
                    )

            elif isinstance(obj, _SireMol.Molecule):
                try:
                    return _SireWrappers.Molecule(obj)
                except:
                    raise _ConversionError(
                        "Unable to convert object to BioSimSpace format!"
                    )

            elif isinstance(obj, _SireMol.Residue):
                try:
                    return _SireWrappers.Residue(obj)
                except:
                    raise _ConversionError(
                        "Unable to convert object to BioSimSpace format!"
                    )

            elif isinstance(obj, _SireMol.Atom):
                try:
                    return _SireWrappers.Atom(obj)
                except:
                    raise _ConversionError(
                        "Unable to convert object to BioSimSpace format!"
                    )

        elif format == "rdkit":
            if isinstance(
                obj, (_NewSireSystem.System, _SireSystem.System, _SireMol.MoleculeGroup)
            ):
                try:
                    return _sire_convert.to(obj.molecules(), "rdkit", map=property_map)
                except:
                    raise _ConversionError("Unable to convert object to RDKit format!")

            elif isinstance(obj, (_SireMol.Molecules, _SireMol.Molecule)):
                try:
                    return _sire_convert.to(obj, "rdkit", map=property_map)
                except:
                    raise _ConversionError("Unable to convert object to RDKit format!")

            elif isinstance(obj, (_SireMol.Residue, _SireMol.Atom)):
                try:
                    return _sire_convert.to(
                        _SireMol.PartialMolecule(obj).extract(),
                        "rdkit",
                        map=property_map,
                    )
                except:
                    raise _ConversionError("Unable to convert object to RDKit format!")

    # RDKit objects.
    elif isinstance(obj, _RDMol):
        try:
            # Convert the object, then return as a single residue, atom,
            # or entire molecule.
            new_obj = _sire_convert.to(obj, format, map=property_map)

            if format == "biosimspace":
                if new_obj.nAtoms() == 1:
                    return new_obj.getAtoms()[0]
                elif new_obj.nResidues() == 1:
                    return new_obj.getResidues()[0]
                else:
                    return new_obj

            elif format == "sire":
                if new_obj.nAtoms() == 1:
                    return new_obj.atoms()[0]
                if new_obj.nResidues() == 1:
                    return new_obj.residues()[0]
                else:
                    return new_obj
        except:
            raise _ConversionError(f"Unable to convert object to {format} format!")

    # Lists or tuples of objects.
    elif isinstance(obj, (list, tuple)):
        initial_type = type(obj[0])
        if not all(isinstance(x, initial_type) for x in obj):
            raise ValueError("Containers must contain objects of the same type!")

        if isinstance(obj[0], (_RDMol, _SireWrappers.Molecule)):
            try:
                new_obj = _sire_convert.to(obj, format, map=property_map)
                if format == "biosimspace":
                    new_obj = _SireWrappers.Molecules(new_obj)
                return new_obj
            except:
                raise _ConversionError(f"Unable to convert object to {format} format!")

        else:
            raise TypeError(f"Cannot convert object of type '{type(obj)}'")

    else:
        raise TypeError(f"Cannot convert from object of type '{type(obj)}'")


def toBioSimSpace(obj, property_map={}):
    """
    Convert an object to BioSimSpace format.

    Parameters
    ----------

    obj :
        The input object to convert.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    converted_obj :
       The object in BioSimSpace format.
    """
    return to(obj, format="biosimspace", property_map=property_map)


def toOpenMM(obj, property_map={}):
    """
    Convert an object to OpenMM format.

    Parameters
    ----------

    obj :
        The input object to convert.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    converted_obj :
       The object in RDKit format.
    """
    if _have_imported(_openmm):
        return to(obj, format="openmm", property_map=property_map)
    else:
        raise ConversionError(
            "Conversion to OpenMM format is currently not supported on this platform."
        )


def toRDKit(obj, force_stereo_inference=False, property_map={}):
    """
    Convert an object to RDKit format.

    Parameters
    ----------

    obj :
        The input object to convert.

    bool : force_stereo_inference
        Whether to force inference of stereochemistry, overriding any
        stereochemistry present in the input object. This is useful when
        the object has been loaded from a file with invalid stereochemistry.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    converted_obj :
       The object in OpenMM format.
    """

    if not isinstance(force_stereo_inference, bool):
        raise TypeError("'force_stereo_inference' must be of type 'bool'.")

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'.")

    property_map["force_stereo_inference"] = _SireBase.wrap(force_stereo_inference)

    return to(obj, format="rdkit", property_map=property_map)


def toSire(obj, property_map={}):
    """
    Convert an object to Sire format.

    Parameters
    ----------

    obj :
        The input object to convert.

    format : str
        The format to convert to.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    converted_obj :
       The object in Sire format.
    """
    return to(obj, format="sire", property_map=property_map)


def _to_rdkit(molecule, work_dir=_os.getcwd(), direct=True, property_map={}):
    """
    Internal function to convert two BioSimSpace molecules to RDKit format.
    This will go via an intermediate file format if direct conversion fails.

    Parameters
    ----------

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
        The BioSimSpace molecule.

    work_dir : str
        The path to the working directory.

    dire : bool
        Whether to use a direct conversion (True) or go via an intermediate
        file format (False).

    property_map : dict
        A dictionary that maps "properties" in 'molecule' to their user
        defined values. This allows the user to refer to properties
        with their own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    rdmol : rdkit.Chem.rdchem.Mol
        The molecule in RDKit format.
    """

    if not isinstance(molecule, _SireWrappers.Molecule):
        raise TypeError(
            "'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'."
        )

    if not isinstance(work_dir, str):
        raise TypeError("'work_dir' must be of type 'str'.")

    if not _os.path.exists(work_dir):
        raise ValueError(f"'work_dir' doesn't exist: {work_dir}")

    if not isinstance(direct, bool):
        raise TypeError("'direct' must be of type 'bool'.")

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'.")

    # First try to convert to RDKit format directly.
    try:
        if not direct:
            raise Exception
        else:
            rdmol = toRDKit(molecule, property_map=property_map)

    # If this fails, then go via an intermediate file format.
    except:
        try:
            is_sdf = False
            filebase = work_dir + "/tmp"

            # Try to go via SDF format to preserve bond orders.
            if molecule._sire_object.hasProperty("fileformat"):
                if "SDF" in molecule._sire_object.property("fileformat").value():
                    _IO.saveMolecules(
                        filebase, molecule, "SDF", property_map=property_map
                    )
                    rdmol = _Chem.MolFromMolFile(filebase + ".sdf", True, False, 0)

                    # Flag that we went via SDF format.
                    is_sdf = True

            # Convert via PDB format..
            if not is_sdf:
                _IO.saveMolecules(filebase, molecule, "PDB", property_map=property_map)
                rdmol = _Chem.MolFromPDBFile(filebase + ".pdb", True, False, 0)

        except:
            raise _ConversionError("Unable to convert molecule to RDKit format!")

    return rdmol
