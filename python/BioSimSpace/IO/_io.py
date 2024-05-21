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

"""Functionality for reading/writing molecular systems."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = [
    "expand",
    "fileFormats",
    "formatInfo",
    "readMolecules",
    "readPDB",
    "readPerturbableSystem",
    "saveMolecules",
    "savePerturbableSystem",
]

from collections import OrderedDict as _OrderedDict
from glob import glob as _glob
from io import StringIO as _StringIO

import json as _json
import os as _os
import shlex as _shlex
import shutil as _shutil
import sys as _sys
import subprocess as _subprocess
import warnings as _warnings

# Flag that we've not yet raised a warning about GROMACS not being installed.
_has_gmx_warned = False

import sire as _sire

from sire.legacy import Base as _SireBase
from sire.legacy import IO as _SireIO
from sire.legacy import Mol as _SireMol
from sire.legacy import System as _SireSystem

from .. import _amber_home
from .. import _gmx_path
from .. import _isVerbose
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import Molecule as _Molecule
from .._SireWrappers import Molecules as _Molecules
from .._SireWrappers import System as _System
from .. import _Utils

from ._file_cache import _check_cache
from ._file_cache import _update_cache
from ._file_cache import _cache_active


# Context manager for capturing stdout.
# Taken from:
# https://stackoverflow.com/questions/16571150/how-to-capture-stdout-output-from-a-python-function-call
class _Capturing(list):
    def __enter__(self):
        self._stdout = _sys.stdout
        _sys.stdout = self._stringio = _StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio
        _sys.stdout = self._stdout


# Capture the supported format information
with _Capturing() as format_info:
    print(r"%s" % _SireIO.MoleculeParser.supportedFormats())

# Create a list of the supported formats.
_formats = []

# Create a dictionary of format-description key:value pairs.
_formats_dict = _OrderedDict()

# Loop over the format information to populate the dictionary.
for index, line in enumerate(format_info):
    if "Parser" in line:
        format = line.split()[2]
        extensions = format_info[index + 1]
        description = format_info[index + 2]

        if format != "SUPPLEMENTARY":
            _formats.append(format)
            _formats_dict[format.replace(" ", "").upper()] = (format, description)

# Delete the redundant variables.
del format_info, index, line, format, extensions, description


def expand(base, path, suffix=None):
    """
    Expand the set of paths with the supplied base.

    Parameters
    ----------

    base : str
        The base to prepend to all paths.

    path : str, [str]
        The filename (or names) that will be prepended with the base.

    suffix : str
        An optional suffix to append to all files, e.g. ".bz2".

    Returns
    -------
    path : [str]
        The list of expanded filenames or URLs.
    """

    if not isinstance(base, str):
        raise TypeError("'base' must be of type 'str'")

    # Convert single values to a list.
    if isinstance(path, str):
        path = [path]

    if not isinstance(path, (list, tuple)) and not all(
        isinstance(x, str) for x in path
    ):
        raise TypeError("'path' must be a list of 'str' types.")

    if suffix is not None and not isinstance(suffix, str):
        raise TypeError("'suffix' must be of type 'str'")

    return _sire.expand(base, path, suffix=suffix)


def fileFormats():
    """
    Return a list of the supported formats.

    Returns
    -------

    formats : [str]
        A list of the support file formats.
    """
    return _formats


def formatInfo(format):
    """
    Return information for the specified file format.

    Parameters
    ----------

    format : str
        The file format.

    Returns
    -------

    info : str
        A description of the named file format.

    Examples
    --------

    Display information regarding the PDB format.

    >>> import BioSimSpace as BSS
    >>> BSS.formatInfo("PDB")

    Print information for each of the supported file formats.

    >>> import BioSimSpace as BSS
    >>> for format in BSS.IO.fileFormats:
    ...     BSS.IO.formatInfo(format)
    """

    try:
        return _formats_dict[format.replace(" ", "").upper()][1]
    except KeyError:
        print("Unsupported format: '%s'" % format)
        return None


def readPDB(id, pdb4amber=False, work_dir=None, show_warnings=False, property_map={}):
    """
    Read a molecular system from a Protein Data Bank (PDBP) ID in the RSCB PDB
    website.

    Parameters
    ----------

    id : str
        The PDB ID string, or path to a PDB file.

    pdb4amber : bool
        Whether to process the PDB file using pdb4amber. This reformats the file
        such that it can be handled by the AMBER suite of tools.

    show_warnings : bool
        Whether to show any warnings raised during parsing of the input files.

    work_dir : str
        The working directory used to run pdb4amber.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        A molecular system.

    Examples
    --------

    Create a molecular system from the deoxy human haemoglobin Protein
    Data Bank (PDB) record.

    >>> import BioSimSpace as BSS
    >>> system = BSS.IO.readPDB("1a3n")

    Create a molecular system from a PDB file on disk and re-format so that
    it is compatible with the AmberTools suite.
    Data Bank (PDB) record.

    >>> import BioSimSpace as BSS
    >>> system = BSS.IO.readPDB("file.pdb", pdb4amber=True)
    """

    if not isinstance(id, str):
        raise TypeError("'id' must be of type 'str'")

    if not isinstance(pdb4amber, bool):
        raise TypeError("'pdb4amber' must be of type 'bool'")

    if work_dir and not isinstance(work_dir, str):
        raise TypeError("'work_dir' must be of type 'str'")

    # Create the working directory.
    work_dir = _Utils.WorkDir(work_dir)

    # Path to a PDB file.
    if _os.path.isfile(id):
        pdb_file = _os.path.abspath(id)

    # This is a URL.
    elif id.startswith(("http", "www")):
        from sire._load import _resolve_path

        try:
            pdb_file = _resolve_path(id, directory=str(work_dir))[0]
        except:
            raise IOError(f"Unable to download PDB file: '{id}'")

    # ID from the Protein Data Bank.
    else:
        # Strip any whitespace from the PDB ID and convert to upper case.
        id = id.replace(" ", "").upper()

        # Use Sire to download the PDB.
        try:
            system = _patch_sire_load(id, directory=str(work_dir))
        except:
            raise IOError("Retrieval failed, invalid PDB ID: %s" % id)

        # Store the absolute path of the file.
        pdb_file = f"{work_dir}/{id}"

        # Save to PDB format.
        saveMolecules(pdb_file, _System(system), "pdb")

        # Add the extension.
        pdb_file += ".pdb"

    # Process the file with pdb4amber.
    if pdb4amber:
        # Check that pdb4amber exists.
        if _amber_home is None:
            raise _MissingSoftwareError(
                "Please install AmberTools for pdb4amber support: http://ambermd.org"
            )
        else:
            _pdb4amber_exe = "%s/bin/pdb4amber" % _amber_home
            if not _os.path.isfile(_pdb4amber_exe):
                raise IOError("Missing pdb4amber executable: '%s'" % _pdb4amber_exe)

                # Create the file prefix.
        prefix = work_dir + "/"

        # Create the pdb4amber command.
        command = "%s %s -o pdb4amber.pdb" % (_pdb4amber_exe, pdb_file)

        # Create files for stdout/stderr.
        stdout = open(prefix + "pdb4amber.out", "w")
        stderr = open(prefix + "pdb4amber.err", "w")

        # Run pdb4amber as a subprocess.
        proc = _subprocess.run(
            _Utils.command_split(command),
            cwd=str(work_dir),
            shell=False,
            stdout=stdout,
            stderr=stderr,
        )
        stdout.close()
        stderr.close()

        # Check that the output PDB file was generated.
        # the expected output was generated.
        if _os.path.isfile("%s/pdb4amber.pdb" % work_dir):
            pdb_file = "%s/pdb4amber.pdb" % work_dir
        else:
            raise IOError("pdb4amber failed!")

    # Read the file and return a molecular system.
    return readMolecules(
        pdb_file, show_warnings=show_warnings, property_map=property_map
    )


def readMolecules(
    files,
    make_whole=False,
    rotate_box=False,
    reduce_box=False,
    show_warnings=False,
    download_dir=None,
    property_map={},
):
    """
    Read a molecular system from file.

    Parameters
    ----------

    files : str, [str]
        A file name, or a list of file names. Note that the file names can
        be URLs, in which case the files will be downloaded and (if necessary)
        extracted before reading.

    make_whole : bool
        Whether to make molecules whole, i.e. unwrap those that are split across
        the periodic boundary.

    rotate_box : bool
        Rotate the box vectors of the system so that the first vector is
        aligned with the x-axis, the second vector lies in the x-y plane,
        and the third vector has a positive z-component. This is a requirement
        of certain molecular dynamics engines and is used for simulation
        efficiency. All vector properties of the system, e.g. coordinates,
        velocities, etc., will be rotated accordingly.

    reduce_box : bool
        Reduce the box vectors.

    show_warnings : bool
        Whether to show any warnings raised during parsing of the input files.

    download_dir : str
        The directory to download files to. If None, then a temporary directory
        will be created for you.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        A molecular system.

    Examples
    --------

    Load a molecular system from AMBER coordinate and topology files.

    >>> import BioSimSpace as BSS
    >>> system = BSS.IO.readMolecules(["ala.rst7", "ala.prm7"])

    Load the same system, but map the "charge" property to the key "my-charge".

    >>> import BioSimSpace as BSS
    >>> system = BSS.IO.readMolecules(["ala.rst7", "ala.prm7"], property_map={"charge" : "my-charge"})

    >>> import BioSimSpace as BSS
    >>> system = BSS.IO.readMolecules(["ala.rst7", "ala.prm7"])

    Load a molecular system from all of the files contained within a directory.

    >>> import BioSimSpace as BSS
    >>> system = BSS.IO.readMolecules("dir/*")

    Load a molecular system from GROMACS coordinate and topology files using
    a custom GROMACS topology directory.

    >>> import BioSimSpace as BSS
    >>> system = BSS.IO.readMolecules(["mol.gro87", "mol.grotop"], property_map={"GROMACS_PATH" : "/path/to/gromacs/topology"})
    """

    global _has_gmx_warned
    if _gmx_path is None and not _has_gmx_warned:
        _warnings.warn(
            "BioSimSpace.IO: Please install GROMACS (http://www.gromacs.org) "
            "for GROMACS topology file support."
        )
        _has_gmx_warned = True

    # Convert a single string to a list.
    if isinstance(files, str):
        files = [files]

    # Check that all arguments are of type 'str'.
    if isinstance(files, (list, tuple)):
        if not all(isinstance(x, str) for x in files):
            raise TypeError("'files' must be a list of 'str' types.")
        if len(files) == 0:
            raise ValueError("The list of input files is empty!")
        # Convert tuple to list.
        if isinstance(files, tuple):
            files = list(files)
    else:
        raise TypeError("'files' must be of type 'str', or a list of 'str' types.")

    # Glob all files to catch wildcards.
    new_files = []
    for file in files:
        if not file.startswith(("http", "www")):
            new_files += _glob(file)
        else:
            new_files.append(file)
    files = new_files

    # Validate the molecule unwrapping flag.
    if not isinstance(make_whole, bool):
        raise TypeError("'make_whole' must be of type 'bool'.")
    # Flag that we want' to unwrap molecules.
    if make_whole:
        property_map["make_whole"] = _SireBase.wrap(make_whole)

    # Validate the box rotation flag.
    if not isinstance(rotate_box, bool):
        raise TypeError("'rotate_box' must be of type 'bool'.")

    # Validate the box reduction flag.
    if not isinstance(reduce_box, bool):
        raise TypeError("'reduce_box' must be of type 'bool'.")

    # Validate the warning message flag.
    if not isinstance(show_warnings, bool):
        raise TypeError("'show_warnings' must be of type 'bool'.")

    # Validate the download directory.
    if download_dir is not None:
        if not isinstance(download_dir, str):
            raise TypeError("'download_dir' must be of type 'str'")

    # Create the download directory.
    download_dir = _Utils.WorkDir(download_dir)

    # Validate the map.
    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    # Add the GROMACS topology file path.
    if _gmx_path is not None and ("GROMACS_PATH" not in property_map):
        property_map["GROMACS_PATH"] = _gmx_path

    # Check that the files exist (if not a URL).
    for file in files:
        if not file.startswith(("http", "www")) and not _os.path.isfile(file):
            raise IOError("Missing input file: '%s'" % file)

    # Try to read the files and return a molecular system.
    try:
        system = _patch_sire_load(
            files,
            directory=str(download_dir),
            property_map=property_map,
            show_warnings=show_warnings,
        )
    except Exception as e0:
        if "There are no lead parsers!" in str(e0):
            # First check to see if the failure was due to the presence
            # of CMAP records.
            for file in files:
                try:
                    amber_prm = _SireIO.AmberPrm(file)
                except Exception as e1:
                    if "CMAP" in str(e1).upper():
                        msg = (
                            "Unable to parse AMBER topology file. CMAP "
                            "records are currently unsupported."
                        )
                        if _isVerbose():
                            raise IOError(msg) from e1
                        else:
                            raise IOError(msg) from None
                    elif "CHAMBER" in str(e1).upper():
                        msg = (
                            "Unable to parse AMBER topology file. "
                            "CHAMBER files are currently unsupported."
                        )
                        if _isVerbose():
                            raise IOError(msg) from e1
                        else:
                            raise IOError(msg) from None

            msg = (
                "Failed to read molecules from %s. "
                "It looks like you failed to include a topology file."
            ) % files
            if _isVerbose():
                raise IOError(msg) from e0
            else:
                raise IOError(msg) from None
        else:
            if "Incompatibility" in str(e0):
                msg = (
                    "Incompatibility between molecular information in files: %s" % files
                )
                if _isVerbose():
                    raise IOError(msg) from e0
                else:
                    raise IOError(msg) from None
            else:
                msg = "Failed to read molecules from: %s" % files
                if _isVerbose():
                    raise IOError(msg) from e0
                else:
                    raise IOError(msg) from None

    # Add a file format shared property.
    prop = property_map.get("fileformat", "fileformat")
    system.addSharedProperty(prop, system.property(prop))

    # Remove "space" and "time" shared properties since this causes incorrect
    # behaviour when extracting molecules and recombining them to make other
    # systems.
    try:
        # Space.
        prop = property_map.get("space", "space")
        space = system.property(prop)
        system.removeSharedProperty(prop)
        system.setProperty(prop, space)

        # Time.
        prop = property_map.get("time", "time")
        time = system.property(prop)
        system.removeSharedProperty(prop)
        system.setProperty(prop, time)
    except:
        pass

    # Wrap the Sire system.
    system = _System(system)

    # Rotate the box vectors.
    if rotate_box:
        system.rotateBoxVectors()

    # Reduce the box vectors.
    if reduce_box:
        system.reduceBoxVectors()

    return system


def saveMolecules(
    filebase, system, fileformat, match_water=True, property_map={}, **kwargs
):
    """
    Save a molecular system to file.

    Parameters
    ----------

    filebase : str
        The base name of the output files.

    system : :class:`System <BioSimSpace._SireWrappers.System>`, \
             :class:`Molecule< BioSimSpace._SireWrappers.Molecule>` \
             :class:`Molecule< BioSimSpace._SireWrappers.Molecules>`
        The molecular system.

    fileformat : str, [str]
        The file format (or formats) to save to.

    match_water : bool
        Whether to update the naming of water molecules to match the expected
        convention for the chosen file format. This is useful when a system
        is being saved to a different file format to that from which it was
        loaded.

    property_map : dict
        A dictionary that maps system "properties" to their user
        defined values. This allows the user to refer to properties
        with their own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    files : [str]
        The list of files that were generated.

    Examples
    --------

    Load a molecular system from AMBER coordinate and topology files then
    try to save it to all supported file formats.

    >>> import BioSimSpace as BSS
    >>> files = BSS.IO.expand(BSS.tutorialUrl(), ["ala.top", "ala.crd"], ".bz2")
    >>> system = BSS.IO.readMolecules(files)
    >>> for format in BSS.IO.fileFormats():
    ...     try:
    ...         BSS.IO.saveMolecules("test", system, format)
    ...     except:
    ...         print("Could not convert to format: '%s'" % format)

    Load a molecular system from AMBER coordinate and topology files then
    try to save it to GROMACS format, mapping and un-mapping the charge
    property along the way.

    >>> import BioSimSpace as BSS
    >>> files = BSS.IO.expand(BSS.tutorialUrl(), ["ala.top", "ala.crd"], ".bz2")
    >>> system = BSS.IO.readMolecules(files, property_map={"charge" : "my-charge"})
    >>> BSS.IO.saveMolecules("test", system, ["gro87", "grotop"], property_map={"charge" : "my-charge"})
    """

    global _has_gmx_warned
    if _gmx_path is None and not _has_gmx_warned:
        _warnings.warn(
            "BioSimSpace.IO: Please install GROMACS (http://www.gromacs.org) "
            "for GROMACS topology file support."
        )
        _has_gmx_warned = True

    # Check that the filebase is a string.
    if not isinstance(filebase, str):
        raise TypeError("'filebase' must be of type 'str'")

    # Convert to absolute path.
    if not _os.path.isabs(filebase):
        filebase = _os.path.abspath(filebase)

    # Check that that the system is of the correct type.

    # A System object.
    if isinstance(system, _System):
        pass
    # A Molecule object.
    elif isinstance(system, _Molecule):
        system = _System(system)
    elif isinstance(system, _Molecules):
        system = system.toSystem()
    # A list of Molecule objects.
    elif isinstance(system, list) and all(isinstance(x, _Molecule) for x in system):
        system = _System(system)
    # Invalid type.
    else:
        raise TypeError(
            "'system' must be of type 'BioSimSpace.SireWrappers.System', "
            "'BioSimSpace._SireWrappers.Molecule, 'BioSimSpace._SireWrappers.Molecules' "
            "or a list of 'BiSimSpace._SireWrappers.Molecule' types."
        )

    # Check that fileformat argument is of the correct type.

    # Convert to a list if a single string is passed.
    # We split on ',' since the user might pass system.fileFormats() as the argument.
    if isinstance(fileformat, str):
        fileformat = fileformat.split(",")
    # Lists and tuples are okay!
    elif isinstance(fileformat, (list, tuple)):
        pass
    else:
        raise TypeError("'fileformat' must be a 'str' or a 'list' of 'str' types.")

    # Make sure all items in list or tuple are strings.
    if not all(isinstance(x, str) for x in fileformat):
        raise TypeError("'fileformat' must be a 'str' or a 'list' of 'str' types.")

    # Validate the match_water flag.
    if not isinstance(match_water, bool):
        raise TypeError("'match_water' must be of type 'bool'.")

    # Make a list of the matched file formats.
    formats = []

    # Make sure that all of the formats are valid.
    for format in fileformat:
        try:
            f = _formats_dict[format.replace(" ", "").upper()][0]
            formats.append(f)
        except KeyError:
            raise ValueError(
                "Unsupported file format '%s'. Supported formats "
                "are: %s." % (format, str(_formats))
            )

    # Validate the map.
    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    # Copy the map.
    _property_map = property_map.copy()

    # Add the GROMACS topology file path.
    if _gmx_path is not None and ("GROMACS_PATH" not in _property_map):
        _property_map["GROMACS_PATH"] = _gmx_path

    # Get the directory name.
    dirname = _os.path.dirname(filebase)

    # Create the directory if it doesn't already exist.
    if not _os.path.isdir(dirname):
        _os.makedirs(dirname, exist_ok=True)

    # A list of the files that have been written.
    files = []

    # Save the system using each file format.
    for format in formats:
        # Copy an existing file if it exists in the cache.
        if _cache_active():
            ext = _check_cache(
                system,
                format,
                filebase,
                match_water=match_water,
                property_map=property_map,
                **kwargs,
            )
        else:
            ext = None
        if ext:
            files.append(_os.path.abspath(filebase + ext))
            continue

        # Add the file format to the property map.
        _property_map["fileformat"] = _SireBase.wrap(format)

        # Warn the user if any molecules are parameterised with a force field
        # that uses geometric combining rules. While we can write this to file
        # the information is lost on read.
        if format.upper() == "PRM7":
            # Get the name of the "forcefield" property.
            forcefield = _property_map.get("forcefield", "forcefield")

            # Loop over all molecules in the system.
            for mol in system.getMolecules():
                if mol._sire_object.hasProperty(forcefield):
                    if (
                        mol._sire_object.property(forcefield).combiningRules()
                        == "geometric"
                    ):
                        _warnings.warn(
                            "AMBER topology files do not support force fields that "
                            "use geometric combining rules, as this cannot be specified "
                            "in the file. When this file is re-read, then arithmetic "
                            "combining rules will be assumed."
                        )
                        # Exit after the first non-arithmetic molecule we encounter.
                        break

        # Write the file.
        try:
            # Make sure AMBER and GROMACS files have the expected water topology
            # and save GROMACS files with an extension such that they can be run
            # directly by GROMACS without needing to be renamed.
            if format.upper() == "PRM7":
                if match_water:
                    system_copy = system.copy()
                    system_copy._set_water_topology("AMBER", property_map=_property_map)
                else:
                    system_copy = system
                file = _SireIO.MoleculeParser.save(
                    system_copy._sire_object, filebase, _property_map
                )
            elif format.upper() == "GROTOP":
                if match_water:
                    system_copy = system.copy()
                    system_copy._set_water_topology(
                        "GROMACS", property_map=_property_map
                    )
                else:
                    system_copy = system
                    _property_map["skip_water"] = _SireBase.wrap(True)
                file = _SireIO.MoleculeParser.save(
                    system_copy._sire_object, filebase, _property_map
                )[0]
                new_file = file.replace("grotop", "top")
                _os.rename(file, new_file)
                file = [new_file]
            elif format.upper() == "GRO87":
                if match_water:
                    system_copy = system.copy()
                    system_copy._set_water_topology(
                        "GROMACS", property_map=_property_map
                    )
                else:
                    system_copy = system
                    _property_map["skip_water"] = _SireBase.wrap(True)
                # Write to 3dp by default, unless greater precision is
                # requested by the user.
                if "precision" not in _property_map:
                    _property_map["precision"] = _SireBase.wrap(3)
                file = _SireIO.MoleculeParser.save(
                    system_copy._sire_object, filebase, _property_map
                )[0]
                new_file = file.replace("gro87", "gro")
                _os.rename(file, new_file)
                file = [new_file]
            else:
                file = _SireIO.MoleculeParser.save(
                    system._sire_object, filebase, _property_map
                )

            files += file

            # If this is a new file, then add it to the cache.
            if _cache_active():
                _update_cache(
                    system, format, file[0], match_water=match_water, **kwargs
                )

        except Exception as e:
            msg = "Failed to save system to format: '%s'" % format
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

    # Return the list of files.
    return files


def savePerturbableSystem(filebase, system, property_map={}):
    """
    Save a system containing a perturbable molecule. This will be written in
    AMBER format, with a topology file for each end state of the perturbation,
    i.e. 'filebase0.prm7' and 'filebase1.prm7'. Coordinates for the end
    states are written to 'filebase0.rst7' and 'filebase1.rst7'.

    Parameters
    ----------

    filebase : str
        The base name of the output files.

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The molecular system.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }
    """

    # Check that the filebase is a string.
    if not isinstance(filebase, str):
        raise TypeError("'filebase' must be of type 'str'")

    # Check that the system is valid.
    if isinstance(system, _System):
        pass
    # A Molecule object.
    elif isinstance(system, _Molecule):
        system = _System(system)
    # A Molecules object.
    elif isinstance(system, _Molecules):
        system = system.toSystem()
    # A list of Molecule objects.
    elif isinstance(system, list) and all(isinstance(x, _Molecule) for x in system):
        system = _System(system)
    # Invalid type.
    else:
        raise TypeError(
            "'system' must be of type 'BioSimSpace.SireWrappers.System', "
            "'BioSimSpace._SireWrappers.Molecule, 'BioSimSpace._SireWrappers.Molecules' "
            "or a list of 'BiSimSpace._SireWrappers.Molecule' types."
        )

    # Validate the map.
    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    # Validate that there is a single perturbable molecule in the system.
    pert_mols = system.getPerturbableMolecules()
    if len(pert_mols) != 1:
        raise ValueError(
            "The 'system' must contain a single perturbable molecule. "
            f"Found {len(pert_mols)}!"
        )

    # Extract the molecule
    pert_mol = pert_mols[0]

    # Create a copy of the system for the lambda=0 and lambda=1 end states.
    system0 = system.copy()
    system1 = system.copy()

    # Update the perturbable molecule in each system.
    system0.updateMolecules(
        pert_mol._toRegularMolecule(property_map=property_map, is_lambda1=False)
    )
    system1.updateMolecules(
        pert_mol._toRegularMolecule(property_map=property_map, is_lambda1=True)
    )

    # Save the topology files.
    saveMolecules(filebase + "0", system0, "prm7")
    saveMolecules(filebase + "1", system1, "prm7")

    # Save the coordinate files.
    saveMolecules(filebase + "0", system0, "rst7")
    saveMolecules(filebase + "1", system1, "rst7")


def readPerturbableSystem(top0, coords0, top1, coords1, property_map={}):
    """
    Read a perturbable system from file.

    Parameters
    ----------

    top0 : str
        The path to the topology file for the lambda=0 end state.

    coords0 : str
        The path to the coordinate file for the lambda=0 end state.

    top1 : str
        The path to the topology file for the lambda=1 end state.

    coords1 : str
        The path to the coordinate file for the lambda=1 end state.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        A molecular system.
    """

    if not isinstance(top0, str):
        raise TypeError("'top0' must be of type 'str'.")

    if not isinstance(coords0, str):
        raise TypeError("'coords0' must be of type 'str'.")

    if not isinstance(top1, str):
        raise TypeError("'top1' must be of type 'str'.")

    if not isinstance(coords1, str):
        raise TypeError("'coords1' must be of type 'str'.")

    # Check that the coordinate and topology files can be parsed.

    prefixes = ("http", "www")

    # Don't validate URLs.
    if (
        not top0.startswith(prefixes)
        and not coords0.startswith(prefixes)
        and not top1.startswith(prefixes)
        and not coords1.startswith(prefixes)
    ):
        # lamba = 0 coordinates.
        try:
            _SireIO.AmberRst7(coords0)
        except Exception as e:
            msg = f"Unable to read lambda=0 coordinate file: {coords0}"
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # lamba = 1 coordinates.
        try:
            _SireIO.AmberRst7(coords1)
        except Exception as e:
            msg = f"Unable to read lambda=1 coordinate file: {coords1}"
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # lamba = 0 topology.
        try:
            parser = _SireIO.AmberPrm(top0)
        except Exception as e:
            msg = f"Unable to read lambda=0 topology file: {top0}"
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None
        if parser.isEmpty():
            raise ValueError(
                f"Unable to read topology file for lamba=0 end state: {top0}"
            )

        # lamba = 1 topology.
        try:
            parser = _SireIO.AmberPrm(top1)
        except Exception as e:
            msg = f"Unable to read lambda=1 topology file: {top1}"
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None
        if parser.isEmpty():
            raise ValueError(
                f"Unable to read topology file for lamba=1 end state: {top1}"
            )

    # Try loading the two end states.
    system0 = readMolecules([coords0, top0], property_map=property_map)
    system1 = readMolecules([coords1, top1], property_map=property_map)

    # Make sure the systems have the same number of molecules.
    if system0.nMolecules() != system1.nMolecules():
        raise ValueError("The two topologies contain a different number of molecules!")

    # Now loop through the molecules in each system to work out which
    # is the perturbable molecule. This will differ in the  'ambertype'
    # 'LJ' or 'charge' property.
    ambertype = property_map.get("ambertype", "ambertype")
    LJ = property_map.get("LJ", "LJ")
    charge = property_map.get("charge", "charge")
    has_pert = False
    for idx, (mol0, mol1) in enumerate(
        zip(system0.getMolecules(), system1.getMolecules())
    ):
        for atom0, atom1 in zip(mol0.getAtoms(), mol1.getAtoms()):
            if (
                atom0._sire_object.property(ambertype)
                != atom1._sire_object.property(ambertype)
                or atom0._sire_object.property(LJ) != atom1._sire_object.property(LJ)
                or atom0._sire_object.property(charge)
                != atom1._sire_object.property(charge)
            ):
                has_pert = True
                break
        if has_pert:
            break

    if not has_pert:
        raise ValueError("No perturbable molecule was found?")

    # Extract the perturbable molecule.
    pert_mol = system0[idx]

    # Extract and copy the Sire molecule.
    mol = pert_mol._sire_object.__deepcopy__()

    # Make the molecule editable.
    mol = mol.edit()

    # Rename all properties in the molecule for the lambda=0 end state,
    # e.g.: "prop" --> "prop0". Then delete all properties named "prop"
    # and "prop1".
    for prop in mol.propertyKeys():
        # See if this property exists in the user map.
        new_prop = property_map.get(prop, prop) + "0"

        # Copy the property using the updated name.
        mol = mol.setProperty(new_prop, mol.property(prop)).molecule()

        # Delete the redundant property.
        mol = mol.removeProperty(prop).molecule()

    # Now add the properties for the lambda=1 end state.
    mol1 = system1[idx]._sire_object
    for prop in mol1.propertyKeys():
        # See if this property exists in the user map.
        new_prop = property_map.get(prop, prop) + "1"

        # Copy the property using the updated name.
        mol = mol.setProperty(new_prop, mol1.property(prop)).molecule()

    # Flag that the molecule is perturbable.
    mol.setProperty("is_perturbable", _SireBase.wrap(True))

    # Get the two molecules.
    mol0 = system0[idx]._sire_object
    mol1 = system1[idx]._sire_object

    # Add the molecule0 and molecule1 properties.
    mol.setProperty("molecule0", mol0)
    mol.setProperty("molecule1", mol1)

    # Get the connectivity property name.
    conn_prop = property_map.get("connectivity", "connectivity")

    # Get the connectivity from the end states.
    conn0 = mol.property(conn_prop + "0")
    conn1 = mol.property(conn_prop + "1")

    # Check whether the connectivity is the same.
    if conn0 == conn1:
        # The connectivity is the same, so we can use the connectivity
        # from the lambda=0 end state.
        mol = mol.setProperty(conn_prop, conn0).molecule()

        # Delete the end state properties.
        mol = mol.removeProperty(conn_prop + "0").molecule()
        mol = mol.removeProperty(conn_prop + "1").molecule()

    # Reconstruct the intrascale matrices using the GroTop parser.
    intra0 = (
        _SireIO.GroTop(_Molecule(mol0).toSystem()._sire_object)
        .toSystem()[0]
        .property("intrascale")
    )
    intra1 = (
        _SireIO.GroTop(_Molecule(mol1).toSystem()._sire_object)
        .toSystem()[0]
        .property("intrascale")
    )

    # Set the "intrascale" properties.
    intrascale_prop = property_map.get("intrascale", "intrascale")
    mol.setProperty(intrascale_prop + "0", intra0)
    mol.setProperty(intrascale_prop + "1", intra0)

    # Commit the changes.
    mol = _Molecule(mol.commit())

    # Update the molecule in the original system.
    system0.updateMolecules(mol)

    # Remove "space" and "time" shared properties since this causes incorrect
    # behaviour when extracting molecules and recombining them to make other
    # systems.
    try:
        # Space.
        prop = property_map.get("space", "space")
        space = system0._sire_object.property(prop)
        system0._sire_object.removeSharedProperty(prop)
        system0._sire_object.setProperty(prop, space)

        # Time.
        prop = property_map.get("time", "time")
        time = system0._sire_object.property(prop)
        system0._sire_object.removeSharedProperty(prop)
        system0._sire_object.setProperty(prop, time)
    except:
        pass

    return system0


def _patch_sire_load(path, *args, show_warnings=True, property_map={}, **kwargs):
    """
    Load the molecular system at 'path'. This can be a filename
    of a URL. If it is a URL, then the file will be downloaded
    to the current directory and loaded from there.

    Parameters
    ----------

    path : str or list[str]
        The filename (or names) or the URL or URLS of the molecular
        system to load. This allows multiple paths to be input
        as some molecular file formats split molecular information
        across multiple files. Multiple paths can also be passed
        as multiple arguments to this function.

    log : (dict)
        Optional dictionary that you can pass in that will be populated
        with any error messages or warnings from the parsers as they
        attempt to load in the molecular data. This can be helpful
        in diagnosing why your file wasn't loaded.

    show_warnings : bool
        Whether or not to print out any warnings that are encountered
        when loading your file(s). This is default True, and may lead
        to noisy output. Set `show_warnings=False` to silence this output.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    directory : str
        Optional directory which will be used when creating any
        files (e.g. as a download from a URL or which unzipping files)

    Returns
    -------

    system : sire.legacy.System.System:
        The molecules that have been loaded are returned as
        a sire.legacy.System.System.
    """

    if type(path) is not list:
        paths = [path]
    else:
        paths = path

    for arg in args:
        paths.append(arg)

    if "log" in kwargs:
        log = kwargs["log"]
    else:
        log = {}

    if "directory" in kwargs:
        directory = kwargs["directory"]
    else:
        directory = "."

    if "silent" in kwargs:
        silent = kwargs["silent"]
    else:
        silent = False

    p = []

    for i in range(0, len(paths)):
        # resolve the paths, downloading as needed
        p += _sire._load._resolve_path(paths[i], directory=directory, silent=silent)

    paths = p

    if len(paths) == 0:
        raise IOError("No valid files specified. Nothing to load?")

    s = _sire.io.load_molecules(paths, map=_sire.base.create_map(property_map))

    return _sire._load._to_legacy_system(s)
