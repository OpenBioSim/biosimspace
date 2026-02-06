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

"""Functionality for caching molecular files to avoid re-writing."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["clearCache", "disableCache", "enableCache"]

import collections as _collections


class _FixedSizeOrderedDict(_collections.OrderedDict):
    """A utility class to implement a fixed-sized cache."""

    def __init__(self, *args, max=2, **kwargs):
        """
        Constructor.

        Parameters
        ----------

        max : float
            The maximum size in GB.
        """
        import sys as _sys

        # Work out the approximate maximum number of atoms.
        if max > 0:
            self._max_atoms = int((max * 1e9) / (9 * _sys.getsizeof(float())))
        else:
            self._max_atoms = 0

        # Store the total number of atoms.
        self._num_atoms = 0

        super().__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        _collections.OrderedDict.__setitem__(self, key, value)
        self._num_atoms += value[0].nAtoms()
        if self._max_atoms > 0:
            if self._num_atoms > self._max_atoms:
                key, value = self.popitem(False)
                self._num_atoms -= value[0].nAtoms()

    def __delitem__(self, key):
        value = self[key]
        self._num_atoms -= value[0].nAtoms()
        _collections.OrderedDict.__delitem__(self, key)


# Initialise a "cache" dictionary. This maps a key of the system UID, file format
# and excluded properties a value of the system and file path. When saving to a
# given format, we can then to see if a matching system has previously been written
# to the same format, allowing us to re-use the existing file.
_cache = _FixedSizeOrderedDict()

# Whether to use the cache.
_use_cache = True


def clearCache():
    """
    Clear the file cache.
    """
    global _cache
    _cache = _FixedSizeOrderedDict()


def disableCache():
    """
    Disable the file cache.
    """
    global _use_cache
    _use_cache = False


def enableCache():
    """
    Enable the file cache.
    """
    global _use_cache
    _use_cache = True


def _cache_active():
    """
    Internal helper function to check whether the cache is active.
    """
    global _use_cache
    return _use_cache


def _check_cache(
    system,
    format,
    filebase,
    match_water=True,
    property_map={},
    excluded_properties=[],
    skip_water=True,
    **kwargs,
):
    """
    Check whether a Sire system has previously been written to the specified format.

    Parameters
    ----------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The system.

    format : str
        The molecular file format.

    filebase : str
        The file base to copy the file to.

    match_water : bool
        Whether to update the naming of water molecules to match the expected
        convention for the chosen file format. This is useful when a system
        is being saved to a different file format to that from which it was
        loaded.

    property_map : dict
        A dictionary that maps system "properties" to their user
        defined values. This allows the user to refer to properties
        with their own naming scheme, e.g. { "charge" : "my-charge" }

    excluded_properties : [str]
        A list of properties to exclude when comparing systems when checking
        the file cache.

    skip_water : bool
        Whether to skip water molecules when comparing systems.

    Returns
    -------

    extension : str
        The extension for cached file. False if no file was found.
    """
    from .._SireWrappers import System as _System
    import shutil as _shutil
    import os as _os

    # Validate input.

    if not isinstance(system, _System):
        raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")

    if not isinstance(format, str):
        raise TypeError("'format' must be of type 'str'")

    if not isinstance(filebase, str):
        raise TypeError("'filebase' must be of type 'str'")

    if not isinstance(match_water, bool):
        raise TypeError("'match_water' must be of type 'bool'")

    if not isinstance(excluded_properties, (list, tuple)):
        raise TypeError("'excluded_properties' must be a list of 'str' types.")

    if not all(isinstance(x, str) for x in excluded_properties):
        raise TypeError("'excluded_properties' must be a list of 'str' types.")

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'.")

    if not isinstance(skip_water, bool):
        raise TypeError("'skip_water' must be of type 'bool'.")

    global _cache

    # Create the key.
    key = (
        system._sire_object.uid().to_string(),
        format,
        _compress_molnum_key(str(system._mol_nums)),
        str(set(excluded_properties)),
        str(match_water),
        str(skip_water),
    )

    # Get the existing file path and MD5 hash from the cache.
    try:
        (prev_system, path, original_hash) = _cache[key]
    except:
        return False

    # Whether the cache entry is still valid.
    cache_valid = True

    # Is this system the same as the previous?
    if not system.isSame(
        prev_system,
        excluded_properties=excluded_properties,
        property_map0=property_map,
        property_map1=property_map,
        skip_water=skip_water,
    ):
        cache_valid = False

    # Make sure the file still exists.
    if not _os.path.exists(path):
        cache_valid = False
    # Make sure the MD5 sum is still the same.
    else:
        current_hash = _get_md5_hash(path)
        if current_hash != original_hash:
            cache_valid = False

    # If the cache isn't valid, delete the entry and return False.
    if not cache_valid:
        if key in _cache:
            del _cache[key]
        return False

    # Copy the old file to the new location.
    else:
        # Get the file extension.
        ext = _os.path.splitext(path)[1]

        # Add the extension to the file base.
        new_path = filebase + ext

        # Copy the file to the new location.
        try:
            _shutil.copyfile(path, new_path)
        except _shutil.SameFileError:
            pass
        except:
            del _cache[key]
            return False

        return ext


def _update_cache(
    system,
    format,
    path,
    excluded_properties=[],
    match_water=True,
    skip_water=True,
    **kwargs,
):
    """
    Update the file cache when a new system is written to a specified format.

    Parameters
    ----------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The system.

    format : str
        The molecular file format.

    path : str
        The path to the file.

    excluded_properties : [str]
        A list of properties to exclude when comparing systems when checking

    match_water : bool
        Whether to update the naming of water molecules to match the expected
        convention for the chosen file format. This is useful when a system
        is being saved to a different file format to that from which it was
        loaded.

    skip_water : bool
        Whether to skip water molecules when comparing systems.
    """
    from .._SireWrappers import System as _System
    import os as _os

    # Validate input.

    if not isinstance(system, _System):
        raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")

    if not isinstance(format, str):
        raise TypeError("'format' must be of type 'str'")

    if not isinstance(excluded_properties, (list, tuple)):
        raise TypeError("'excluded_properties' must be a list of 'str' types.")

    if not isinstance(path, str):
        raise TypeError("'path' must be of type 'str'")

    if not _os.path.exists(path):
        raise IOError(f"File does not exist: '{path}'")

    if not all(isinstance(x, str) for x in excluded_properties):
        raise TypeError("'excluded_properties' must be a list of 'str' types.")

    if not isinstance(match_water, bool):
        raise TypeError("'match_water' must be of type 'bool'")

    if not isinstance(skip_water, bool):
        raise TypeError("'skip_water' must be of type 'bool'.")

    global _cache

    # Convert to an absolute path.
    path = _os.path.abspath(path)

    # Get the MD5 checksum for the file.
    hash = _get_md5_hash(path)

    # Create the key.
    key = (
        system._sire_object.uid().to_string(),
        format,
        _compress_molnum_key(str(system._mol_nums)),
        str(set(excluded_properties)),
        str(match_water),
        str(skip_water),
    )

    # Update the cache.
    _cache[key] = (system.copy(), path, hash)


def _get_md5_hash(path):
    """
    Internal helper function to return the MD5 checksum for a file.

    Returns
    -------

    hash : hashlib.HASH
    """
    import hashlib as _hashlib

    # Get the MD5 hash of the file. Process in chunks in case the file is too
    # large to process.
    hash = _hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash.update(chunk)

    return hash.hexdigest()


def _compress_molnum_key(str):
    """
    Internal helper function to compress the MolNum list section of the key.
    """
    return str.replace("MolNum(", "").replace(")", "")
