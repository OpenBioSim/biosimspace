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

"""
Making biomolecular simulation a breeze!

A collection of tools that makes it easy to write robust and interoperable
molecular workflow components.

www.biosimspace.org
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = [
    "Align",
    "Box",
    "Convert",
    "FreeEnergy",
    "Gateway",
    "IO",
    "Metadynamics",
    "MD",
    "Node",
    "Notebook",
    "Parameters",
    "Process",
    "Protocol",
    "Solvent",
    "Stream",
    "Trajectory",
    "Types",
    "Units",
]

# Disable NumPy warnings.
import warnings as _warnings

_warnings.filterwarnings("ignore", module="numpy")

# Determine whether we're being imported from a Jupyter notebook.
try:
    _shell = get_ipython().__class__.__name__
    if _shell == "ZMQInteractiveShell":
        _is_notebook = True  # Jupyter notebook or qtconsole
    elif _shell == "TerminalInteractiveShell":
        _is_notebook = False  # Terminal running IPython
    else:
        _is_notebook = False  # Other type (?)
    del _shell
except NameError:
    _is_notebook = False  # Probably standard Python interpreter

# Determine whether we're being run interactively.
try:
    _shell = get_ipython().__class__.__name__
    if _shell == "ZMQInteractiveShell":
        _is_interactive = True  # Jupyter notebook or qtconsole
    elif _shell == "TerminalInteractiveShell":
        _is_interactive = True  # Terminal running IPython
    else:
        _is_interactive = False  # Other type (?)
    del _shell
except NameError:
    _is_interactive = False  # Probably standard Python interpreter

# Default to non-verbose error messages, unless the 'BSS_VERBOSE_ERRORS'
# environment variable is set to '1' (this allows verbose to be set before
# import, so that we can see verbose messages if there are any problems
# while importing BioSimSpace)
from os import environ as _environ

_is_verbose = "BSS_VERBOSE_ERRORS" in _environ and _environ["BSS_VERBOSE_ERRORS"] == "1"


def setVerbose(verbose):
    """
    Set verbosity of error messages.

    Parameters
    ----------

    verbose : bool
        Whether to print verbose error messages.
    """
    if type(verbose) is not bool:
        raise TypeError("'verbose' must be of type 'bool'.")

    global _is_verbose
    _is_verbose = verbose


def tutorialUrl():
    """Return the base URL for example files used in the tests and tutorials."""
    return "https://biosimspace.openbiosim.org/m"


def _isVerbose():
    """
    Whether verbose error messages are active.

    Returns
    -------

    is_verbose : bool
        Whether verbose error messages are active.
    """
    global _is_verbose
    return _is_verbose


# Check to see if AMBERHOME is set.
if "AMBERHOME" in _environ:
    _amber_home = _environ.get("AMBERHOME")
else:
    _amber_home = None

# Check to see if GROMACS is installed.
from sire.legacy import Base as _SireBase
from os import path as _path

# First, let the user tell us where to find GROMACS. This
# assumes that gromacs is installed in $GROMACSHOME/bin/gmx.
_gmx_exe = None
if "GROMACSHOME" in _environ:
    try:
        _gmx_exe = _SireBase.findExe(
            "%s/bin/gmx" % _environ.get("GROMACSHOME")
        ).absoluteFilePath()
    except:
        try:
            _gmx_exe = _SireBase.findExe(
                "%s/bin/gmx_mpi" % _environ.get("GROMACSHOME")
            ).absoluteFilePath()
        except:
            pass

if _gmx_exe is None:
    # The user has not told us where it is, so need to look in $PATH.
    try:
        _gmx_exe = _SireBase.findExe("gmx").absoluteFilePath()
    except:
        try:
            _gmx_exe = _SireBase.findExe("gmx_mpi").absoluteFilePath()
        except:
            pass

del _environ
del _SireBase

_gmx_path = None
_gmx_version = None

# Try using the GROMACS exe to get the version and data directory.
if _gmx_exe is not None:
    import shlex as _shlex
    import subprocess as _subprocess

    # Generate the shell command. (Run gmx -version.)
    _command = "%s -version" % _gmx_exe

    from ._Utils import command_split

    # Run the command.
    _proc = _subprocess.run(
        command_split(_command),
        shell=False,
        text=True,
        stdout=_subprocess.PIPE,
        stderr=_subprocess.PIPE,
    )

    del _command
    del command_split

    # Get the data prefix.
    if _proc.returncode == 0:
        for _line in _proc.stdout.split("\n"):
            # Extract the "Data prefix" from the output.
            if "Data prefix" in _line:
                _gmx_path = _line.split(":")[1].strip() + "/share/gromacs/top"
            # Extract the version from the output.
            elif "GROMACS version" in _line:
                _gmx_version = float(_line.split(":")[1].split("-")[0])
                break
        del _line

        # Check that the topology file directory exists.
        if _gmx_path is not None:
            if not _path.isdir(_gmx_path):
                _gmx_path = None

    del _path
    del _proc
    del _shlex
    del _subprocess

import lazy_import as _lazy_import

Align = _lazy_import.lazy_module("BioSimSpace.Align")
Box = _lazy_import.lazy_module("BioSimSpace.Box")
Convert = _lazy_import.lazy_module("BioSimSpace.Convert")
FreeEnergy = _lazy_import.lazy_module("BioSimSpace.FreeEnergy")
Gateway = _lazy_import.lazy_module("BioSimSpace.Gateway")
IO = _lazy_import.lazy_module("BioSimSpace.IO")
Metadynamics = _lazy_import.lazy_module("BioSimSpace.Metadynamics")
MD = _lazy_import.lazy_module("BioSimSpace.MD")
Node = _lazy_import.lazy_module("BioSimSpace.Node")
Notebook = _lazy_import.lazy_module("BioSimSpace.Notebook")
Parameters = _lazy_import.lazy_module("BioSimSpace.Parameters")
Process = _lazy_import.lazy_module("BioSimSpace.Process")
Protocol = _lazy_import.lazy_module("BioSimSpace.Protocol")
Solvent = _lazy_import.lazy_module("BioSimSpace.Solvent")
Stream = _lazy_import.lazy_module("BioSimSpace.Stream")
Trajectory = _lazy_import.lazy_module("BioSimSpace.Trajectory")
Types = _lazy_import.lazy_module("BioSimSpace.Types")
Units = _lazy_import.lazy_module("BioSimSpace.Units")

del _lazy_import

from . import _version

__version__ = _version.get_versions()["version"]
del _version
