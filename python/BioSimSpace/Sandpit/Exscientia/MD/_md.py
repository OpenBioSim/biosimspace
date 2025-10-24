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

"""Functionality for configuring and driving molecular dynamics simulations."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["run"]

import os as _os

from sire.legacy import Base as _SireBase

from .. import _amber_home, _gmx_exe
from .._Exceptions import IncompatibleError as _IncompatibleError
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import System as _System
from .. import Process as _Process
from .. import Protocol as _Protocol

# A dictionary mapping MD engines to their executable names and GPU support.
#                engine        EXE               GPU
_md_engines = {
    "AMBER": {
        "pmemd.cuda.MPI": True,
        "pmemd.cuda": True,
        "pmemd.MPI": False,
        "pmemd": False,
        "sander": False,
    },
    "GROMACS": {"gmx": True, "gmx_mpi": True},
    "NAMD": {"namd2": False},
    "OPENMM": {"sire_python": True},
    "SOMD": {"somd": True},
}

# A dictionary reverse mapping MD engines to their supported Sire file extensions.
# Use SOMD as a fall-back where possible. Since we can't guarantee interconversion
# of potentials for CHARMM-PSF format input files, we restrict such simulations to
# only run using NAMD.
#                    EXTENSION        ENGINES
_file_extensions = {
    "PRM7,RST7": ["AMBER", "GROMACS", "OPENMM", "SOMD"],
    "PRM7,RST": ["AMBER", "GROMACS", "OPENMM", "SOMD"],
    "GroTop,Gro87": ["GROMACS", "AMBER", "OPENMM", "SOMD"],
    "GROTOP,GRO87": ["GROMACS", "AMBER", "OPENMM", "SOMD"],
    "PSF,PDB": ["NAMD"],
}

# Whether each engine supports free energy simulations. This dictionary needs to
# be updated as support for different engines is added.
_free_energy = {
    "AMBER": True,
    "GROMACS": True,
    "NAMD": False,
    "OPENMM": False,
    "SOMD": True,
}

# Whether each engine supports metadynamics simulations. This dictionary needs to
# be updated as support for different engines is added.
_metadynamics = {
    "AMBER": True,
    "GROMACS": True,
    "NAMD": False,
    "OPENMM": True,
    "SOMD": False,
}

# Whether each engine supports steered molecular dynamics simulations. This
# dictionary needs to # be updated as support for different engines is added.
_steering = {
    "AMBER": True,
    "GROMACS": True,
    "NAMD": False,
    "OPENMM": False,
    "SOMD": False,
}


def _find_md_engines(system, protocol, engine="AUTO", gpu_support=False):
    """
    Find molecular dynamics engines on the system that
    support the given protocol and GPU requirements.

    Parameters
    ----------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The molecular system.

    protocol : :class:`Protocol <BioSimSpace.Protocol>`
        The simulation protocol.

    engine : str
        The molecular dynamics engine to use. If "AUTO", then a matching
        engine will automatically be chosen.

    gpu_support : bool
        Whether the engine must have GPU support.

    Returns
    -------

    engines, exes : [ str ], [ str ]
       Lists containing the supported MD engines and executables.
    """

    # The input has already been validated in the run method, so no need
    # to re-validate here.

    # Get the file format of the molecular system.
    fileformat = system.fileFormat()

    # Make sure that this format is supported.
    if fileformat is None:
        if engine == "AUTO":
            raise ValueError("Cannot automatically choose an MD engine")
        engines = list(_md_engines.keys())
    elif not fileformat in _file_extensions:
        raise ValueError(
            "Cannot find an MD engine that supports format: %s" % fileformat
        )
    else:
        engines = _file_extensions[fileformat]

    # If engine != "AUTO", then check the chosen engine supports the file
    # format.
    md_engine = engine
    if md_engine != "AUTO":
        if md_engine not in _md_engines.keys():
            raise ValueError(f"The {engine} MD engine isn't supported!")
        if md_engine not in engines:
            raise ValueError(f"The {engine} MD engine doesn't format {fileformat}")
        else:
            # Just search for the chosen engine.
            engines = [md_engine]

    is_free_energy = False
    is_metadynamics = False
    is_steering = False

    if isinstance(protocol, _Protocol._FreeEnergyMixin):
        is_free_energy = True
    elif isinstance(protocol, _Protocol.Metadynamics):
        is_metadynamics = True
    elif isinstance(protocol, _Protocol.Steering):
        is_steering = True

    # Create a list to store all of the engines and executables.
    found_engines = []
    found_exes = []

    # Loop over each engine that supports the file format.
    for engine in engines:
        # Don't continue if the engine doesn't support the protocol.
        if (
            (not is_free_energy or _free_energy[engine])
            and (not is_metadynamics or _metadynamics[engine])
            and (not is_steering or _steering[engine])
        ):
            # Check whether this engine exists on the system and has the desired
            # GPU support.
            for exe, gpu in _md_engines[engine].items():
                # If the user has requested GPU support make sure the engine
                # supports it.
                if not gpu_support or gpu:
                    # AMBER
                    if engine == "AMBER":
                        # Search AMBERHOME, if set.
                        if _amber_home is not None:
                            _exe = "%s/bin/%s" % (_amber_home, exe)
                            if _os.path.isfile(_exe):
                                found_engines.append(engine)
                                found_exes.append(_exe)
                        # Search system PATH.
                        else:
                            try:
                                exe = _SireBase.findExe(exe).absolute_file_path()
                                found_engines.append(engine)
                                found_exes.append(exe)
                            except:
                                pass
                    # GROMACS
                    elif engine == "GROMACS":
                        if _gmx_exe is not None and _os.path.basename(_gmx_exe) == exe:
                            found_engines.append(engine)
                            found_exes.append(_gmx_exe)
                    # OPENMM
                    elif engine == "OPENMM":
                        found_engines.append(engine)
                        found_exes.append(_SireBase.getBinDir() + "/sire_python")
                    # SOMD
                    elif engine == "SOMD":
                        found_engines.append(engine)
                        if is_free_energy:
                            found_exes.append(_SireBase.getBinDir() + "/somd-freenrg")
                        else:
                            found_exes.append(_SireBase.getBinDir() + "/somd")
                    # Search system PATH.
                    else:
                        try:
                            exe = _SireBase.findExe(exe).absolute_file_path()
                            found_engines.append(engine)
                            found_exes.append(exe)
                        except:
                            pass

    # No engine was found.
    if len(found_engines) == 0:
        if md_engine == "AUTO":
            raise _MissingSoftwareError(
                "Couldn't find an engine that supports the protocol!"
            )
        else:
            raise _MissingSoftwareError(
                "The chosen engine doesn't support the protocol!"
            )

    return found_engines, found_exes


def run(
    system,
    protocol,
    engine="AUTO",
    gpu_support=False,
    auto_start=True,
    name="md",
    work_dir=None,
    seed=None,
    property_map={},
    **kwargs,
):
    """
    Auto-configure and run a molecular dynamics process.

    Parameters
    ----------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The molecular system.

    protocol : :class:`Protocol <BioSimSpace.Protocol>`
        The simulation protocol.

    engine : str
        The molecular dynamics engine to use. If "AUTO", then a matching
        engine will automatically be chosen. Supported engines can be
        found using 'BioSimSpace.MD.engines()'.

    gpu_support : bool
        Whether to choose an engine with GPU support.

    auto_start : bool
        Whether to start the process automatically.

    name : str
        The name of the process.

    work_dir : str
        The working directory for the process.

    seed : int
        A random number seed.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    kwargs : dict
        A dictionary of optional keyword arguments neeeded by the engine.

    Returns
    -------

    process : :class:`Process <BioSimSpace.Process>`
        A process to run the molecular dynamics protocol.
    """

    # Check that the system is valid.
    if not isinstance(system, _System):
        raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")

    # Make sure the system is parameterised.
    if not system._isParameterised():
        raise _IncompatibleError(
            "Cannot execute a Process for this System since it appears "
            "to contain molecules that are not parameterised. Consider "
            "using the 'BioSimSpace.Parameters' engine."
        )

    # Check that the protocol is valid.
    if not isinstance(protocol, _Protocol._protocol.Protocol):
        if isinstance(protocol, str):
            protocol = _Protocol.Custom(protocol)
        else:
            raise TypeError(
                "'protocol' must be of type 'BioSimSpace.Protocol' "
                "or the path to a custom configuration file."
            )

    # Validate optional arguments.

    if not isinstance(engine, str):
        raise TypeError("'engine' must be of type 'str'.")
    md_engine = engine.upper().replace(" ", "")

    if not isinstance(gpu_support, bool):
        raise TypeError("'gpu_support' must be of type 'bool'")

    if not isinstance(auto_start, bool):
        raise TypeError("'auto_start' must be of type 'bool'")

    if not isinstance(name, str):
        raise TypeError("'name' must be of type 'str'")

    if work_dir is not None:
        if not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'")

    if seed is not None:
        if not type(seed) is int:
            raise TypeError("'seed' must be of type 'int'")

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    # Find a molecular dynamics engine and executable.
    engines, exes = _find_md_engines(system, protocol, md_engine, gpu_support)

    # Create the process object, return the first supported engine that can
    # instantiate a process.

    for engine, exe in zip(engines, exes):
        try:
            # AMBER.
            if engine == "AMBER":
                process = _Process.Amber(
                    system,
                    protocol,
                    exe=exe,
                    name=name,
                    work_dir=work_dir,
                    seed=seed,
                    property_map=property_map,
                    **kwargs,
                )

            # GROMACS.
            elif engine == "GROMACS":
                process = _Process.Gromacs(
                    system,
                    protocol,
                    exe=exe,
                    name=name,
                    work_dir=work_dir,
                    seed=seed,
                    property_map=property_map,
                    **kwargs,
                )

            # NAMD.
            elif engine == "NAMD":
                process = _Process.Namd(
                    system,
                    protocol,
                    exe=exe,
                    name=name,
                    work_dir=work_dir,
                    seed=seed,
                    property_map=property_map,
                    **kwargs,
                )

            # OPENMM.
            elif engine == "OPENMM":
                if gpu_support:
                    platform = "CUDA"
                else:
                    platform = "CPU"
                # Don't pass the executable name through so that this works on Windows too.
                process = _Process.OpenMM(
                    system,
                    protocol,
                    exe=None,
                    name=name,
                    work_dir=work_dir,
                    seed=seed,
                    property_map=property_map,
                    platform=platform,
                    **kwargs,
                )

            # SOMD.
            elif engine == "SOMD":
                if gpu_support:
                    platform = "CUDA"
                else:
                    platform = "CPU"
                # Don't pass the executable name through so that this works on Windows too.
                process = _Process.Somd(
                    system,
                    protocol,
                    exe=None,
                    name=name,
                    work_dir=work_dir,
                    seed=seed,
                    property_map=property_map,
                    platform=platform,
                    **kwargs,
                )

            # Start the process.
            if auto_start:
                return process.start()
            else:
                return process

        except:
            pass

    # If we got here, then we couldn't create a process.
    if md_engine == "AUTO":
        raise Exception(
            f"Unable to create a process using any supported engine: {engines}"
        )
    else:
        raise Exception(
            f"Unable to create a process using the chosen engine: {md_engine}"
        )
