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

"""Functionality for relative and absolute alchemical free-energy simulations."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["AlchemicalFreeEnergy", "getData"]


import copy as _copy
import math as _math
import os as _os
import shutil as _shutil
import subprocess as _subprocess
import sys as _sys
import warnings as _warnings
import zipfile as _zipfile
from glob import glob as _glob

from .._Utils import _assert_imported, _have_imported, _try_import

# alchemlyb isn't available on all variants of Python that we support, so we
# need to try_import it.
_alchemlyb = _try_import("alchemlyb")

if _have_imported(_alchemlyb):
    from alchemlyb.workflows import ABFE
    from alchemlyb.postprocessors.units import R_kJmol as _R_kJmol
    from alchemlyb.postprocessors.units import kJ2kcal as _kJ2kcal
    from alchemlyb.preprocessing.subsampling import (
        statistical_inefficiency as _statistical_inefficiency,
    )

    try:
        from alchemlyb.estimators import AutoMBAR as _AutoMBAR
    except ImportError:
        from alchemlyb.estimators import MBAR as _AutoMBAR
    from alchemlyb.estimators import TI as _TI
    from alchemlyb.postprocessors.units import to_kcalmol as _to_kcalmol

import numpy as _np
import pandas as _pd
from sire.legacy.Base import getBinDir as _getBinDir
from sire.legacy.Base import getShareDir as _getShareDir

from ._restraint import Restraint as _Restraint

from .._Exceptions import AnalysisError as _AnalysisError
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import System as _System
from .._Utils import cd as _cd
from .. import _gmx_exe
from .. import _is_notebook
from .. import Process as _Process
from .. import Protocol as _Protocol
from .. import Types as _Types
from .. import Units as _Units
from .. import _Utils

from ..MD._md import _find_md_engines

if _is_notebook:
    from IPython.display import FileLink as _FileLink

# Check that the analyse_freenrg script exists.
if _sys.platform != "win32":
    _analyse_freenrg = _os.path.join(_getBinDir(), "analyse_freenrg")
else:
    _analyse_freenrg = _os.path.join(
        _os.path.normpath(_getShareDir()), "scripts", "analyse_freenrg.py"
    )
if not _os.path.isfile(_analyse_freenrg):
    raise _MissingSoftwareError(
        "Cannot find free energy analysis script in expected location: '%s'"
        % _analyse_freenrg
    )
if _sys.platform == "win32":
    _analyse_freenrg = "%s %s" % (
        _os.path.join(_os.path.normpath(_getBinDir()), "sire_python.exe"),
        _analyse_freenrg,
    )


class AlchemicalFreeEnergy:
    """
    Class for configuring and running relative and absolute
    free-energy perturbation simulations.
    """

    # Create a list of supported molecular dynamics engines.
    _engines = ["AMBER", "GROMACS", "SOMD"]

    def __init__(
        self,
        system,
        protocol=None,
        work_dir=None,
        engine=None,
        gpu_support=False,
        setup_only=False,
        ignore_warnings=False,
        show_errors=True,
        extra_options=None,
        extra_lines=None,
        estimator="MBAR",
        restraint=None,
        property_map={},
    ):
        """
        Constructor.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system for the perturbation leg. This must contain
            a single perturbable molecule and is assumed to have already
            been equilibrated.

        protocol : :class:`Protocol.FreeEnergy <BioSimSpace.Protocol.FreeEnergy>`, \
            The simulation protocol.

        work_dir : str
            The working directory for the free-energy perturbation
            simulation.

        engine : str
            The molecular dynamics engine used to run the simulation. Available
            options are "AMBER", "GROMACS", or "SOMD". If this argument is omitted
            then BioSimSpace will choose an appropriate engine for you.

        gpu_support : bool
            Whether the engine must have GPU support.

        setup_only : bool
            Whether to only support simulation setup. If True, then no
            simulation processes objects will be created, only the directory
            hierarchy and input files to run a simulation externally. This
            can be useful when you don't intend to use BioSimSpace to run
            the simulation. Note that a 'work_dir' must also be specified.

        ignore_warnings : bool
            Whether to ignore warnings when generating the binary run file.
            This option is specific to GROMACS and will be ignored when a
            different molecular dynamics engine is chosen.

        show_errors : bool
            Whether to show warning/error messages when generating the binary
            run file. This option is specific to GROMACS and will be ignored
            when a different molecular dynamics engine is chosen.

        extra_options : dict
            A dictionary containing extra options. Overrides the ones generated from the protocol.

        extra_lines : list
            A list of extra lines to be put at the end of the script.

        estimator : str
            Estimator used for the analysis - must be either 'MBAR' or 'TI'.

        restraint : :class:`Restraint <BioSimSpace.FreeEnergy.Restraint>`
            A Restraint object specifying a receptor-ligand restraint for
            ABFE calculations.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Validate the input.

        if not isinstance(system, _System):
            raise TypeError(
                "'system' must be of type 'BioSimSpace._SireWrappers.System'"
            )
        else:
            # Store a copy of solvated system.
            self._system = system.copy()

        if protocol is not None:
            if isinstance(protocol, _Protocol._FreeEnergyMixin):
                self._protocol = protocol
            else:
                raise TypeError(
                    "'protocol' must be of type 'BioSimSpace.Protocol.FreeEnergy'"
                )
        else:
            # Use a default protocol.
            self._protocol = _Protocol.FreeEnergy()

        self._extra_options = extra_options if extra_options is not None else {}
        self._extra_lines = extra_lines if extra_lines is not None else []

        if not isinstance(setup_only, bool):
            raise TypeError("'setup_only' must be of type 'bool'.")
        else:
            self._setup_only = setup_only

        if work_dir is None and setup_only:
            raise ValueError(
                "A 'work_dir' must be specified when 'setup_only' is True!"
            )

        # Create the working directory.
        self._work_dir = _Utils.WorkDir(work_dir)

        # Validate the user specified molecular dynamics engine.
        self._exe = None
        if engine is not None:
            if not isinstance(engine, str):
                raise TypeError("'engine' must be of type 'str'.")

            # Strip whitespace from engine and convert to upper case.
            engine = engine.replace(" ", "").upper()

            # Check that the engine is supported.
            if engine not in self._engines:
                raise ValueError(
                    f"Unsupported molecular dynamics engine {engine}. "
                    "Supported engines are: %r." % ", ".join(self._engines)
                )

            # Make sure GROMACS is installed if GROMACS engine is selected.
            if engine == "GROMACS":
                if _gmx_exe is None:
                    raise _MissingSoftwareError(
                        "Cannot use GROMACS engine as GROMACS is not installed!"
                    )

                # The system must have a perturbable molecule.
                if (
                    system.nPerturbableMolecules() == 0
                    and system.nDecoupledMolecules() == 0
                ):
                    raise ValueError(
                        "The system must contain a perturbable or decoupled molecule! "
                        "Use the 'BioSimSpace.Align' package to map and merge molecules."
                    )

                pert_type = self._protocol.getPerturbationType()
                if pert_type not in ["full", "release_restraint"]:
                    raise NotImplementedError(
                        "GROMACS currently only supports the 'full' perturbation "
                        "type. Please use engine='SOMD' when running multistep "
                        "perturbation types."
                    )
                if pert_type == "release_restraint":
                    restraint_err = ValueError(
                        "The 'release_restraint' perturbation type requires a multiple "
                        "distance restraint restraint type."
                    )
                    if not restraint:
                        raise restraint_err
                    if restraint._restraint_type != "multiple_distance":
                        raise restraint_err

                self._exe = _gmx_exe
            elif engine == "AMBER":
                # Find a molecular dynamics engine and executable.
                engines, exes = _find_md_engines(system, protocol, engine, gpu_support)
                if not exes:
                    raise _MissingSoftwareError(
                        "Cannot use AMBER engine as AMBER is not installed!"
                    )
                elif len(exes) > 1:
                    _warnings.warn(
                        f"Multiple AMBER engines were found. Proceeding with {exes[0]}..."
                    )
                self._exe = exes[0]

                # The system must have a perturbable molecule.
                if system.nPerturbableMolecules() == 0:
                    raise ValueError(
                        "The system must contain a perturbable molecule! "
                        "Use the 'BioSimSpace.Align' package to map and merge molecules."
                    )

                if self._protocol.getPerturbationType() != "full":
                    raise NotImplementedError(
                        "AMBER currently only supports the 'full' perturbation "
                        "type. Please use engine='SOMD' when running multistep "
                        "perturbation types."
                    )

                # AMBER currently can not handle restraints
                if restraint is not None:
                    raise NotImplementedError(
                        "AMBER currently does not support restraints within BSS"
                    )
        else:
            # Use SOMD as a default.
            engine = "SOMD"

            # The system must have a single perturbable molecule.
            if (
                system.nPerturbableMolecules() != 1
                and system.nDecoupledMolecules() != 1
            ):
                raise ValueError(
                    "The system must contain a single perturbable or decoupled molecule! "
                    "Use the 'BioSimSpace.Align' package to map and merge molecules, or "
                    "to mark them for decoupling."
                )
            # The system must not have a decoupled and perturbable molecule.
            if (
                system.nPerturbableMolecules() == 1
                and system.nDecoupledMolecules() == 1
            ):
                raise ValueError(
                    "The system must not contain a perturbable molecule and a decoupled molecule!"
                )

        # Set the engine.
        self._engine = engine

        if not isinstance(ignore_warnings, bool):
            raise ValueError("'ignore_warnings' must be of type 'bool'.")
        self._ignore_warnings = ignore_warnings

        if not isinstance(show_errors, bool):
            raise ValueError("'show_errors' must be of type 'bool'.")
        self._show_errors = show_errors

        # Check that the estimator is either MBAR or TI.
        if not isinstance(estimator, str):
            raise ValueError("'estimator' must be of type 'str'.")
        if estimator not in ["MBAR", "TI"]:
            raise ValueError("'estimator' must be either 'MBAR' or 'TI'.")
        self._estimator = estimator

        # Check that the map is valid.
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'.")
        self._property_map = property_map

        # Check that if a restraint is passed (bound leg simulation) it is valid.
        # For free leg simulations, the restraint will be None.
        if restraint is not None:
            if not isinstance(restraint, _Restraint):
                raise TypeError(
                    "'restraint' must be of type 'BioSimSpace.FreeEnergy.Restraint'."
                )
            else:
                # Ensure that the system is compatible with the restraint
                restraint.system = self._system

            # Warn the user about instabilities with multiple distance restraints in SOMD.
            if restraint._restraint_type == "multiple_distance" and engine == "SOMD":
                _warnings.warn(
                    "SOMD simulations with some non-interacting ligands can be unstable. This "
                    "affects some systems with multiple distance restraints during the release "
                    "restraint state. If you experience problems, consider using multiple distance "
                    "restraints with GROMACS instead."
                )

        self._restraint = restraint

        # Create fake instance methods for 'analyse' and 'difference'. These
        # pass instance data through to the staticmethod versions.
        self.analyse = self._analyse
        # self.analyse_all_repeats = self._analyse_all_repeats
        self.difference = self._difference

        # Initialise the process runner.
        self._initialise_runner(self._system)

    def run(self, serial=True):
        """
        Run the simulation.

        Parameters
        ----------

        serial : bool
            Whether to run the individual processes for the lambda windows
            in serial.
        """
        if not isinstance(serial, bool):
            raise TypeError("'serial' must be of type 'bool'.")

        if self._setup_only:
            _warnings.warn("No processes exist! Object created in 'setup_only' mode.")
        else:
            self._runner.startAll(serial=serial)

    def wait(self):
        """Wait for the simulation to finish."""
        if self._setup_only:
            _warnings.warn("No processes exist! Object created in 'setup_only' mode.")
        else:
            self._runner.wait()

    def kill(self, index):
        """
        Kill a process for a specific lambda window.

        Parameters
        ----------

        index : int
            The index of the lambda window.
        """
        self._runner.kill(index)

    def killAll(self):
        """Kill any running processes for all lambda windows."""

        self._runner.killAll()

    def workDir(self):
        """
        Return the working directory.

        Returns
        -------

        work_dir : str
            The path of the working directory.
        """
        return str(self._work_dir)

    def getData(self, name="data", file_link=False, work_dir=None):
        """
        Return a link to a zip file containing the data files required for
        post-simulation analysis.

        Parameters
        ----------

        name : str
            The name of the zip file.

        file_link : bool
            Whether to return a FileLink when working in Jupyter.

        work_dir : str
            The working directory for the free-energy perturbation
            simulation.

        Returns
        -------

        output : str, IPython.display.FileLink
            A path, or file link, to an archive of the process input.
        """

        if self._work_dir is None:
            raise ValueError("'work_dir' must be set!")
        else:
            if not isinstance(work_dir, str):
                raise TypeError("'work_dir' must be of type 'str'.")
            if not _os.path.isdir(work_dir):
                raise ValueError("'work_dir' doesn't exist!")

        if not isinstance(name, str):
            raise TypeError("'name' must be of type 'str'")

        # Generate the zip file name.
        zipname = "%s.zip" % name

        # Get the current working directory.
        cwd = _os.getcwd()

        # Change into the working directory.
        with _cd(work_dir):
            # Glob all of the analysis files.

            # First try SOMD data.
            files = _glob("*/*/gradients.dat")

            if len(files) == 0:
                files = _glob("*/*/gromacs.xvg")

                if len(files) == 0:
                    raise ValueError(
                        f"Couldn't find any analysis files in '{work_dir}'"
                    )

            # Write to the zip file.
            with _zipfile.ZipFile(cwd + f"/{zipname}", "w") as zip:
                for file in files:
                    zip.write(file)

        # Return a link to the archive.
        if _is_notebook:
            if file_link:
                # Create a FileLink to the archive.
                f_link = _FileLink(zipname)

                # Set the download attribute so that JupyterLab doesn't try to open the file.
                f_link.html_link_str = (
                    f"<a href='%s' target='_blank' download='{zipname}'>%s</a>"
                )

                # Return a link to the archive.
                return f_link
            else:
                return zipname
        # Return the path to the archive.
        else:
            return zipname

    @staticmethod
    def analyse(work_dir, temperature=None, estimator="MBAR", **kwargs):
        """
        Analyse existing free-energy data from a simulation working directory.

        Parameters
        ----------

        work_dir : str
            The working directory for the simulation.

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The temperature.

        estimator : str
            The estimator ('MBAR' or 'TI') used. Default is MBAR.

        kwargs :
            For non SOMD engines, keyword arguments passed to
            alchemlyb.workflows.ABFE.

        Returns
        -------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF). The data is a list of tuples,
            where each tuple contains the lambda value, the PMF, and the
            standard error.

        overlap : [ [ float, float, ... ] ]
            The overlap matrix. This gives the overlap between each lambda
            window.
        """

        _assert_imported(_alchemlyb)

        if not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError("'work_dir' doesn't exist!")

        if estimator not in ["MBAR", "TI"]:
            raise ValueError("'estimator' must be either 'MBAR' or 'TI'.")

        mask_dict = {
            "SOMD": "/lambda_*/gradients.dat",
            "GROMACS": "/lambda_*/gromacs.xvg",
            "AMBER": "/lambda_*/amber.out",
            "PARQUET": "/lambda_*/*.parquet",
        }

        for engine, mask in mask_dict.items():
            data = _glob(work_dir + mask)
            if data:
                if engine == "SOMD":
                    return AlchemicalFreeEnergy._analyse_somd(work_dir, estimator)
                else:
                    if not isinstance(temperature, _Types.Temperature):
                        raise TypeError(
                            "'temperature' must be of type 'BioSimSpace.Types.Temperature'"
                        )
                    return AlchemicalFreeEnergy._analyse_noSOMD(
                        engine=engine,
                        work_dir=work_dir,
                        estimator=estimator,
                        temperature=temperature,
                        **kwargs,
                    )

        raise ValueError(
            "Couldn't find any SOMD, GROMACS, AMBER or PARQUET free-energy output?"
        )

    def _analyse(self):
        """
        Analyse free-energy data for this object using MBAR.

        Returns
        -------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF). The data is a list of tuples,
            where each tuple contains the lambda value, the PMF, and the
            standard error.

        overlap : [ [ float, float, ... ] ]
            The overlap matrix. This gives the overlap between each lambda
            window.
        """

        # Return the result of calling the staticmethod, passing in the working
        # directory and estimator of this object.
        if isinstance(self._protocol, _Protocol.Production):
            temperature = self._protocol.getTemperature()
        else:
            raise TypeError(
                "'protocol' must be of type 'BioSimSpace.Protocol.Production'"
            )
        return AlchemicalFreeEnergy.analyse(
            str(self._work_dir), self._estimator, temperature=temperature
        )

    @staticmethod
    def _analyse_noSOMD(
        engine=None, work_dir=None, estimator=None, temperature=None, **kwargs
    ):
        """
        Analyse the free energy data, that is not from SOMD (e.g. GROMACS,
        AMBER, NAMD, GOMC).

        Parameters
        ----------

        engine : str
            The molecular dynamics engine used to run the simulation. Available
            options are "AMBER", "GROMACS".

        work_dir : str
            The path to the working directory.

        estimator : str
            The estimator ('MBAR' or 'TI') used.

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The temperature.

        kwargs :
            Keyword arguments passed to alchemlyb.workflows.ABFE.

        Returns
        -------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF). The data is a list of tuples,
            where each tuple contains the lambda value, the PMF, and the
            standard error.

        overlap : [ [ float, float, ... ] ]
            The overlap matrix. This gives the overlap between each lambda
            window. For TI, this gives the dhdl.
        """

        if not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError("'work_dir' doesn't exist!")

        if estimator in ["MBAR", "BAR"]:
            prefix = "u_nk"
        elif estimator == "TI":
            prefix = "dHdl"
        else:
            raise ValueError("'estimator' must be either 'MBAR' or 'TI'.")

        try:
            # Use the parquet files if they are available.
            workflow = ABFE(
                units="kcal/mol",
                software="PARQUET",
                dir=work_dir,
                prefix="/lambda_*/" + prefix,
                suffix="parquet",
                T=temperature / _Units.Temperature.kelvin,
                outdirectory=work_dir,
                **kwargs,
            )
            workflow.run(estimators=estimator, breakdown=None, forwrev=None, **kwargs)
        except ValueError:
            if engine == "AMBER":
                prefix = "amber"
                suffix = "out"
            elif engine == "GROMACS":
                prefix = "gromacs"
                suffix = "xvg"
            else:
                raise ValueError(f"{engine} has to be either 'AMBER' or " f"'GROMACS'.")

            workflow = ABFE(
                units="kcal/mol",
                software=engine,
                dir=work_dir,
                prefix="lambda_*/" + prefix,
                suffix=suffix,
                T=temperature / _Units.Temperature.kelvin,
                outdirectory=work_dir,
                **kwargs,
            )
            workflow.run(estimators=estimator, breakdown=None, forwrev=None, **kwargs)

        # Extract the data from the mbar results.
        data = []
        # convert the data frames to kcal/mol
        delta_f_ = _to_kcalmol(workflow.estimator[estimator].delta_f_)
        d_delta_f_ = _to_kcalmol(workflow.estimator[estimator].d_delta_f_)
        for lambda_, mbar_value, mbar_error in zip(
            delta_f_.index, delta_f_.iloc[0, :], d_delta_f_.iloc[0, :]
        ):
            # Append the data.
            data.append(
                (
                    lambda_,
                    (mbar_value) * _Units.Energy.kcal_per_mol,
                    (mbar_error) * _Units.Energy.kcal_per_mol,
                )
            )
        if estimator == "MBAR":
            return data, workflow.estimator[estimator].overlap_matrix
        elif estimator == "TI":
            # In the original implementation when TI is used,
            # The TI estimator object is returned as overlap matrix, I
            # don't know why this is the case but I will preserve this
            # behaviour.
            return data, workflow.estimator[estimator]

    @staticmethod
    def _analyse_somd(work_dir=None, estimator=None):
        """
        Analyse the SOMD free energy data.

        Parameters
        ----------

        work_dir : str
            The path to the working directory.

        estimator : str
            The estimator ('MBAR' or 'TI') used.

        Returns
        -------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF). The data is a list of tuples,
            where each tuple contains the lambda value, the PMF, and the
            standard error.

        overlap : [ [ float, float, ... ] ]
            The overlap matrix. This gives the overlap between each lambda
            window. For TI, this gives the dhdl.
        """

        if not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'.")
        if not _os.path.isdir(work_dir):
            raise ValueError("'work_dir' doesn't exist!")

        if estimator not in ["MBAR", "TI"]:
            raise ValueError(
                "'estimator' must be either 'MBAR' or 'TI' for SOMD output."
            )

        files = sorted(_glob(work_dir + "/lambda_*/simfile.dat"))
        if not files:
            raise FileNotFoundError(
                f"No simfile.dat files found in lambda_* directories in {work_dir}"
            )
        lambdas = [float(x.split("/")[-2].split("_")[-1]) for x in files]

        temperatures = []
        for file in files:
            found_temperature = False
            with open(file, "r") as f:
                for line in f.readlines():
                    temp = None
                    start = "#Generating temperature is"
                    if start in line:
                        split_line = line.split()
                        temp = split_line[3]
                        try:
                            unit = split_line[4]
                        except IndexError:
                            # Must be °C
                            temp, unit = temp.split("°")
                        if unit == "C":
                            temp = float(temp) + 273.15  # Convert to K
                        else:
                            temp = float(temp)
                        temperatures.append(temp)
                        if temp is not None:
                            found_temperature = True
                            break

            if not found_temperature:
                raise ValueError(
                    f"The temperature was not detected in the SOMD output file, {file}"
                )

        if temperatures[0] != temperatures[-1]:
            raise ValueError("The temperatures at the endstates don't match!")

        def _somd_extract_u_nk(simfile, T):
            """
            Return reduced potentials `u_nk` from Somd outputfile.

            Parameters
            ----------
            outfile : str
                Path to simfile.dat file to extract data from.
            T : float
                Temperature in Kelvin at which the simulations were performed;
                needed to generated the reduced potential (in units of kT)

            Returns
            -------
            u_nk : DataFrame
                Reduced potential for each alchemical state (k) for each frame (n).
            """
            # open the file - check if it is okay, if not raise an error
            file = simfile

            # # beta vale for calc kT - not needed as already reduced in simfile.dat
            # T = T
            # k_b = R_kJmol # * kJ2kcal ??
            # beta = 1/(k_b * T)

            # find out which lambda window
            found_lambda = False
            found_array = False
            found_time = False
            with open(file, "r") as f:
                lambda_win = None
                lambda_array = None
                sim_length = None
                for line in f.readlines():
                    start_w = "#Generating lambda is"
                    start_a = "#Alchemical array is"
                    start_t = " and "
                    end_t = " ps"
                    if start_w in line:
                        lambda_win = float(line.replace(start_w, "").strip())
                        if lambda_win is not None:
                            found_lambda = True
                    if start_a in line:
                        lambda_array = (
                            (line.replace(start_a, ""))
                            .strip()
                            .replace("(", "")
                            .replace(")", "")
                            .replace(" ", "")
                        ).split(
                            ","
                        )  # list
                        if lambda_array is not None:
                            found_array = True
                    if start_t and end_t in line:
                        sim_length = float(
                            ((line.split(start_t)[1]).split(end_t)[0]).strip()
                        )
                        if sim_length is not None:
                            found_time = True
                    if found_lambda:
                        if found_array:
                            if found_time:
                                break

            if not found_lambda:
                raise ValueError(
                    f"The lambda window was not detected in the SOMD output file, {file}"
                )

            if not found_array:
                raise ValueError(
                    f"The lambda array was not detected in the SOMD output file, {file}"
                )

            if not found_time:
                raise ValueError(
                    f"The simulation time was not detected in the SOMD output file, {file}"
                )

            # get header from things instead of like this
            header = [
                "step",
                "potential_kcal/mol",
                "gradient_kcal/mol",
                "forward_Metropolis",
                "backward_Metropolis",
            ]
            header.extend(lambda_array)

            file_df = _pd.read_fwf(
                file, skipinitialspace=True, skiprows=13, header=None, names=header
            )
            # print(file_df)

            time_step = sim_length / len(file_df["step"])
            time_rows = _np.arange(0, len(file_df["step"]), 1)
            time = _np.arange(0, sim_length, time_step)

            mbar_energies = (
                []
            )  # results in list of lists where each list is 0 to 1 window values

            # # so for the energies for each lambda, append the kt to the data list of vals for all lambda wins
            # then trun into df
            for t in time_rows:
                row = file_df.loc[t][lambda_array].to_numpy()
                E_ref = row[lambda_array.index(str(lambda_win))]
                energies = []
                for lam in lambda_array:
                    E_ = row[lambda_array.index(lam)]
                    energies.append((E_ - E_ref))
                mbar_energies.append(energies)

            df = _pd.DataFrame(
                mbar_energies,
                columns=_np.array(lambda_array, dtype=_np.float64),
                index=_pd.MultiIndex.from_arrays(
                    [time, _np.repeat(lambda_win, len(time))], names=["time", "lambdas"]
                ),
            )
            df.attrs["temperature"] = T
            df.attrs["energy_unit"] = "kT"

            return df

        def _somd_extract_dHdl(simfile, T):
            """
            Return gradients ``dH/dl`` from Somd outputfile.

            Parameters
            ----------
            outfile : str
                Path to simfile.dat file to extract data from.
            T : float
                Temperature in Kelvin at which the simulations were performed.

            Returns
            -------
            dH/dl : Series
                dH/dl as a function of time for this lambda window.
            """
            # open the file
            file = simfile

            # for dhdl need to consider the T, as the gradient is in kcal/mol in the simfile.dat
            T = 300
            k_b = _R_kJmol * _kJ2kcal
            beta = 1 / (k_b * T)

            found_lambda = False
            found_array = False
            found_time = False
            with open(file, "r") as f:
                lambda_win = None
                lambda_array = None
                sim_length = None
                for line in f.readlines():
                    start_w = "#Generating lambda is"
                    start_a = "#Alchemical array is"
                    start_t = " and "
                    end_t = " ps"
                    if start_w in line:
                        lambda_win = float(line.replace(start_w, "").strip())
                        if lambda_win is not None:
                            found_lambda = True
                    if start_a in line:
                        lambda_array = (
                            (line.replace(start_a, ""))
                            .strip()
                            .replace("(", "")
                            .replace(")", "")
                            .replace(" ", "")
                        ).split(
                            ","
                        )  # list
                        if lambda_array is not None:
                            found_array = True
                    if start_t and end_t in line:
                        sim_length = float(
                            ((line.split(start_t)[1]).split(end_t)[0]).strip()
                        )
                        if sim_length is not None:
                            found_time = True
                    if found_lambda:
                        if found_array:
                            if found_time:
                                break

            if not found_lambda:
                raise ValueError(
                    f"The lambda window was not detected in the SOMD output file, {file}"
                )

            if not found_array:
                raise ValueError(
                    f"The lambda array was not detected in the SOMD output file, {file}"
                )

            if not found_time:
                raise ValueError(
                    f"The simulation time was not detected in the SOMD output file, {file}"
                )

            # get header
            header = [
                "step",
                "potential_kcal/mol",
                "gradient_kcal/mol",
                "forward_Metropolis",
                "backward_Metropolis",
            ]
            header.extend(lambda_array)

            file_df = _pd.read_fwf(
                file, skipinitialspace=True, skiprows=13, header=None, names=header
            )

            time_step = sim_length / len(file_df["step"])
            time_rows = _np.arange(0, len(file_df["step"]), 1)
            time = _np.arange(0, sim_length, time_step)

            gradient_energies = []  # results in list of the gradients at that lambda

            # turn gradient into list of reduced gradients
            for t in time_rows:
                gradient = file_df.loc[t]["gradient_kcal/mol"]
                red_gradient = gradient * beta
                gradient_energies.append(red_gradient)

            # df in the format needed for alchemlyb
            df = _pd.DataFrame(
                gradient_energies,
                columns=["fep"],
                index=_pd.MultiIndex.from_arrays(
                    [time, _np.repeat(lambda_win, len(time))],
                    names=["time", "fep-lambda"],
                ),
            )

            df.attrs["temperature"] = T
            df.attrs["energy_unit"] = "kT"

            return df

        # Process the data files using the function defined above.
        if estimator == "MBAR":
            # Process the data files using the alchemlyb library.
            # Subsample according to statistical inefficiency and then calculate the MBAR.
            sample_okay = False
            try:
                u_nk = [_somd_extract_u_nk(x, T=t) for x, t in zip(files, temperatures)]
                sampled_u_nk = _alchemlyb.concat(
                    [_statistical_inefficiency(i, i.iloc[:, 0]) for i in u_nk]
                )
                mbar = _AutoMBAR().fit(sampled_u_nk)
                sample_okay = True
            except:
                print("Could not calculate statistical inefficiency.")

            if not sample_okay:
                print(
                    "Running without calculating the statistical inefficiency and without subsampling..."
                )
                try:
                    u_nk = _alchemlyb.concat(
                        [
                            _somd_extract_u_nk(x, T=t)
                            for x, t in zip(files, temperatures)
                        ]
                    )
                    mbar = _AutoMBAR().fit(u_nk)
                except:
                    raise _AnalysisError("MBAR free-energy analysis failed!")

            # Extract the data from the mbar results.
            data = []
            # convert the data frames to kcal/mol
            delta_f_ = _to_kcalmol(mbar.delta_f_)
            d_delta_f_ = _to_kcalmol(mbar.d_delta_f_)
            for lambda_, t in zip(lambdas, temperatures):
                x = lambdas.index(lambda_)
                mbar_value = delta_f_.iloc[0, x]
                mbar_error = d_delta_f_.iloc[0, x]

                # Append the data.
                data.append(
                    (
                        lambda_,
                        (mbar_value) * _Units.Energy.kcal_per_mol,
                        (mbar_error) * _Units.Energy.kcal_per_mol,
                    )
                )

            # Calculate overlap matrix.
            overlap = mbar.overlap_matrix

            return (data, overlap)

        if estimator == "TI":
            # Process the data files using the alchemlyb library.
            # Subsample according to statistical inefficiency and then calculate the TI.
            sample_okay = False
            try:
                dhdl = [_somd_extract_dHdl(x, T=t) for x, t in zip(files, temperatures)]
                sampled_dhdl = _alchemlyb.concat(
                    [_statistical_inefficiency(i, i.iloc[:, 0]) for i in dhdl]
                )
                ti = _TI().fit(sampled_dhdl)
                sample_okay = True
            except:
                print("Could not calculate statistical inefficiency.")

            if not sample_okay:
                print(
                    "Running without calculating the statistical inefficiency and without subsampling..."
                )
                try:
                    dhdl = _alchemlyb.concat(
                        [
                            _somd_extract_dHdl(x, T=t)
                            for x, t in zip(files, temperatures)
                        ]
                    )
                    ti = _TI().fit(dhdl)
                except:
                    raise _AnalysisError("TI free-energy analysis failed!")

            # Extract the data from the ti results.
            data = []
            # convert the data frames to kcal/mol
            delta_f_ = _to_kcalmol(ti.delta_f_)
            d_delta_f_ = _to_kcalmol(ti.d_delta_f_)
            for lambda_ in lambdas:
                x = lambdas.index(lambda_)
                ti_value = delta_f_.iloc[0, x]
                ti_error = d_delta_f_.iloc[1, x]

                # Append the data.
                data.append(
                    (
                        lambda_,
                        (ti_value) * _Units.Energy.kcal_per_mol,
                        (ti_error) * _Units.Energy.kcal_per_mol,
                    )
                )

            # For TI, dHdl graph.
            overlap = ti

            return (data, overlap)

    @staticmethod
    def difference(pmf, pmf_ref):
        """
        Compute the relative free-energy difference between two perturbation
        legs.

        Parameters
        ----------

        pmf : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The potential of mean force (PMF) of interest. The data is a list
            of tuples, where each tuple contains the lambda value, the PMF,
            and the standard error.

        pmf_ref : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The reference potential of mean force (PMF). The data is a list
            of tuples, where each tuple contains the lambda value, the PMF,
            and the standard error.

        Returns
        -------

        free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
            The relative free-energy difference and its associated error.
        """

        if not isinstance(pmf, list):
            raise TypeError("'pmf' must be of type 'list'.")
        if not isinstance(pmf_ref, list):
            raise TypeError("'pmf_ref' must be of type 'list'.")

        for rec in pmf:
            if not isinstance(rec, tuple):
                raise TypeError(
                    "'pmf1' must contain tuples containing lambda "
                    "values and the associated free-energy and error."
                )
            else:
                if len(rec) != 3:
                    raise ValueError(
                        "Each tuple in 'pmf1' must contain three items: "
                        "a lambda value and the associated free energy "
                        "and error."
                    )
                for val in rec[1:]:
                    if not isinstance(val, _Types.Energy):
                        raise TypeError(
                            "'pmf' must contain 'BioSimSpace.Types.Energy' types."
                        )

        for rec in pmf_ref:
            if not isinstance(rec, tuple):
                raise TypeError(
                    "'pmf_ref' must contain tuples containing lambda "
                    "values and the associated free-energy and error."
                )
            else:
                if len(rec) != 3:
                    raise ValueError(
                        "Each tuple in 'pmf_ref' must contain three items: "
                        "a lambda value and the associated free energy "
                        "and error."
                    )
                for val in rec[1:]:
                    if not isinstance(val, _Types.Energy):
                        raise TypeError(
                            "'pmf_ref' must contain 'BioSimSpace.Types.Energy' types."
                        )

        # Work out the difference in free energy.
        free_energy = (pmf[-1][1] - pmf[0][1]) - (pmf_ref[-1][1] - pmf_ref[0][1])

        # Propagate the errors. (These add in quadrature.)

        # Measure.
        error0 = _math.sqrt(
            (pmf[-1][2].value() * pmf[-1][2].value())
            + (pmf[0][2].value() * pmf[0][2].value())
        )

        # Reference.
        error1 = _math.sqrt(
            (pmf_ref[-1][2].value() * pmf_ref[-1][2].value())
            + (pmf_ref[0][2].value() * pmf_ref[0][2].value())
        )

        # Error for free-energy difference.
        error = (
            _math.sqrt((error0 * error0) + (error1 * error1))
            * _Units.Energy.kcal_per_mol
        )

        return (free_energy, error)

    def _difference(self, pmf_ref):
        """
        Compute the relative free-energy difference between two perturbation
        legs.

        Parameters
        ----------

        pmf_ref : [(float, :class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)]
            The reference potential of mean force (PMF). The data is a list
            of tuples, where each tuple contains the lambda value, the PMF,
            and the standard error.

        Returns
        -------

        free_energy : (:class:`Energy <BioSimSpace.Types.Energy>`, :class:`Energy <BioSimSpace.Types.Energy>`)
            The relative free-energy difference and its associated error.
        """

        # Calculate the PMF for this object.
        pmf, _ = self.analyse()

        # Now call the staticmethod passing in both PMFs.
        return AlchemicalFreeEnergy.difference(pmf, pmf_ref)

    def _initialise_runner(self, system):
        """
        Internal helper function to initialise the process runner.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.
        """

        # Initialise list to store the processes
        processes = []

        # Convert to an appropriate water topology.
        if self._engine in ["AMBER", "SOMD"]:
            system._set_water_topology("AMBER", property_map=self._property_map)
        elif self._engine == "GROMACS":
            system._set_water_topology("GROMACS", property_map=self._property_map)

        # Setup all of the simulation processes for each leg.

        # Get the lambda values from the protocol for the first leg.
        lam_vals = self._protocol.getLambdaValues(type="dataframe")

        # Create a process for the first lambda value.
        lam = lam_vals.loc[0]

        # Update the protocol lambda values.
        self._protocol.setLambdaValues(lam=lam, lam_vals=lam_vals)

        # Create and append the required processes for each leg.
        # Nest the working directories inside self._work_dir.

        # Name the first directory.
        if self._engine == "SOMD":
            first_dir = f"{self._work_dir}/lambda_{lam.values[0]:5.4f}"
        else:
            first_dir = f"{self._work_dir}/lambda_{self._protocol.getLambdaIndex()}"

        # SOMD.
        if self._engine == "SOMD":
            # Check for GPU support.
            if "CUDA_VISIBLE_DEVICES" in _os.environ:
                platform = "CUDA"
            else:
                platform = "CPU"

            # Check against passing multiple sets of lam vals to SOMD
            if lam_vals.shape[1] != 1:
                raise ValueError(
                    "SOMD can only handle a single set of lambda values for a given perturbation."
                )

            first_process = _Process.Somd(
                system,
                self._protocol,
                platform=platform,
                work_dir=first_dir,
                property_map=self._property_map,
                extra_options=self._extra_options,
                extra_lines=self._extra_lines,
                restraint=self._restraint,
            )

        # GROMACS.
        elif self._engine == "GROMACS":
            first_process = _Process.Gromacs(
                system,
                self._protocol,
                work_dir=first_dir,
                ignore_warnings=self._ignore_warnings,
                show_errors=self._show_errors,
                extra_options=self._extra_options,
                extra_lines=self._extra_lines,
                restraint=self._restraint,
            )

        # AMBER.
        # TODO: Make the restraint valid for AMBER
        elif self._engine == "AMBER":
            first_process = _Process.Amber(
                system,
                self._protocol,
                exe=self._exe,
                work_dir=first_dir,
                extra_options=self._extra_options,
                extra_lines=self._extra_lines,
            )

        if self._setup_only and not self._engine == "GROMACS":
            del first_process
        else:
            processes.append(first_process)

        # Loop over the rest of the lambda values.
        for x, lam in lam_vals.iterrows():
            if x == 0:
                # Skip the zero-th window
                continue
            self._protocol.setLambda(lam)
            # Name the directory.
            if self._engine == "SOMD":
                new_dir = f"{self._work_dir}/lambda_{lam.values[0]:5.4f}"
            # If using GROMACS, doesn't make sense to use single value of lambda
            # when using lambda array.
            else:
                new_dir = f"{self._work_dir}/lambda_{self._protocol.getLambdaIndex()}"

            # Use absolute path.
            if not _os.path.isabs(new_dir):
                new_dir = _os.path.abspath(new_dir)

            # Delete any existing directories.
            if _os.path.isdir(new_dir):
                _shutil.rmtree(new_dir, ignore_errors=True)

            # Copy the first directory to that of the current lambda value.
            _shutil.copytree(first_dir, new_dir)

            # Now update the lambda values in the config files.

            # SOMD.
            if self._engine == "SOMD":
                new_config = []
                with open(new_dir + "/somd.cfg", "r") as f:
                    for line in f:
                        if "lambda_val" in line:
                            # Get lam val from the Series
                            new_config.append(f"lambda_val = {lam.values[0]}\n")
                        else:
                            new_config.append(line)
                with open(new_dir + "/somd.cfg", "w") as f:
                    for line in new_config:
                        f.write(line)

                # Create a copy of the process and update the working
                # directory.
                if not self._setup_only:
                    process = _copy.copy(first_process)
                    process._system = first_process._system.copy()
                    process._protocol = self._protocol
                    process._work_dir = new_dir
                    process._std_out_file = new_dir + "/somd.out"
                    process._std_err_file = new_dir + "/somd.err"
                    process._rst_file = new_dir + "/somd.rst7"
                    process._top_file = new_dir + "/somd.prm7"
                    process._traj_file = new_dir + "/traj000000001.dcd"
                    process._restart_file = new_dir + "/latest.rst"
                    process._config_file = new_dir + "/somd.cfg"
                    process._pert_file = new_dir + "/somd.pert"
                    process._gradients_file = new_dir + "/gradients.dat"
                    process._input_files = [
                        process._config_file,
                        process._rst_file,
                        process._top_file,
                        process._pert_file,
                    ]
                    processes.append(process)

            # GROMACS.
            elif self._engine == "GROMACS":
                # Delete the existing binary run file.
                _os.remove(new_dir + "/gromacs.tpr")
                new_config = []
                with open(new_dir + "/gromacs.mdp", "r") as f:
                    for line in f:
                        if "init-lambda-state" in line:
                            new_config.append(
                                "init-lambda-state = %d\n"
                                % (self._protocol.getLambdaIndex())
                            )
                        else:
                            new_config.append(line)
                with open(new_dir + "/gromacs.mdp", "w") as f:
                    for line in new_config:
                        f.write(line)

                mdp = new_dir + "/gromacs.mdp"
                gro = new_dir + "/gromacs.gro"
                top = new_dir + "/gromacs.top"
                tpr = new_dir + "/gromacs.tpr"

                # Use grompp to generate the portable binary run input file.
                _Process.Gromacs._generate_binary_run_file(
                    mdp,
                    gro,
                    top,
                    gro,
                    tpr,
                    first_process._exe,
                    ignore_warnings=self._ignore_warnings,
                    show_errors=self._show_errors,
                )

                # Create a copy of the process and update the working
                # directory.
                if not self._setup_only:
                    process = _copy.copy(first_process)
                    process._system = first_process._system.copy()
                    process._protocol = self._protocol
                    process._work_dir = new_dir
                    process._std_out_file = new_dir + "/gromacs.out"
                    process._std_err_file = new_dir + "/gromacs.err"
                    process._gro_file = new_dir + "/gromacs.gro"
                    process._top_file = new_dir + "/gromacs.top"
                    process._traj_file = new_dir + "/gromacs.trr"
                    process._config_file = new_dir + "/gromacs.mdp"
                    process._tpr_file = new_dir + "/gromacs.tpr"
                    process._input_files = [
                        process._config_file,
                        process._gro_file,
                        process._top_file,
                        process._tpr_file,
                    ]
                    processes.append(process)

            # AMBER.
            elif self._engine == "AMBER":
                new_config = []
                with open(new_dir + "/amber.cfg", "r") as f:
                    for line in f:
                        if "clambda" in line:
                            new_config.append("   clambda=%s,\n" % lam)
                        else:
                            new_config.append(line)
                with open(new_dir + "/amber.cfg", "w") as f:
                    for line in new_config:
                        f.write(line)

                # Create a copy of the process and update the working
                # directory.
                if not self._setup_only:
                    process = _copy.copy(first_process)
                    process._system = first_process._system.copy()
                    process._protocol = self._protocol
                    process._work_dir = new_dir
                    process._std_out_file = new_dir + "/amber.out"
                    process._std_err_file = new_dir + "/amber.err"
                    process._rst_file = new_dir + "/amber.rst7"
                    process._top_file = new_dir + "/amber.prm7"
                    process._traj_file = new_dir + "/amber.nc"
                    process._config_file = new_dir + "/amber.cfg"
                    process._nrg_file = new_dir + "/amber.nrg"
                    process._input_files = [
                        process._config_file,
                        process._rst_file,
                        process._top_file,
                    ]
                    processes.append(process)

        if not self._setup_only:
            # Initialise the process runner. All processes have already been nested
            # inside the working directory so no need to re-nest.
            self._runner = _Process.ProcessRunner(processes)

    def _update_run_args(self, args):
        """
        Internal function to update run arguments for all subprocesses.

        Parameters
        ----------

        args : dict, collections.OrderedDict
            A dictionary which contains the new command-line arguments
            for the process executable.
        """

        if not isinstance(args, dict):
            raise TypeError("'args' must be of type 'dict'")

        for process in self._runner.processes():
            process.setArgs(args)


def getData(name="data", file_link=False, work_dir=None):
    """
    Return a link to a zip file containing the data files required for
    post-simulation analysis.

    Parameters
    ----------

    name : str
        The name of the zip file.

    file_link : bool
        Whether to return a FileLink when working in Jupyter.

    work_dir : str
        The working directory for the simulation.

    Returns
    -------

    output : str, IPython.display.FileLink
        A path, or file link, to an archive of the process input.
    """
    return AlchemicalFreeEnergy.getData(
        name=name, file_link=file_link, work_dir=work_dir
    )
