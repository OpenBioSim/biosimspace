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

"""Functionality for running simulations using AMBER."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Amber"]

from .._Utils import _try_import

_pygtail = _try_import("pygtail")


from . import _process


class Amber(_process.Process):
    """A class for running simulations using AMBER."""

    def __init__(
        self,
        system,
        protocol,
        reference_system=None,
        explicit_dummies=False,
        exe=None,
        is_gpu=False,
        name="amber",
        work_dir=None,
        seed=None,
        extra_options={},
        extra_lines=[],
        extra_args={},
        property_map={},
        **kwargs,
    ):
        """
        Constructor.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.

        protocol : :class:`Protocol <BioSimSpace.Protocol>`
            The protocol for the AMBER process.

        reference_system : :class:`System <BioSimSpace._SireWrappers.System>` or None
            An optional system to use as a source of reference coordinates for position
            restraints. It is assumed that this system has the same topology as "system".
            If this is None, then "system" is used as a reference.

        explicit_dummies : bool
            Whether to keep dummy atoms explicit at alchemical end states, or remove them.
            This option is provided for legacy support of alchemical free energy calculations
            using the old PMEMD CPU implementation. The default is False, which should be
            used for any recent AMBER version or for GPU accelerated PMEMD.

        exe : str
            The full path to the AMBER executable.

        is_gpu : bool
            Whether to use the GPU accelerated version of AMBER.

        name : str
            The name of the process.

        work_dir :
            The working directory for the process.

        seed : int
            A random number seed.

        extra_options : dict
            A dictionary containing extra options. Overrides the defaults generated
            by the protocol.

        extra_lines : [str]
            A list of extra lines to put at the end of the configuration file.

        extra_args : dict
            A dictionary of extra command-line arguments to pass to the AMBER executable.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        kwargs : dict
            Additional keyword arguments.
        """
        import os as _os
        from ..Protocol._free_energy_mixin import _FreeEnergyMixin
        from .._Config import Amber as _AmberConfig

        # Call the base class constructor.
        super().__init__(
            system,
            protocol,
            reference_system=reference_system,
            name=name,
            work_dir=work_dir,
            seed=seed,
            extra_options=extra_options,
            extra_lines=extra_lines,
            extra_args=extra_args,
            property_map=property_map,
        )

        # Set the package name.
        self._package_name = "AMBER"

        # This process can generate trajectory data.
        self._has_trajectory = True

        if not isinstance(is_gpu, bool):
            raise TypeError("'is_gpu' must be of type 'bool'")

        # Check whether this is a vacuum simulation.
        self._is_vacuum = not (
            _AmberConfig.hasBox(self._system, self._property_map)
            or _AmberConfig.hasWater(self._system)
        )

        # Flag to indicate whether the original system has a box.
        self._has_box = _AmberConfig.hasBox(self._system, self._property_map)

        # Take note of whether the original reference system was None
        # This will be used later to avoid duplication
        if reference_system is not None:
            self._is_real_reference = True
        else:
            self._is_real_reference = False

        # If the path to the executable wasn't specified, then search
        # for it in AMBERHOME and the PATH.
        if exe is None:
            if isinstance(protocol, _FreeEnergyMixin):
                is_free_energy = True
            else:
                is_free_energy = False

            self._exe = _findExe(
                is_gpu=is_gpu, is_free_energy=is_free_energy, is_vacuum=self._is_vacuum
            )
        else:
            # Make sure executable exists.
            if _os.path.isfile(exe):
                self._exe = exe
            else:
                raise IOError("AMBER executable doesn't exist: '%s'" % exe)

        # Is this a CUDA enabled version of AMBER?
        if "cuda" in self._exe.lower():
            self._is_pmemd_cuda = True
            self._is_pmemd = False
        else:
            self._is_pmemd_cuda = False
            if "pmemd" in self._exe.lower():
                self._is_pmemd = True
            else:
                self._is_pmemd = False

        if not isinstance(explicit_dummies, bool):
            raise TypeError("'explicit_dummies' must be of type 'bool'")
        self._explicit_dummies = explicit_dummies

        # Initialise the energy dictionary and header.
        self._stdout_dict = _process._MultiDict()

        # Initialise dictionaries to hold stdout records for all possible
        # regions. For regular simulations there will be one, for free-energy
        # simulations there can be up to four, i.e. one for each of the TI regions
        # and one for the soft-core part of the system in each region, if present.
        # The order of the dictionaries is:
        #  - TI region 1
        #  - TI region 1 (soft-core part)
        #  - TI region 2
        #  - TI region 2 (soft-core part)
        self._stdout_dict = [
            _process._MultiDict(),
            _process._MultiDict(),
            _process._MultiDict(),
            _process._MultiDict(),
        ]

        # Initialise mappings between "universal" stdout keys, and the actual
        # record key used for the different regions (and soft-core parts) from
        # in the AMBER output. Ordering is the same as for the stdout_dicts above.
        self._stdout_key = [{}, {}, {}, {}]

        # Flag for the current record region in the AMBER output file.
        self._current_region = 0

        # Initialise log file parsing flags.
        self._has_results = False
        self._finished_results = False
        self._is_header = False

        # The names of the input files.
        self._rst_file = _os.path.join(str(self._work_dir), f"{name}.rst7")
        self._top_file = _os.path.join(str(self._work_dir), f"{name}.prm7")
        self._ref_file = _os.path.join(str(self._work_dir), f"{name}_ref.rst7")

        # The name of the trajectory file.
        self._traj_file = _os.path.join(str(self._work_dir), f"{name}.nc")

        # Set the path for the AMBER configuration file.
        self._config_file = _os.path.join(str(self._work_dir), f"{name}.cfg")

        # Create the list of input files.
        self._input_files = [self._config_file, self._rst_file, self._top_file]

        # Add the reference file if there are position restraints.
        if self._protocol.getRestraint() is not None:
            self._input_files.append(self._ref_file)

        # Now set up the working directory for the process.
        self._setup(**kwargs)

    def _setup(self, **kwargs):
        """Setup the input files and working directory ready for simulation."""
        from .. import IO as _IO
        import shutil as _shutil
        from ..Align._squash import _squash
        import os as _os
        from .. import _isVerbose
        from .. import Protocol as _Protocol
        from ..Protocol._free_energy_mixin import _FreeEnergyMixin

        # Create the input files...

        # Create a copy of the system.
        system = self._system.copy()
        reference_system = self._reference_system.copy()

        # Convert the water model topology so that it matches the AMBER naming convention.
        system._set_water_topology("AMBER", property_map=self._property_map)
        reference_system._set_water_topology("AMBER", property_map=self._property_map)

        # Create the squashed system.
        if isinstance(self._protocol, _FreeEnergyMixin):
            # Check that the system contains a perturbable molecule.
            if self._system.nPerturbableMolecules() == 0:
                raise ValueError(
                    "'BioSimSpace.Protocol.FreeEnergy' requires a "
                    "perturbable molecule!"
                )

            # Make sure the protocol is valid.
            if self._protocol.getPerturbationType() != "full":
                raise NotImplementedError(
                    "AMBER currently only supports the 'full' perturbation "
                    "type. Please use engine='SOMD' when running multistep "
                    "perturbation types."
                )

            # If this is vacuum simulation with pmemd.cuda then
            # we need to add a simulation box.
            if self._is_vacuum and self._is_pmemd_cuda:
                # Get the existing box information.
                box, _ = system.getBox()

                # We need to add a box.
                if box is None:
                    from ..Box import cubic as _cubic
                    from ..Units.Length import angstrom as _angstrom

                    # Get the bounding box of the system.
                    box_min, box_max = system.getAxisAlignedBoundingBox()

                    # Work out the box size from the difference in the coordinates.
                    box_size = [y - x for x, y in zip(box_min, box_max)]

                    # Work out the size of the box assuming an 8 Angstrom non-bonded cutoff.
                    padding = 8 * _angstrom
                    box_length = max(box_size) + padding
                    # Small box fix. Should be patched in future versions of pmemd.cuda.
                    if box_length < 30 * _angstrom:
                        box_length = 30 * _angstrom

                    # Set the simulation box.
                    system.setBox(*_cubic(box_length))
                    reference_system.setBox(*_cubic(box_length))

            # Apply SOMD1 compatibility to the perturbation.
            if (
                "somd1_compatibility" in kwargs
                and kwargs.get("somd1_compatibility") is True
            ):
                from ._somd import _somd1_compatibility

                system = _somd1_compatibility(system)

            system, self._mapping = _squash(
                system, explicit_dummies=self._explicit_dummies
            )
            self._squashed_system = system

        else:
            # Check for perturbable molecules and convert to the chosen end state.
            system = self._checkPerturbable(system)

        # RST file (coordinates).
        try:
            file = _os.path.splitext(self._rst_file)[0]
            _IO.saveMolecules(file, system, "rst7", property_map=self._property_map)
        except Exception as e:
            msg = "Failed to write system to 'RST7' format."
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Reference file for position restraints.
        try:
            if self._is_real_reference:
                reference_system, _ = _squash(
                    system, explicit_dummies=self._explicit_dummies
                )
                reference_system = self._checkPerturbable(reference_system)
                file = _os.path.splitext(self._ref_file)[0]
                _IO.saveMolecules(
                    file, reference_system, "rst7", property_map=self._property_map
                )
            else:
                _shutil.copy(self._rst_file, self._ref_file)
        except Exception as e:
            msg = "Failed to write reference system to 'RST7' format."
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # PRM file (topology).
        try:
            file = _os.path.splitext(self._top_file)[0]
            _IO.saveMolecules(
                file, system, "prm7", match_water=False, property_map=self._property_map
            )
        except Exception as e:
            msg = "Failed to write system to 'PRM7' format."
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Generate the AMBER configuration file.
        # Skip if the user has passed a custom config.
        if isinstance(self._protocol, _Protocol.Custom):
            self.setConfig(self._protocol.getConfig())
        else:
            self._generate_config()
        self.writeConfig(self._config_file)

        # Generate the dictionary of command-line arguments.
        self._generate_args()

        # Return the list of input files.
        return self._input_files

    def _generate_config(self):
        """Generate AMBER configuration file strings."""
        import shutil as _shutil
        import os as _os
        from .. import Protocol as _Protocol
        from .._Config import Amber as _AmberConfig
        from ._plumed import Plumed as _Plumed

        extra_options = self._extra_options.copy()
        extra_lines = self._extra_lines.copy()

        # Set the random number seed.
        if self._seed is None:
            extra_options["ig"] = -1
        else:
            extra_options["ig"] = self._seed

        # Add configuration variables for a metadynamics simulation.
        if isinstance(self._protocol, (_Protocol.Metadynamics, _Protocol.Steering)):
            extra_options["plumed"] = 1
            extra_options["plumedfile"] = "'plumed.dat'"

            # Create the PLUMED input file and copy auxiliary files to the working directory.
            self._plumed = _Plumed(str(self._work_dir))
            plumed_config, auxiliary_files = self._plumed.createConfig(
                self._system, self._protocol, self._property_map
            )
            self._setPlumedConfig(plumed_config)
            if auxiliary_files is not None:
                for file in auxiliary_files:
                    file_name = _os.path.basename(file)
                    _shutil.copyfile(
                        file, _os.path.join(str(self._work_dir), file_name)
                    )
            self._input_files.append(self._plumed_config_file)

            # Expose the PLUMED specific member functions.
            setattr(self, "getPlumedConfig", self._getPlumedConfig)
            setattr(self, "getPlumedConfigFile", self._getPlumedConfigFile)
            setattr(self, "setPlumedConfig", self._setPlumedConfig)
            setattr(self, "getFreeEnergy", self._getFreeEnergy)
            setattr(self, "getCollectiveVariable", self._getCollectiveVariable)
            setattr(self, "sampleConfigurations", self._sampleConfigurations)
            setattr(self, "getTime", self._getTime)

        # Instantiate the AMBER configuration generator.
        amber_config = _AmberConfig(self._system, self._protocol)

        # Create the configuration.
        self.setConfig(
            amber_config.createConfig(
                is_pmemd=self._is_pmemd,
                is_pmemd_cuda=self._is_pmemd_cuda,
                explicit_dummies=self._explicit_dummies,
                extra_options=extra_options,
                extra_lines=extra_lines,
            )
        )

        # Flag that this isn't a custom protocol.
        if not self._extra_options and not self._extra_lines:
            self._protocol._setCustomised(False)

    def _generate_args(self):
        """Generate the dictionary of command-line arguments."""
        from .. import Protocol as _Protocol
        from ..Protocol._position_restraint_mixin import _PositionRestraintMixin

        # Clear the existing arguments.
        self.clearArgs()

        # Add the default arguments.
        self.setArg("-O", True)  # Overwrite.
        self.setArg("-i", "%s.cfg" % self._name)  # Input file.
        self.setArg("-p", "%s.prm7" % self._name)  # Topology file.
        self.setArg("-c", "%s.rst7" % self._name)  # Coordinate file.
        self.setArg("-o", "%s.out" % self._name)  # Redirect stdout to file.
        self.setArg("-r", "%s.crd" % self._name)  # Restart file.
        self.setArg("-inf", "%s.nrg" % self._name)  # Energy info file.

        # Skip if the user has passed a custom protocol.
        if not isinstance(self._protocol, _Protocol.Custom):
            # Append a reference file if a position restraint is specified.
            if isinstance(self._protocol, _PositionRestraintMixin):
                if self._protocol.getRestraint() is not None:
                    self.setArg("-ref", "%s_ref.rst7" % self._name)

            # Append a trajectory file if this anything other than a minimisation.
            if not isinstance(self._protocol, _Protocol.Minimisation):
                self.setArg("-x", "%s.nc" % self._name)

        # Add the extra arguments.
        for key, value in self._extra_args.items():
            self.setArg(key, value)

    def start(self):
        """
        Start the AMBER process.

        Returns
        -------

        process : :class:`Process.Amber <BioSimSpace.Process.Amber>`
            The process object.
        """
        from sire.legacy import Base as _SireBase
        import timeit as _timeit
        from .. import _Utils

        # The process is currently queued.
        if self.isQueued():
            return

        # Process is already running.
        if self._process is not None:
            if self._process.is_running():
                return

        # Run the process in the working directory.
        with _Utils.cd(self._work_dir):
            # Create the arguments string list.
            args = self.getArgStringList()

            # Write the command-line process to a README.txt file.
            with open("README.txt", "w") as file:
                # Set the command-line string.
                self._command = "%s " % self._exe + self.getArgString()

                # Write the command to file.
                file.write("# AMBER was run with the following command:\n")
                file.write("%s\n" % self._command)

            # Start the timer.
            self._timer = _timeit.default_timer()

            # Start the simulation. Pass a null string for the stdout file
            # since we've explicitly redirected AMBER output to file since
            # pmemd doesn't write to standard output.
            self._process = _SireBase.Process.run(
                self._exe, args, "", "%s.err" % self._name
            )

        return self

    def getSystem(self, block="AUTO"):
        """
        Get the latest molecular system.

        Parameters
        ----------

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The latest molecular system.
        """
        from .._SireWrappers import System as _System
        from sire.legacy import Mol as _SireMol
        import shutil as _shutil
        from ..Align._squash import _unsquash
        from sire.legacy import IO as _SireIO
        import os as _os
        from ..Protocol._free_energy_mixin import _FreeEnergyMixin
        import warnings as _warnings
        import tempfile as _tempfile

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        # Create the name of the restart CRD file.
        restart = _os.path.join(str(self._work_dir), "%s.crd" % self._name)

        # Check that the file exists.
        if _os.path.isfile(restart):
            # Do we need to get coordinates for the lambda=1 state.
            if "is_lambda1" in self._property_map:
                is_lambda1 = True
            else:
                is_lambda1 = False

            # Copy the restart file to a temporary location. Sire streams from
            # binary files with the same path, so we need to ensure that a new
            # stream is created each time.
            with _tempfile.TemporaryDirectory() as tmp_dir:
                tmp_file = f"{tmp_dir}/{self._name}.crd"
                _shutil.copyfile(restart, tmp_file)

                # Create a new molecular system from the restart file.
                new_system = _System(
                    _SireIO.MoleculeParser.read(
                        [tmp_file, self._top_file], self._property_map
                    )
                )

            # Create a copy of the existing system object.
            old_system = self._system.copy()

            if isinstance(self._protocol, _FreeEnergyMixin):
                # Udpate the coordinates and velocities and return a mapping between
                # the molecule indices in the two systems.
                mapping = {
                    _SireMol.MolIdx(x): _SireMol.MolIdx(x)
                    for x in range(0, self._squashed_system.nMolecules())
                }
                (
                    self._squashed_system._sire_object,
                    _,
                ) = _SireIO.updateCoordinatesAndVelocities(
                    self._squashed_system._sire_object,
                    new_system._sire_object,
                    mapping,
                    is_lambda1,
                    self._property_map,
                    self._property_map,
                )

                # Update the unsquashed system based on the updated squashed system.
                old_system = _unsquash(
                    old_system,
                    self._squashed_system,
                    self._mapping,
                    explicit_dummies=self._explicit_dummies,
                )

            else:
                # Update the coordinates and velocities and return a mapping between
                # the molecule indices in the two systems.
                sire_system, mapping = _SireIO.updateCoordinatesAndVelocities(
                    old_system._sire_object,
                    new_system._sire_object,
                    self._mapping,
                    is_lambda1,
                    self._property_map,
                    self._property_map,
                )

                # Update the underlying Sire object.
                old_system._sire_object = sire_system

                # Store the mapping between the MolIdx in both systems so we don't
                # need to recompute it next time.
                self._mapping = mapping

            # Update the box information in the original system.
            if self._has_box:
                if "space" in new_system._sire_object.property_keys():
                    box = new_system._sire_object.property("space")
                    if box.is_periodic():
                        old_system._sire_object.set_property(
                            self._property_map.get("space", "space"), box
                        )

            return old_system

        else:
            return None

    def getCurrentSystem(self):
        """
        Get the latest molecular system.

        Returns
        -------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The latest molecular system.
        """
        return self.getSystem(block=False)

    def getTrajectory(self, backend="AUTO", block="AUTO"):
        """
        Return a trajectory object.

        Parameters
        ----------

        backend : str
            The backend to use for trajectory parsing. To see supported backends,
            run BioSimSpace.Trajectory.backends(). Using "AUTO" will try each in
            sequence.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        trajectory : :class:`Trajectory <BioSimSpace.Trajectory.Trajectory>`
            The latest trajectory object.
        """
        from .. import Trajectory as _Trajectory
        import warnings as _warnings

        if not isinstance(backend, str):
            raise TypeError("'backend' must be of type 'str'")

        if not isinstance(block, (bool, str)):
            raise TypeError("'block' must be of type 'bool' or 'str'")

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        try:
            return _Trajectory.Trajectory(process=self, backend=backend)

        except:
            return None

    def getFrame(self, index):
        """
        Return a specific trajectory frame.

        Parameters
        ----------

        index : int
            The index of the frame.

        Returns
        -------

        frame : :class:`System <BioSimSpace._SireWrappers.System>`
            The System object of the corresponding frame.
        """
        from ..Align._squash import _unsquash
        from sire.legacy import IO as _SireIO
        from .. import Trajectory as _Trajectory
        from .. import Protocol as _Protocol
        from sire.legacy import Mol as _SireMol

        if not type(index) is int:
            raise TypeError("'index' must be of type 'int'")

        max_index = int(
            (self._protocol.getRunTime() / self._protocol.getTimeStep())
            / self._protocol.getRestartInterval()
        )

        if index < 0 or index > max_index:
            raise ValueError(f"'index' must be in range [0, {max_index}].")

        try:
            # Do we need to get coordinates for the lambda=1 state.
            if "is_lambda1" in self._property_map:
                is_lambda1 = True
            else:
                is_lambda1 = False

            # Get the latest trajectory frame.
            new_system = _Trajectory.getFrame(self._traj_file, self._top_file, index)

            # Create a copy of the existing system object.
            old_system = self._system.copy()

            if isinstance(self._protocol, _Protocol._FreeEnergyMixin):
                # Udpate the coordinates and velocities and return a mapping between
                # the molecule indices in the two systems.
                mapping = {
                    _SireMol.MolIdx(x): _SireMol.MolIdx(x)
                    for x in range(0, self._squashed_system.nMolecules())
                }
                (
                    self._squashed_system._sire_object,
                    _,
                ) = _SireIO.updateCoordinatesAndVelocities(
                    self._squashed_system._sire_object,
                    new_system._sire_object,
                    mapping,
                    is_lambda1,
                    self._property_map,
                    self._property_map,
                )

                # Update the unsquashed system based on the updated squashed system.
                old_system = _unsquash(
                    old_system,
                    self._squashed_system,
                    self._mapping,
                    explicit_dummies=self._explicit_dummies,
                )

            else:
                # Update the coordinates and velocities and return a mapping between
                # the molecule indices in the two systems.
                sire_system, mapping = _SireIO.updateCoordinatesAndVelocities(
                    old_system._sire_object,
                    new_system._sire_object,
                    self._mapping,
                    is_lambda1,
                    self._property_map,
                    self._property_map,
                )

                # Update the underlying Sire object.
                old_system._sire_object = sire_system

                # Store the mapping between the MolIdx in both systems so we don't
                # need to recompute it next time.
                self._mapping = mapping

            # Update the box information in the original system.
            if "space" in new_system._sire_object.property_keys():
                box = new_system._sire_object.property("space")
                old_system._sire_object.set_property(
                    self._property_map.get("space", "space"), box
                )

            return old_system

        except:
            return None

    def getRecordKey(self, record, region=0, soft_core=False):
        """
        Parameters
        ----------

        record : str
            The record used in the AMBER standard output, e.g. 'TEMP(K)'.
            Please consult the current AMBER manual for details:
            https://ambermd.org/Manuals.php

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        key : str
            The universal record key that can be used with getRecord.
        """

        # Validate the record string.
        if not isinstance(record, str):
            raise TypeError("'record' must be of type 'str'")

        # Validate the region.
        if not isinstance(region, int):
            raise TypeError("'region' must be of type 'int'")
        else:
            if region < 0 or region > 1:
                raise ValueError("'region' must be in range [0, 1]")

        # Validate the soft-core flag.
        if not isinstance(soft_core, bool):
            raise TypeError("'soft_core' must be of type 'bool'.")

        # Convert to the full index.
        idx = 2 * region + int(soft_core)

        # Strip whitespace from the beginning and end of the record and convert
        # to upper case.
        cleaned_record = record.strip().upper()

        # Make sure the record exists in the key mapping.
        if not cleaned_record in self._stdout_key[idx].values():
            raise ValueError(f"No key found for record '{record}'")

        return list(self._stdout_key[idx].keys())[
            list(self._stdout_key[idx].values()).index(cleaned_record)
        ]

    def getRecord(
        self, key, time_series=False, unit=None, region=0, soft_core=False, block="AUTO"
    ):
        """
        Get a record from the stdout dictionary.

        Parameters
        ----------

        key : str
            A universal record key based on the key used in the AMBER standard
            output. Use 'getRecordKey(record)` to generate the key. The records
            are those used in the AMBER standard output, e.g. 'TEMP(K)'. Please
            consult the current AMBER manual for details: https://ambermd.org/Manuals.php

        time_series : bool
            Whether to return a list of time series records.

        unit : :class:`Unit <BioSimSpace.Units>`
            The unit to convert the record to.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        record : :class:`Type <BioSimSpace.Types>`
            The matching record.
        """
        import warnings as _warnings

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        return self._get_stdout_record(
            key.strip().upper(),
            time_series=time_series,
            unit=unit,
            region=region,
            soft_core=soft_core,
        )

    def getCurrentRecord(
        self, key, time_series=False, unit=None, region=0, soft_core=False
    ):
        """
        Get a current record from the stdout dictionary.

        Parameters
        ----------

        key : str
            A universal record key based on the key used in the AMBER standard
            output. Use 'getRecordKey(record)` to generate the key. The records
            are those used in the AMBER standard output, e.g. 'TEMP(K)'. Please
            consult the current AMBER manual for details: https://ambermd.org/Manuals.php

        time_series : bool
            Whether to return a list of time series records.

        unit : :class:`Unit <BioSimSpace.Units>`
            The unit to convert the record to.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        record : :class:`Type <BioSimSpace.Types>`
            The matching record.
        """
        import warnings as _warnings

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        return self._get_stdout_record(
            key.strip().upper(),
            time_series=time_series,
            unit=unit,
            region=region,
            soft_core=soft_core,
        )

    def getRecords(self, region=0, soft_core=False, block="AUTO"):
        """
        Return the dictionary of stdout time-series records.

        Parameters
        ----------

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        records : :class:`MultiDict <BioSimSpace.Process._process._MultiDict>`
           The dictionary of time-series records.
        """
        import warnings as _warnings

        # Validate the region.
        if not isinstance(region, int):
            raise TypeError("'region' must be of type 'int'")
        else:
            if region < 0 or region > 1:
                raise ValueError("'region' must be in range [0, 1]")

        # Validate the soft-core flag.
        if not isinstance(soft_core, bool):
            raise TypeError("'soft_core' must be of type 'bool'.")

        # Convert to the full index, region + soft_core.
        idx = 2 * region + int(soft_core)

        # Wait for the process to finish.
        if block is True:
            self.wait()
        elif block == "AUTO" and self._is_blocked:
            self.wait()

        # Warn the user if the process has exited with an error.
        if self.isError():
            _warnings.warn("The process exited with an error!")

        self.stdout(0)

        return self._stdout_dict[idx].copy()

    def getCurrentRecords(self, region=0, soft_core=False):
        """
        Return the current dictionary of stdout time-series records.

        Parameters
        ----------

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        records : :class:`MultiDict <BioSimSpace.Process._process._MultiDict>`
           The dictionary of time-series records.
        """
        return self.getRecords(region=region, soft_core=soft_core, block=False)

    def getTime(self, time_series=False, region=0, soft_core=False, block="AUTO"):
        """
        Get the simulation time.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The current simulation time in nanoseconds.
        """
        from .. import Protocol as _Protocol
        from .. import Units as _Units

        # No time records for minimisation protocols.
        if isinstance(self._protocol, _Protocol.Minimisation):
            return None

        # Get the list of time steps.
        time_steps = self.getRecord(
            "TIME(PS)",
            time_series=time_series,
            unit=None,
            region=region,
            soft_core=soft_core,
            block=block,
        )

        # Convert from picoseconds to nanoseconds.
        if time_steps is not None:
            if time_series:
                return [
                    (x * _Units.Time.picosecond)._to_default_unit() for x in time_steps
                ]
            else:
                return (time_steps * _Units.Time.picosecond)._to_default_unit()

    def getCurrentTime(self, time_series=False, region=0, soft_core=False):
        """
        Get the current simulation time.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The current simulation time in nanoseconds.
        """
        return self.getTime(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getStep(self, time_series=False, region=0, soft_core=False, block="AUTO"):
        """
        Get the number of integration steps.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        step : int
            The current number of integration steps.
        """
        return self.getRecord(
            "NSTEP",
            time_series=time_series,
            unit=None,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentStep(self, time_series=False, region=0, soft_core=False):
        """
        Get the current number of integration steps.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        step : int
            The current number of integration steps.
        """
        return self.getStep(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getBondEnergy(self, time_series=False, region=0, soft_core=False, block="AUTO"):
        """
        Get the bond energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The bond energy.
        """
        from .. import Units as _Units

        return self.getRecord(
            "BOND",
            time_series=time_series,
            unit=_Units.Energy.kcal_per_mol,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentBondEnergy(self, time_series=False, region=0, soft_core=False):
        """
        Get the current bond energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The bond energy.
        """
        return self.getBondEnergy(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getAngleEnergy(
        self, time_series=False, region=0, soft_core=False, block="AUTO"
    ):
        """
        Get the angle energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The angle energy.
        """
        from .. import Units as _Units

        return self.getRecord(
            "ANGLE",
            time_series=time_series,
            unit=_Units.Energy.kcal_per_mol,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentAngleEnergy(self, time_series=False, region=0, soft_core=False):
        """
        Get the current angle energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The angle energy.
        """
        return self.getAngleEnergy(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getDihedralEnergy(
        self, time_series=False, region=0, soft_core=False, block="AUTO"
    ):
        """
        Get the total dihedral energy (proper + improper).

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The total dihedral energy.
        """
        from .. import Units as _Units

        return self.getRecord(
            "DIHED",
            time_series=time_series,
            unit=_Units.Energy.kcal_per_mol,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentDihedralEnergy(self, time_series=False, region=0, soft_core=False):
        """
        Get the current total dihedral energy (proper + improper).

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The total dihedral energy.
        """
        return self.getDihedralEnergy(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getElectrostaticEnergy(
        self, time_series=False, region=0, soft_core=False, block="AUTO"
    ):
        """
        Get the electrostatic energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The electrostatic energy.
        """
        from .. import Units as _Units

        return self.getRecord(
            "EEL",
            time_series=time_series,
            unit=_Units.Energy.kcal_per_mol,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentElectrostaticEnergy(
        self, time_series=False, region=0, soft_core=False
    ):
        """
        Get the current dihedral energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The electrostatic energy.
        """
        return self.getElectrostaticEnergy(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getElectrostaticEnergy14(
        self, time_series=False, region=0, soft_core=False, block="AUTO"
    ):
        """
        Get the electrostatic energy between atoms 1 and 4.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The electrostatic energy between atoms 1 and 4.
        """
        from .. import Units as _Units

        return self.getRecord(
            "14EEL",
            time_series=time_series,
            unit=_Units.Energy.kcal_per_mol,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentElectrostaticEnergy14(
        self, time_series=False, region=0, soft_core=False
    ):
        """
        Get the current electrostatic energy between atoms 1 and 4.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The electrostatic energy between atoms 1 and 4.
        """
        return self.getElectrostaticEnergy14(
            time_series=time_series, region=region, soft_core=False, block=False
        )

    def getVanDerWaalsEnergy(
        self, time_series=False, region=0, soft_core=False, block="AUTO"
    ):
        """
        Get the Van der Vaals energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The Van der Vaals energy.
        """
        from .. import Units as _Units

        return self.getRecord(
            "VDW",
            time_series=time_series,
            unit=_Units.Energy.kcal_per_mol,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentVanDerWaalsEnergy(self, time_series=False, region=0, soft_core=False):
        """
        Get the current Van der Vaals energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The Van der Vaals energy.
        """
        return self.getVanDerWaalsEnergy(
            time_series=time_series, block=False, region=region, soft_core=soft_core
        )

    def getVanDerWaalsEnergy14(
        self, time_series=False, region=0, soft_core=False, block="AUTO"
    ):
        """
        Get the Van der Vaals energy between atoms 1 and 4.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The Van der Vaals energy between atoms 1 and 4.
        """
        from .. import Units as _Units

        return self.getRecord(
            "14VDW",
            time_series=time_series,
            unit=_Units.Energy.kcal_per_mol,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentVanDerWaalsEnergy14(
        self, time_series=False, region=0, soft_core=False
    ):
        """
        Get the current Van der Vaals energy between atoms 1 and 4.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The Van der Vaals energy between atoms 1 and 4.
        """
        return self.getVanDerWaalsEnergy(
            time_series=time_series, block=False, region=region, soft_core=soft_core
        )

    def getHydrogenBondEnergy(
        self, time_series=False, region=0, soft_core=False, block="AUTO"
    ):
        """
        Get the hydrogen bond energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The hydrogen bond energy.
        """
        from .. import Units as _Units

        return self.getRecord(
            "EHBOND",
            time_series=time_series,
            unit=_Units.Energy.kcal_per_mol,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentHydrogenBondEnergy(
        self, time_series=False, region=0, soft_core=False
    ):
        """
        Get the current hydrogen bond energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The hydrogen bond energy.
        """
        return self.getHydrogenBondEnergy(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getRestraintEnergy(
        self, time_series=False, region=0, soft_core=False, block="AUTO"
    ):
        """
        Get the restraint energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The restraint energy.
        """
        from .. import Units as _Units

        return self.getRecord(
            "RESTRAINT",
            time_series=time_series,
            unit=_Units.Energy.kcal_per_mol,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentRestraintEnergy(self, time_series=False, region=0, soft_core=False):
        """
        Get the current restraint energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The restraint energy.
        """
        return self.getRestraintEnergy(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getPotentialEnergy(
        self, time_series=False, region=0, soft_core=False, block="AUTO"
    ):
        """
        Get the potential energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The potential energy.
        """
        from .. import Units as _Units

        return self.getRecord(
            "EPTOT",
            time_series=time_series,
            unit=_Units.Energy.kcal_per_mol,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentPotentialEnergy(self, time_series=False, region=0, soft_core=False):
        """
        Get the current potential energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The potential energy.
        """
        return self.getPotentialEnergy(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getKineticEnergy(
        self, time_series=False, region=0, soft_core=False, block="AUTO"
    ):
        """
        Get the kinetic energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The kinetic energy.
        """
        from .. import Units as _Units

        return self.getRecord(
            "EKTOT",
            time_series=time_series,
            unit=_Units.Energy.kcal_per_mol,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentKineticEnergy(self, time_series=False, region=0, soft_core=False):
        """
        Get the current kinetic energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The kinetic energy.
        """
        return self.getKineticEnergy(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getNonBondedEnergy14(
        self, time_series=False, region=0, soft_core=False, block="AUTO"
    ):
        """
        Get the non-bonded energy between atoms 1 and 4.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The non-bonded energy between atoms 1 and 4.
        """
        from .. import Units as _Units

        return self.getRecord(
            "14NB",
            time_series=time_series,
            unit=_Units.Energy.kcal_per_mol,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentNonBondedEnergy14(self, time_series=False, region=0, soft_core=False):
        """
        Get the current non-bonded energy between atoms 1 and 4.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The non-bonded energy between atoms 1 and 4.
        """
        return self.getNonBondedEnergy14(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getTotalEnergy(
        self, time_series=False, region=0, soft_core=False, block="AUTO"
    ):
        """
        Get the total energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The total energy.
        """
        from .. import Protocol as _Protocol
        from .. import Units as _Units

        if not isinstance(region, int):
            raise TypeError("'region' must be of type 'int'")
        else:
            if region < 0 or region > 1:
                raise ValueError("'region' must be in range [0, 1]")

        # Validate the soft-core flag.
        if not isinstance(soft_core, bool):
            raise TypeError("'soft_core' must be of type 'bool'.")

        # Convert to the full index, region + soft_core.
        idx = 2 * region + int(soft_core)

        if isinstance(self._protocol, _Protocol.Minimisation) and not soft_core:
            return self.getRecord(
                "ENERGY",
                time_series=time_series,
                unit=_Units.Energy.kcal_per_mol,
                region=region,
                soft_core=soft_core,
                block=block,
            )
        else:
            return self.getRecord(
                "ETOT",
                time_series=time_series,
                unit=_Units.Energy.kcal_per_mol,
                region=region,
                soft_core=soft_core,
                block=block,
            )

    def getCurrentTotalEnergy(self, time_series=False, region=0, soft_core=False):
        """
        Get the current total energy.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The total energy.
        """
        return self.getTotalEnergy(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getCentreOfMassKineticEnergy(
        self, time_series=False, region=0, soft_core=False, block="AUTO"
    ):
        """
        Get the kinetic energy of the centre of mass in translation.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The centre of mass kinetic energy.
        """
        from .. import Units as _Units

        return self.getRecord(
            "EKCMT",
            time_series=time_series,
            unit=_Units.Energy.kcal_per_mol,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentCentreOfMassKineticEnergy(
        self, time_series=False, region=0, soft_core=False
    ):
        """
        Get the current kinetic energy of the centre of mass in translation.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
           The centre of mass kinetic energy.
        """
        return self.getCentreOfMassKineticEnergy(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getVirial(self, time_series=False, region=0, soft_core=False, block="AUTO"):
        """
        Get the virial.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        virial : float
           The virial.
        """
        return self.getRecord(
            "VIRIAL",
            time_series=time_series,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentVirial(self, time_series=False, region=0, soft_core=False):
        """
        Get the current virial.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        virial : float
           The virial.
        """
        return self.getVirial(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getTemperature(
        self, time_series=False, region=0, soft_core=False, block="AUTO"
    ):
        """
        Get the temperature.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
           The temperature.
        """
        from .. import Units as _Units

        return self.getRecord(
            "TEMP(K)",
            time_series=time_series,
            unit=_Units.Temperature.kelvin,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentTemperature(self, time_series=False, region=0, soft_core=False):
        """
        Get the current temperature.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
           The temperature.
        """
        return self.getTemperature(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getPressure(self, time_series=False, region=0, soft_core=False, block="AUTO"):
        """
        Get the pressure.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
           The pressure.
        """
        from .. import Units as _Units

        return self.getRecord(
            "PRESS",
            time_series=time_series,
            unit=_Units.Pressure.bar,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentPressure(self, time_series=False, region=0, soft_core=False):
        """
        Get the current pressure.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
           The pressure.
        """
        return self.getPressure(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getVolume(self, time_series=False, region=0, soft_core=False, block="AUTO"):
        """
        Get the volume.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        volume : :class:`Volume <BioSimSpace.Types.Volume>`
           The volume.
        """
        from .. import Units as _Units

        return self.getRecord(
            "VOLUME",
            time_series=time_series,
            unit=_Units.Volume.angstrom3,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentVolume(self, time_series=False, region=0, soft_core=False):
        """
        Get the current volume.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        volume : :class:`Volume <BioSimSpace.Types.Volume>`
           The volume.
        """
        return self.getVolume(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getDensity(self, time_series=False, region=0, soft_core=False, block="AUTO"):
        """
        Get the density.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        density : float
           The density.
        """
        return self.getRecord(
            "DENSITY",
            time_series=time_series,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentDensity(self, time_series=False, region=0, soft_core=False):
        """
        Get the current density.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        density : float
           The density.
        """
        return self.getDensity(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def getDVDL(self, time_series=False, region=0, soft_core=False, block="AUTO"):
        """
        Get the gradient of the total energy with respect to lambda.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        block : bool
            Whether to block until the process has finished running.

        Returns
        -------

        dv_dl : float
            The gradient of the total energy with respect to lambda.
        """
        return self.getRecord(
            "DVDL",
            time_series=time_series,
            region=region,
            soft_core=soft_core,
            block=block,
        )

    def getCurrentDVDL(self, time_series=False, region=0, soft_core=False):
        """
        Get the current gradient of the total energy with respect to lambda.

        Parameters
        ----------

        time_series : bool
            Whether to return a list of time series records.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        dv_dl : float
            The current gradient of the total energy with respect to lambda.
        """
        return self.getDVDL(
            time_series=time_series, region=region, soft_core=soft_core, block=False
        )

    def stdout(self, n=10):
        """
        Print the last n lines of the stdout buffer.

        Parameters
        ----------

        n : int
            The number of lines to print.
        """
        from .. import Protocol as _Protocol
        from ..Protocol._free_energy_mixin import _FreeEnergyMixin
        import re as _re

        # Ensure that the number of lines is positive.
        if n < 0:
            raise ValueError("The number of lines must be positive!")

        # Flag that this isn't a header line.
        self._is_header = False

        # Append any new lines to the stdout list.
        for line in _pygtail.Pygtail(self._stdout_file):
            self._stdout.append(line.rstrip())
            line = line.strip()

            # Swap dictionary based on the protocol and the degre of freedom to
            # which the next block of records correspond.
            if isinstance(self._protocol, _FreeEnergyMixin):
                if "TI region  1" in line:
                    self._current_region = 0
                elif "TI region  2" in line:
                    self._current_region = 2
                elif "Softcore part" in line and self._current_region == 0:
                    self._current_region = 1
                elif "Softcore part" in line and self._current_region == 2:
                    self._current_region = 3
                elif "Detailed TI info" in line:
                    # This flags that we should skip records until the start of
                    # the next set for the first TI region.
                    self._current_region = 4
            # Default stdout dictionary.
            else:
                self._current_region = 0

            # Continue if we are ignoring this record block.
            if self._current_region == 4:
                continue

            stdout_dict = self._stdout_dict[self._current_region]
            stdout_key = self._stdout_key[self._current_region]

            # Skip empty lines and summary reports.
            if len(line) > 0 and line[0] != "|" and line[0] != "-":
                # Skip EAMBER records.
                if "EAMBER (non-restraint)" in line:
                    continue
                # Flag that we've started recording results.
                elif not self._has_results and line.startswith("NSTEP"):
                    self._has_results = True
                    self._finished_results = False
                # Flag that we've finished recording results.
                elif "A V E R A G E S" in line:
                    self._finished_results = True

                # Parse the results.
                if self._has_results and not self._finished_results:
                    # The first line of output has different formatting for minimisation protocols.
                    if isinstance(self._protocol, _Protocol.Minimisation):
                        # No equals sign in the line.
                        if "NSTEP" in line and "=" not in line:
                            # Split the line using whitespace.
                            data = line.upper().split()

                            # If we find a header, jump to the top of the loop.
                            if len(data) > 0:
                                if data[0] == "NSTEP":
                                    self._is_header = True
                                    continue

                        # Process the header record.
                        if self._is_header:
                            # Split the line using whitespace.
                            data = line.upper().split()

                            # The file hasn't been updated.
                            if (
                                "NSTEP" in stdout_dict
                                and data[0] == stdout_dict["NSTEP"][-1]
                            ):
                                self._finished_results = True
                                continue

                            # Add the timestep and energy records to the dictionary.
                            stdout_dict["NSTEP"] = data[0]
                            stdout_dict["ENERGY"] = data[1]

                            # Add the keys to the mapping
                            stdout_key["NSTEP"] = "NSTEP"
                            stdout_key["ENERGY"] = "ENERGY"

                            # Turn off the header flag now that the data has been recorded.
                            self._is_header = False

                    # All other records are formatted as RECORD = VALUE.

                    # Use a regex search to split the line into record names and values.
                    records = _re.findall(
                        r"([SC_]*[EEL_]*[RES_]*[VDW_]*\d*\-*\d*\s*[A-Z/]+\(*[A-Z]*\)*)\s*=\s*(\-*\d+\.?\d*|\**)",
                        line.upper(),
                    )

                    # Append each record to the dictionary.
                    for key, value in records:
                        # Strip whitespace from beginning and end.
                        key = key.strip()

                        # Format key so it can be re-used for records corresponding to
                        # different regions, which use different abbreviations.
                        universal_key = (
                            key.replace("SC_", "")
                            .replace(" ", "")
                            .replace("-", "")
                            .replace("EELEC", "EEL")
                            .replace("VDWAALS", "VDW")
                        )

                        # Handle missing values, which will appear as asterisks, e.g.
                        # PRESS=********
                        try:
                            tmp = float(value)
                        except:
                            value = None

                        # Store the record using the original key.
                        stdout_dict[key] = value

                        # Map the universal key to the original.
                        stdout_key[universal_key] = key

        # Get the current number of lines.
        num_lines = len(self._stdout)

        # Set the line from which to start printing.
        if num_lines < n:
            start = 0
        else:
            start = num_lines - n

        # Print the lines.
        for x in range(start, num_lines):
            print(self._stdout[x])

    def kill(self):
        """Kill the running process."""

        # Kill the process.
        if not self._process is None and self._process.is_running():
            self._process.kill()

    def _get_stdout_record(
        self, key, time_series=False, unit=None, region=0, soft_core=False
    ):
        """
        Helper function to get a stdout record from the dictionary.

        Parameters
        ----------

        key : str
            The universal record key.

        time_series : bool
            Whether to return a time series of records.

        unit : :class:`Type <BioSimSpace.Types._type.Type>`
            The unit to convert the record to.

        region : int
            The region to which the record corresponds. There will only be more
            than one region for FreeEnergy protocols, where 1 indicates the second
            TI region.

        soft_core : bool
            Whether to get the record for the soft-core part of the system for the
            chosen region.

        Returns
        -------

        record :
            The matching stdout record.
        """
        from ..Types._type import Type as _Type
        import warnings as _warnings

        # Update the standard output dictionary.
        self.stdout(0)

        # No data!
        if len(self._stdout_dict) == 0:
            return None

        if not isinstance(time_series, bool):
            _warnings.warn("Non-boolean time-series flag. Defaulting to False!")
            time_series = False

        # Validate the unit.
        if unit is not None:
            if not isinstance(unit, _Type):
                raise TypeError("'unit' must be of type 'BioSimSpace.Types'")

        # Validate the region.
        if not isinstance(region, int):
            raise TypeError("'region' must be of type 'int'")
        else:
            if region < 0 or region > 1:
                raise ValueError("'region' must be in range [0, 1]")

        # Validate the soft-core flag.
        if not isinstance(soft_core, bool):
            raise TypeError("'soft_core' must be of type 'bool'.")

        # Convert to the full index, region + soft_core.
        idx = 2 * region + int(soft_core)

        # Extract the dictionary of stdout records for the specified region and soft-core flag.
        stdout_dict = self._stdout_dict[idx]

        # Map the universal key to the original key used for this region.
        try:
            key = self._stdout_key[idx][key]
        except:
            return None

        # Return the list of dictionary values.
        if time_series:
            try:
                if key == "NSTEP":
                    return [int(x) for x in stdout_dict[key]]
                else:
                    if unit is None:
                        return [float(x) if x else None for x in stdout_dict[key]]
                    else:
                        return [
                            (float(x) * unit)._to_default_unit() if x else None
                            for x in stdout_dict[key]
                        ]

            except KeyError:
                return None

        # Return the most recent dictionary value.
        else:
            try:
                if key == "NSTEP":
                    return int(stdout_dict[key][-1])
                else:
                    if unit is None:
                        try:
                            return float(stdout_dict[key][-1])
                        except:
                            return None
                    else:
                        try:
                            return (
                                float(stdout_dict[key][-1]) * unit
                            )._to_default_unit()
                        except:
                            return None

            except KeyError:
                return None


def _findExe(is_gpu=False, is_free_energy=False, is_vacuum=False):
    """
    Helper function to search for an AMBER executable.

    Parameters
    ----------

    is_gpu : bool
        Whether to search for a GPU-enabled executable.

    is_free_energy : bool
        Whether this is a free energy simulation.

    is_vacuum : bool
        Whether this is a vacuum simulation.

    Returns
    -------

    exe : str
        The path to the executable.
    """
    from .. import _amber_home
    from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
    import os as _os

    if not isinstance(is_gpu, bool):
        raise TypeError("'is_gpu' must be of type 'bool'.")

    if not isinstance(is_free_energy, bool):
        raise TypeError("'is_free_energy' must be of type 'bool'.")

    if not isinstance(is_vacuum, bool):
        raise TypeError("'is_vacuum' must be of type 'bool'.")

    if is_gpu:
        targets = ["pmemd.cuda"]
    else:
        if is_free_energy:
            targets = ["pmemd"]
        else:
            targets = ["pmemd", "sander"]

    # Search for the executable.

    import os as _os
    import pathlib as _pathlib

    from glob import glob as _glob

    # Get the current path.
    path = _os.environ["PATH"].split(_os.pathsep)

    # If AMBERHOME is set, then prepend to the path.
    if _amber_home is not None:
        path = [_amber_home + "/bin"] + path

    # Helper function to check whether a file is executable.
    def is_exe(fpath):
        import os as _os

        return _os.path.isfile(fpath) and _os.access(fpath, _os.X_OK)

    # Loop over each directory in the path and search for the executable.
    for p in path:
        # Loop over each target.
        for t in targets:
            # Glob for the executable.
            results = _glob(f"{t}*", root_dir=p)
            # If we find a match, check that it's executable and return the path.
            # Note that this returns the first match, not the best match. If a
            # user requires a specific version of the executable, they should
            # order their path accordingly, or use the exe keyword argument.
            if results:
                for exe in results:
                    # Exclude "locally enhanced sampling" executables.
                    if "LES" not in exe:
                        exe = _pathlib.Path(p) / exe
                        if is_exe(exe):
                            return str(exe)

    msg = (
        "'BioSimSpace.Process.Amber' is not supported. "
        "Unable to find AMBER executable in AMBERHOME or PATH. "
        "Please install AMBER (http://ambermd.org)."
    )

    if is_free_energy:
        msg += " Free energy simulations require 'pmemd' or 'pmemd.cuda'."

    # If we don't find the executable, raise an error.
    raise _MissingSoftwareError(msg)
