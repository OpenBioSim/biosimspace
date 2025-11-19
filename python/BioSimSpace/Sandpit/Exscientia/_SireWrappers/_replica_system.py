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
A thin wrapper around Sire.System.System and sire.mol._trjajectory.TrajectoryIterator.
This is an internal package and should not be directly exposed to the user.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["ReplicaSystem"]


class ReplicaSystem:
    """
    A container class for storing molecular systems for replica dynamics simulations.
    Here a system has a single topology, but multiple coordinate sets (replicas).
    """

    def __init__(self, system, trajectory=None, num_replicas=None):
        """
        Constructor.

        Parameters
        ----------

        system : sire.system.System, sire.legacy.System.System, :class:`System <BioSimSpace._SireWrappers.System>`
            A Sire or BioSimSpace System object.

        trajectory : str, optional
            The path to a trajectory file containing multiple coordinate sets.

        num_replicas : int, optional
            The number of replicas (coordinate sets) to create. If provided, the
            coordinate set in `system` will be duplicated this many times. This is
            only used if the system does not already contain multiple frames or a
            trajectory is not provided.
        """

        from sire.mol._trajectory import TrajectoryIterator as _TrajectoryIterator
        from sire.system import System as _NewSireSystem
        from sire.legacy.System import System as _SireSystem
        from ._system import System as _System

        # Check that the system is valid.
        if isinstance(system, _NewSireSystem):
            self._new_sire_object = system
            self._sire_object = system._system
        elif isinstance(system, _SireSystem):
            self._new_sire_object = _NewSireSystem(system)
            self._sire_object = system
        elif isinstance(system, _System):
            self._new_sire_object = _NewSireSystem(system._sire_object)
            self._sire_object = system._sire_object
        else:
            raise TypeError(
                "'system' must be of type 'sire.system.System', "
                "'sire.legacy.System.System' or 'BioSimSpace._SireWrappers.System'"
            )

        # Check if the system is perturbable.
        try:
            perturbable = self._new_sire_object["perturbable"]
            self._is_perturbable = True

            from sire.morph import link_to_reference as _link_to_reference

            self._new_sire_object = _link_to_reference(self._new_sire_object)
        except:
            self._is_perturbable = False

        # If a trajectory is provided, make sure the file exists.
        if trajectory is not None:
            if self._new_sire_object.num_frames() > 1:
                raise ValueError(
                    "Cannot provide a trajectory when 'system' already contains multiple frames."
                )

            is_traj_iterator = False
            if isinstance(trajectory, _TrajectoryIterator):
                self._trajectory = trajectory
                is_traj_iterator = True

            if not is_traj_iterator:
                import os as _os
                from tempfile import NamedTemporaryFile as _NamedTemporaryFile
                from .. import _isVerbose

                if not _os.path.isfile(trajectory):
                    raise IOError(f"'trajectory' does not exist: {trajectory}")

                # Get the extension of the file so that we can guess the format.
                _, ext = _os.path.splitext(trajectory)

                # Use a GROMACS topology for XTC files.
                if self._is_perturbable or ext == ".xtc":
                    tmp_top = _NamedTemporaryFile(delete=False, suffix=".top")
                # For now, use AMBER PRM7 for all other formats.
                else:
                    tmp_top = _NamedTemporaryFile(delete=False, suffix=".prm7")

                # Try to write out a temporary topology file.
                try:
                    from sire import save as _save

                    _save(self._new_sire_object, tmp_top.name, silent=True)
                except Exception as e:
                    tmp_top.close()
                    _os.unlink(tmp_top.name)

                    msg = "Failed to write temporary topology file for trajectory loading."

                    if _isVerbose():
                        raise IOError(msg) from e
                    else:
                        raise IOError(msg)

                # Now try to load the trajectory.
                try:
                    from sire import load as _load

                    sire_system = _load(tmp_top.name, trajectory, silent=True)
                except Exception as e:
                    msg = "Failed to load trajectory file: %s" % trajectory

                    if _isVerbose():
                        raise IOError(msg) from e
                    else:
                        raise IOError(msg)

                # Update the internal trajectory.
                self._trajectory = sire_system.trajectory()

                # Clean up the temporary topology file.
                tmp_top.close()
                _os.unlink(tmp_top.name)

        # Store the trajectory if there are multiple frames in the system.
        elif self._new_sire_object.num_frames() > 1:
            self._trajectory = self._new_sire_object.trajectory()
        # Otherwise, duplicate the system if num_replicas is provided.
        elif num_replicas is not None:
            if not isinstance(num_replicas, int) or num_replicas < 1:
                raise ValueError("'num_replicas' must be a positive integer.")

            import os as _os
            from tempfile import TemporaryDirectory as _TemporaryDirectory
            from sire import load as _load
            from sire import save as _save
            from .. import _isVerbose

            with _TemporaryDirectory() as tmp_dir:
                filenames = []
                # Write out the current coordinates the specified number of times.
                try:
                    # Write the first file only.
                    tmp_file = _os.path.join(tmp_dir, f"replica_0000.gro")
                    _save(self._new_sire_object, tmp_file, silent=True)
                    filenames.append(tmp_file)

                    # Copy the first file for the remaining replicas.
                    for i in range(1, num_replicas):
                        tmp_file = _os.path.join(tmp_dir, f"replica_{i:04d}.gro")
                        _os.link(_os.path.join(tmp_dir, f"replica_0000.gro"), tmp_file)
                        filenames.append(tmp_file)
                except Exception as e:
                    msg = "Failed to write temporary files for replica duplication."

                    if _isVerbose():
                        raise IOError(msg) from e
                    else:
                        raise IOError(msg)

                if self._is_perturbable:
                    top_ext = "top"
                else:
                    top_ext = "prm7"

                # Write out the topology once.
                try:
                    tmp_top = _os.path.join(tmp_dir, f"topology.{top_ext}")
                    filenames.append(tmp_top)
                    _save(self._new_sire_object, tmp_top, silent=True)
                except Exception as e:
                    msg = "Failed to write temporary topology file for replica duplication."

                    if _isVerbose():
                        raise IOError(msg) from e
                    else:
                        raise IOError(msg)

                # Now load them all back in as a single system with multiple frames.
                try:
                    sire_system = _load(filenames, silent=True)
                except Exception as e:
                    msg = "Failed to load temporary files for replica duplication."

                    if _isVerbose():
                        raise IOError(msg) from e
                    else:
                        raise IOError(msg)

            # Update the internal trajectory.
            self._trajectory = sire_system.trajectory()

        else:
            raise ValueError(
                "Either 'trajectory' or 'num_replicas' must be provided "
                "if 'system' does not already contain multiple frames."
            )

        # Whether the system is linked to the perturbed state.
        self._is_perturbed = False

    def __str__(self):
        """Return a human readable string representation of the object."""
        return "<BioSimSpace.ReplicaSystem: nReplicas=%d>" % self.nReplicas()

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return "<BioSimSpace.ReplicaSystem: nReplicas=%d>" % self.nReplicas()

    def __getitem__(self, index):
        """
        Get a specific replica (coordinate set) from the system.

        Parameters
        ----------

        index : int
            The index of the replica to retrieve (0-based).

        Returns
        -------

        :class:`System <BioSimSpace._SireWrappers.System>`
            The requested replica as a BioSimSpace System object.
        """

        return self.getReplica(index, is_lambda1=self._is_perturbed)

    def nReplicas(self):
        """
        Get the number of replicas (coordinate sets) in the system.

        Returns
        -------

        int
            The number of replicas.
        """

        return self._trajectory.num_frames()

    def save(self, filename, traj_format="dcd", save_velocities=False):
        """
        Save the replica system to a stream and trajectory file.

        Parameters
        ----------

        filename : str
            The base filename.

        traj_format : str, optional
            The format to save the trajectory in. Default is 'dcd'.
            Options are "dcd" or "xtc".

        save_velocities : bool, optional
            Whether to save velocities along with coordinates. Default is False.

        Returns
        -------

        stream : str
            The path to the saved stream file.

        trajectory : str
            The path to the saved trajectory file.
        """

        from sire import save as _save
        from sire.stream import save as _save_stream
        from .. import _isVerbose

        if not isinstance(filename, str):
            raise TypeError("'filename' must be of type 'str'.")

        # Validate the trajectory format.
        valid_formats = ["dcd", "xtc"]

        if not isinstance(traj_format, str):
            raise TypeError("'traj_format' must be of type 'str'.")

        traj_format = traj_format.lower().replace(" ", "")
        if traj_format not in valid_formats:
            raise ValueError(
                f"'traj_format' must be one of: {', '.join(valid_formats)}"
            )

        if not isinstance(save_velocities, bool):
            raise TypeError("'save_velocities' must be of type 'bool'.")

        # Save the trajectory first.
        try:
            traj_filename = f"{filename}.{traj_format}"
            _save(
                self._trajectory,
                traj_filename,
                save_velocities=save_velocities,
                silent=True,
            )
        except Exception as e:
            msg = "Failed to save trajectory file."

            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg)

        # Now clone the new Sire system so that we can remove the trajectory.
        system = self._new_sire_object.clone()
        system.delete_all_frames()

        # Stream the system to file.
        try:
            stream_filename = f"{filename}.bss"
            _save_stream(system, stream_filename)
        except Exception as e:
            msg = "Failed to save stream file."

            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg)

        return stream_filename, traj_filename

    @staticmethod
    def load(stream, trajectory):
        """
        Load a replica system from a stream and trajectory file.

        Parameters
        ----------

        stream : str
            The path to the stream file.

        trajectory : str
            The path to the trajectory file.

        Returns
        -------

        :class:`ReplicaSystem <BioSimSpace._SireWrappers.ReplicaSystem>`
            The loaded replica system.
        """

        from sire.stream import load as _load
        from .. import _isVerbose

        if not isinstance(stream, str):
            raise TypeError("'stream' must be of type 'str'.")
        if not isinstance(trajectory, str):
            raise TypeError("'trajectory' must be of type 'str'.")

        # Try to load the system.
        try:
            sire_system = _load(stream)
        except Exception as e:
            msg = "Failed to load replica system from stream and trajectory."

            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg)

        return ReplicaSystem(sire_system, trajectory=trajectory)

    def getReplica(self, index, is_lambda1=False, property_map={}):
        """
        Get a specific replica (coordinate set) from the system.

        Parameters
        ----------

        index : int
            The index of the replica to retrieve (0-based).

        is_lambda1 : bool, optional
            Whether to return the system at lambda = 1. Default is False.

        property_map : dict, optional
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        :class:`System <BioSimSpace._SireWrappers.System>`
            The requested replica as a BioSimSpace System object.
        """

        from ._system import System as _System
        from sire.morph import link_to_perturbed as _link_to_perturbed

        if not isinstance(index, int):
            raise TypeError("'index' must be an integer.")

        if index < 0:
            index += self.nReplicas()

        if index < 0 or index >= self.nReplicas():
            raise IndexError(
                f"'index' {index} is out of range [0, {self.nReplicas()-1}]."
            )

        # Clone the new Sire system to avoid modifying the original.
        replica_system = self._new_sire_object.clone()

        # Link to the perturbed state, if requested.
        if is_lambda1:
            try:
                replica_system = _link_to_perturbed(replica_system)
            except:
                pass

        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'.")

        # Get the specific frame.
        frame = self._trajectory[index].current()

        from sire.io import get_coords_array as _get_coords_array
        from sire.legacy.IO import setCoordinates as _setCoordinates

        # Copy the frame coordinates into the system.
        frame = _setCoordinates(
            replica_system._system,
            _get_coords_array(frame).tolist(),
            is_lambda1=is_lambda1,
            map=property_map,
        )

        return _System(frame)

    def saveReplicas(self, filenames, save_velocities=False, is_lambda1=False):
        """
        Save each replica (coordinate set) to individual files.

        Parameters
        ----------

        filenames : list of str
            A list of filenames to save each replica to. The extension will
            determine the format, which must be the same for all replicas.

        save_velocities : bool, optional
            Whether to save velocities along with coordinates. Default is False.
        """

        from sire import save as _save
        from .. import _isVerbose

        if not isinstance(filenames, list):
            raise TypeError("'filenames' must be a list of strings.")
        if not all(isinstance(f, str) for f in filenames):
            raise TypeError("All elements in 'filenames' must be of type 'str'.")

        if len(filenames) != self.nReplicas():
            raise ValueError("'filenames' length must match the number of replicas.")

        if not isinstance(save_velocities, bool):
            raise TypeError("'save_velocities' must be of type 'bool'.")

        try:
            _save(
                self._trajectory,
                filenames,
                save_velocities=save_velocities,
                silent=True,
            )
        except Exception as e:
            msg = "Failed to save replicas to individual files."

            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg)

    @staticmethod
    def loadReplicas(replica_system, filenames):
        """
        Load multiple replicas (coordinate sets) from individual files.

        Parameters
        ----------

        replica_system : :class:`ReplicaSystem <BioSimSpace._SireWrappers.ReplicaSystem>`
            A ReplicaSystem object to use as a template for loading the replicas.

        filenames : list of str
            A list of filenames to load each replica from. The extension will
            determine the format, which must be the same for all replicas.

        Returns
        -------

        :class:`ReplicaSystem <BioSimSpace._SireWrappers.ReplicaSystem>`
            The loaded replica system.
        """

        import os as _os
        from tempfile import NamedTemporaryFile as _NamedTemporaryFile
        from sire import load as _load
        from sire import save as _save
        from .. import _isVerbose

        if not isinstance(replica_system, ReplicaSystem):
            raise TypeError("'replica_system' must be of type 'ReplicaSystem'.")

        if not isinstance(filenames, list):
            raise TypeError("'filenames' must be a list of strings.")
        if not all(isinstance(f, str) for f in filenames):
            raise TypeError("All elements in 'filenames' must be of type 'str'.")

        if len(filenames) == 0:
            raise ValueError("'filenames' must contain at least one filename.")

        # Get the extension of the first file to determine format.
        _, ext = _os.path.splitext(filenames[0])

        # If this is a GRO file, then write a GROMACS topology.
        if replica_system._is_perturbable or ext == ".gro":
            top_suffix = ".top"
        # Otherwise, use AMBER PRM7.
        else:
            top_suffix = ".prm7"

        # Write a temporary topology file using the first file.
        try:
            tmp_top = _NamedTemporaryFile(delete=False, suffix=top_suffix)
            _save(replica_system._new_sire_object, tmp_top.name, silent=True)
        except Exception as e:
            msg = "Failed to write temporary topology file for replica loading."
            tmp_top.close()

            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg)

        # Prepend the temporary topology file to the list of filenames.
        filenames = [tmp_top.name] + filenames

        # Now try to load all the replicas.
        try:
            sire_system = _load(filenames, silent=True)
        except Exception as e:
            msg = "Failed to load replicas from individual files."

            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg)

        # Clean up the temporary topology file.
        tmp_top.close()
        _os.unlink(tmp_top.name)

        return ReplicaSystem(
            replica_system._new_sire_object, trajectory=sire_system.trajectory()
        )

    def setPerturbed(self):
        """
        Set the system to the perturbed (lambda = 1) state.
        """
        self._is_perturbed = True

    def setReference(self):
        """
        Set the system to the reference (lambda = 0) state.
        """
        self._is_perturbed = False
