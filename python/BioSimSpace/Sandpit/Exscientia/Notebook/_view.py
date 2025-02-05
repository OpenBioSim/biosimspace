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

"""Tools for visualising molecular systems."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["View"]

import glob as _glob
import os as _os
import shutil as _shutil
import tempfile as _tempfile
import warnings as _warnings

from sire.legacy import IO as _SireIO
from sire.legacy import Mol as _SireMol
from sire.legacy import System as _SireSystem

from .. import _is_notebook, _isVerbose
from .. import IO as _IO
from ..Process._process import Process as _Process
from .._SireWrappers import System as _System


class View:
    """A class for handling interactive molecular visualisations."""

    def __init__(self, handle, property_map={}, is_lambda1=False):
        """
        Constructor.

        Parameters
        ----------

        handle : :class:`Process <BioSimSpace.Process>`, \
                 :class:`System <BioSimSpace._SireWrappers.System>` \
                 :class:`System <BioSimSpace._SireWrappers.Molecule>` \
                 :class:`System <BioSimSpace._SireWrappers.Molecules>` \
                 str, [str]
            A handle to a process, system, molecule, or molecule container,
            or the path to molecular input file(s).

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        is_lambda1 : bool
            Whether to use the lambda = 1 end state when visualising
            perturbable molecules. By default, the state at lambda = 0
            is used.
        """

        # Make sure we're running inside a Jupyter notebook.
        if not _is_notebook:
            _warnings.warn(
                "You can only use BioSimSpace.Notebook.View from within a Jupyter notebook."
            )
            return None

        # Validate the map.
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")
        self._property_map = property_map

        # Check the end state flag.
        if not isinstance(is_lambda1, bool):
            raise TypeError("'is_lambda1' must be of type 'bool'")
        self._is_lamba1 = is_lambda1

        # Check the handle.

        # Convert tuple to list.
        if isinstance(handle, tuple):
            handle = list(handle)

        # Convert single string to list.
        if isinstance(handle, str):
            handle = [handle]

        # List of strings (file paths).
        if isinstance(handle, list) and all(isinstance(x, str) for x in handle):
            system = _IO.readMolecules(handle, property_map=property_map)
            self._handle = system._getSireObject()
            self._is_process = False

        # BioSimSpace process.
        elif isinstance(handle, _Process):
            self._handle = handle
            self._is_process = True

        # BioSimSpace system.
        elif isinstance(handle, _System):
            self._handle = handle.copy()._getSireObject()
            self._is_process = False

        else:
            try:
                handle = handle.toSystem()
                self._handle = handle._getSireObject()
                self._is_process = False
            except:
                raise TypeError(
                    "The handle must be of type 'BioSimSpace.Process', "
                    "'BioSimSpace._SireWrappers.System', "
                    "'BioSimSpace._SireWrappers.Molecule', "
                    "'BioSimSpace._SireWrappers.Molecules', "
                    "'str', or a list of 'str' types."
                )

        # Create a temporary workspace for the view object.
        self._tmp_dir = _tempfile.TemporaryDirectory()
        self._work_dir = self._tmp_dir.name

        # Zero the number of views.
        self._num_views = 0

        # Reconstruct a system representing the chosen lambda end state.
        if not self._is_process:
            self._handle = self._reconstruct_system(self._handle, self._is_lamba1)

    def system(self, gui=True):
        """
        View the entire molecular system.

        Parameters
        ----------

        gui : bool
            Whether to display the gui.
        """

        # Make sure we're running inside a Jupyter notebook.
        if not _is_notebook:
            return None

        # Get the latest system from the process.
        if self._is_process:
            system = self._handle.getSystem()

            # No system.
            if system is None:
                return
            else:
                system = system._getSireObject()

            # Reconstruct the chosen lambda end state if perturbable molecules
            # are present.
            system = self._reconstruct_system(system, self._is_lamba1)

        else:
            system = self._handle

        # Create and return the view.
        return self._create_view(system, gui=gui)

    def molecules(self, indices=None, gui=True):
        """
        View specific molecules.

        Parameters
        ----------

        indices : [int], range
            A list of molecule indices.

        gui : bool
            Whether to display the gui.
        """

        # Make sure we're running inside a Jupyter notebook.
        if not _is_notebook:
            return None

        # Return a view of the entire system.
        if indices is None:
            return self.system(gui=gui)

        # Convert range or tuple to list.
        if isinstance(indices, (range, tuple)):
            indices = list(indices)
        # Convert single indices to a list.
        elif not isinstance(indices, list):
            indices = [indices]

        # Check that the indices is a list of integers.
        if not all(type(x) is int for x in indices):
            raise TypeError("'indices' must be a 'list' of type 'int'")

        # Get the latest system from the process.
        if self._is_process:
            system = self._handle.getSystem()._getSireObject()

            # No system.
            if system is None:
                return

            # Reconstruct the chosen lambda end state if perturbable molecules
            # are present.
            system = self._reconstruct_system(system, self._is_lamba1)

        else:
            system = self._handle

        # Extract the molecule numbers.
        molnums = system.molNums()

        # Create a new system.
        s = _SireSystem.System("BioSimSpace_System")
        m = _SireMol.MoleculeGroup("all")

        # Loop over all of the indices.
        for index in indices:
            if index < 0:
                index += len(molnums)
            if index < 0 or index > len(molnums):
                raise ValueError("Molecule index is out of range!")

            # Add the molecule.
            m.add(system[molnums[index]])

        # Add all of the molecules to the system.
        s.add(m)

        # Create and return the view.
        return self._create_view(s, gui=gui)

    def molecule(self, index=0, gui=True):
        """
        View a specific molecule.

        Parameters
        ----------

        index : int
            The molecule index.

        gui : bool
            Whether to display the gui.
        """

        # Make sure we're running inside a Jupyter notebook.
        if not _is_notebook:
            return None

        # Check that the index is an integer.
        if not type(index) is int:
            raise TypeError("'index' must be of type 'int'")

        # Get the latest system from the process.
        if self._is_process:
            system = self._handle.getSystem()._getSireObject()

            # No system.
            if system is None:
                return

            # Reconstruct the chosen lambda end state if perturbable molecules
            # are present.
            system = self._reconstruct_system(system, self._is_lamba1)

        else:
            system = self._handle

        # Extract the molecule numbers.
        molnums = system.molNums()

        # Make sure the index is valid.
        if index < 0:
            index += len(molnums)
        if index < 0 or index > len(molnums):
            raise ValueError("Molecule index is out of range!")

        # Create a new system and add a single molecule.
        s = _SireSystem.System("BioSimSpace_System")
        m = _SireMol.MoleculeGroup("all")
        m.add(system[molnums[index]])
        s.add(m)

        # Create and return the view.
        return self._create_view(s, gui=gui)

    def reload(self, index=None, gui=True):
        """
        Reload a particular view.

        Parameters
        ----------

        index : int
            The view index.

        gui : bool
            Whether to display the gui.
        """

        # Make sure we're running inside a Jupyter notebook.
        if not _is_notebook:
            return None

        # Return if there are no views.
        if self._num_views == 0:
            return

        # Default to the most recent view.
        if index is None:
            index = self._num_views - 1

        # Check that the index is an integer.
        elif not type(index) is int:
            raise TypeError("'index' must be of type 'int'")

        # Make sure the view index is valid.
        if index < 0 or index >= self._num_views:
            raise ValueError(
                "View index (%d) is out of range: [0-%d]" % (index, self._num_views - 1)
            )

        # Create and return the view.
        return self._create_view(view=index, gui=gui)

    def nViews(self):
        """
        Return the number of views.

        Return
        ------

        num_views : int
            The number of views.
        """
        return self._num_views

    def savePDB(self, file, index=None):
        """
        Save a specific view as a PDB file.

        Parameters
        ----------

        file : str
            The name of the file to write to.

        index : int
            The view index.
        """

        # Make sure we're running inside a Jupyter notebook.
        if not _is_notebook:
            return None

        # Default to the most recent view.
        if index is None:
            index = self._num_views - 1

        # Check that the index is an integer.
        elif not type(index) is int:
            raise TypeError("'index' must be of type 'int'")

        # Make sure the view index is valid.
        if index < 0 or index >= self._num_views:
            raise ValueError(
                "View index (%d) is out of range: [0-%d]" % (index, self._num_views - 1)
            )

        # Copy the file to the chosen location.
        _shutil.copyfile("%s/view_%04d.pdb" % (self._work_dir, index), file)

    def reset(self):
        """Reset the object, clearing all view files."""

        # Glob all of the view PDB structure files.
        files = _glob.glob("%s/*.pdb" % self._work_dir)

        # Remove all of the files.
        for file in files:
            _os.remove(file)

        # Reset the number of views.
        self._num_views = 0

    def _create_view(self, system=None, view=None, gui=True):
        """
        Helper function to create the NGLview object.

        Parameters
        ----------

        system : Sire.System.System
            A Sire molecular system.

        view : int
            The index of an existing view.

        gui : bool
            Whether to display the gui.
        """

        if system is None and view is None:
            raise ValueError("Both 'system' and 'view' cannot be 'None'.")

        elif system is not None and view is not None:
            raise ValueError("One of 'system' or 'view' must be 'None'.")

        # Make sure gui flag is valid.
        if gui not in [True, False]:
            gui = True

        # Default to the most recent view.
        if view is None:
            index = self._num_views
        else:
            index = view

        # Create the file name.
        filename = "%s/view_%04d.pdb" % (self._work_dir, index)

        # Increment the number of views.
        if view is None:
            self._num_views += 1

        # Create a PDB object and write to file.
        if system is not None:
            try:
                pdb = _SireIO.PDB2(system, self._property_map)
                pdb.writeToFile(filename)
            except Exception as e:
                msg = "Failed to write system to 'PDB' format."
                if _isVerbose():
                    print(msg)
                    raise IOError(e) from None
                else:
                    raise IOError(msg) from None

        # Import NGLView when it is used for the first time.
        import nglview as _nglview

        # Create the NGLview object.
        view = _nglview.show_file(filename)

        # Return the view and display it.
        return view.display(gui=gui)

    def _reconstruct_system(self, system, is_lambda1=False):
        """
        Helper function to reconstruct a lambda end state for a perturbable
        system.

        Parameters
        ----------

        system : Sire.System.System
            A Sire molecular system.

        is_lambda1 : bool
            Whether to use the lambda = 1 end state for reconstructing
            perturbable molecules. By default, the state at lambda = 0
        """

        # Convert to a BioSimSpace system.
        system = _System(system)

        # Get a list of perturbable molecules.
        pert_mols = system.getPerturbableMolecules()

        # If there are no perturbable molecules, simply return the original system.
        if len(pert_mols) == 0:
            return system._sire_object

        # Convert all perturbable molecules to the chosen end state, updating
        # them in the system.
        else:
            for mol in pert_mols:
                system.updateMolecules(mol._toRegularMolecule(is_lambda1=is_lambda1))

        return system._sire_object
