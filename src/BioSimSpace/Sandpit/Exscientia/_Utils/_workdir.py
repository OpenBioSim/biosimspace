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

"""Custom context managers."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["WorkDir"]


class WorkDir:
    """A utility class to create a working directory."""

    def __init__(self, work_dir=None):
        """
        Constructor.

        Parameters
        ----------

        work_dir : str
            The working directory for the context. If None, then a temporary
            working directory will be created.
        """
        import os as _os
        import tempfile as _tempfile

        # Validate the input.
        if work_dir and not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'")

        # Temporary directory.
        if work_dir is None:
            self._tmp_dir = _tempfile.TemporaryDirectory()
            self._work_dir = self._tmp_dir.name

        # User specified directory.
        else:
            self._tmp_dir = None

            # Use absolute path.
            if not _os.path.isabs(work_dir):
                work_dir = _os.path.abspath(work_dir)

            # Create the directory if it doesn't already exist.
            if not _os.path.isdir(work_dir):
                _os.makedirs(work_dir, exist_ok=True)

            self._work_dir = work_dir

    def __str__(self):
        return self._work_dir

    def __repr__(self):
        return self._work_dir

    def __add__(self, other):
        if not isinstance(other, str):
            raise TypeError(
                f"unsupported operand type(s) for +: 'WorkDir' and '{type(other)}'"
            )
        return str(self._work_dir) + other

    def is_temp_dir(self):
        """Whether this is a temporary directory."""
        return self._tmp_dir is not None
