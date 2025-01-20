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

__all__ = ["cd"]

from contextlib import contextmanager as _contextmanager

import os as _os

from ._workdir import WorkDir as _WorkDir


# Adapted from: http://ralsina.me/weblog/posts/BB963.html
@_contextmanager
def cd(work_dir):
    """
    Execute the context in the directory "work_dir".

    Parameters
    ----------

    work_dir : str
        The working directory for the context.
    """

    # Validate the input.
    if not isinstance(work_dir, (str, _WorkDir)):
        raise TypeError(
            "'work_dir' must be of type 'str' or 'BioSimSpace._Utils.WorkDir'"
        )

    # Get path as a string.
    if isinstance(work_dir, _WorkDir):
        work_dir = str(work_dir)

    # Store the current directory.
    old_dir = _os.getcwd()

    # Create the working directory if it doesn't exist.
    if not _os.path.isdir(work_dir):
        work_dir = str(_WorkDir(work_dir))

    # Change to the new directory.
    _os.chdir(work_dir)

    # Execute the context.
    try:
        yield

    # Return to original directory.
    finally:
        _os.chdir(old_dir)
