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

"""Incremental reading of files that are being appended to."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Tail"]


class Tail:
    """
    Read a file incrementally, yielding only the lines appended since the
    previous read. Only complete (newline-terminated) lines are returned;
    a partial trailing line is left for the next read.
    """

    def __init__(self, filename):
        """
        Constructor.

        Parameters
        ----------

        filename : str
            The path to the file.
        """
        if not isinstance(filename, str):
            raise TypeError("'filename' must be of type 'str'")

        self._filename = filename
        self._offset = 0
        self._inode = None

    def __iter__(self):
        import os as _os

        try:
            stat = _os.stat(self._filename)
        except FileNotFoundError:
            return

        # Start again if the file has been replaced or truncated.
        if stat.st_ino != self._inode or stat.st_size < self._offset:
            self._inode = stat.st_ino
            self._offset = 0

        with open(self._filename, "rb") as file:
            file.seek(self._offset)

            while True:
                line = file.readline()

                if not line.endswith(b"\n"):
                    break

                self._offset = file.tell()

                # Normalise Windows line endings.
                if line.endswith(b"\r\n"):
                    line = line[:-2] + b"\n"

                yield line.decode(errors="replace")
