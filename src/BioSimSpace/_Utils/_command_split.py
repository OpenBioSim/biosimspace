######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2025
#
# Authors: Christopher Woods <chryswoods@hey.com>
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

__author__ = "Christopher Woods"
__email__ = "chryswoods@hey.com"

__all__ = ["command_split"]


def command_split(command):
    """
    Cross platform version of 'shlex.split'. This will split the passed
    command into parts, doing the right thing on Linux, MacOS and Windows.
    """
    import sys

    if sys.platform != "win32":
        # We can just use shlex.split - it is only windows that is annoying!
        import shlex

        return shlex.split(command)

    # thanks for inspiration to kxr on this stackoverflow post
    # https://stackoverflow.com/questions/33560364/python-windows-parsing-command-lines-with-shlex

    import re

    regex = r""""((?:""|\\["\\]|[^"])*)"?()|(\\\\(?=\\*")|\\")|(&&?|\|\|?|\d?>|[<])|([^\s"&|<>]+)|(\s+)|(.)"""

    args = []
    accumulated = None

    for qs, qss, esc, pipe, word, white, fail in re.findall(regex, command):
        if word:
            pass  # most frequent
        elif esc:
            word = esc[1]
        elif white or pipe:
            if accumulated is not None:
                args.append(accumulated)
            if pipe:
                args.append(pipe)
            accumulated = None
            continue
        elif fail:
            raise ValueError("invalid or incomplete shell string")
        elif qs:
            word = qs.replace('\\"', '"').replace("\\\\", "\\")
            word = word.replace('""', '"')
        else:
            word = qss  # may be even empty; must be last

        accumulated = (accumulated or "") + word

    if accumulated is not None:
        args.append(accumulated)

    return args
