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
.. currentmodule:: BioSimSpace.Box

Functions
=========

.. autosummary::
    :toctree: generated/

    boxTypes
    generateBoxParameters
    cubic
    rhombicDodecahedronSquare
    rhombicDodecahedronHexagon
    truncatedOctahedron

Examples
========

Generate the lattice vectors and angles for a truncated octahedron
of 10 nanometer lattice separation.

.. code-block:: python

   import BioSimSpace as BSS

   box, angles = BSS.Box.truncatedOctahedron(10 * BSS.Units.Length.nanometer)
"""

import sire as _sr

_sr.use_new_api()
del _sr

from ._box import *
