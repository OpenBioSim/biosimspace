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
.. currentmodule:: BioSimSpace.MD

Functions
=========

.. autosummary::
    :toctree: generated/

    engines
    run

Examples
========

.. code-block:: python

   import BioSimSpace as BSS

   # Load a molecular system from file.
   system = BSS.IO.readMolecules(
       BSS.IO.expand(
           BSS.tutorialUrl(),
           ["ala.top, "ala.crd"],
           ".bz2")
       )
   )

   # Create a default minimisation protocol.
   protocol = BSS.Protocol.Minimisation()

   # Find a molecular dynamics package on the host system that supports the
   # system and protocol defined above. If a package exists, BioSimSpace
   # will auto-generate all of the required input files and return a handle
   # to a background process that is running the simulation.
   process = BSS.MD.run(system, protocol)
"""

from ._md import *
from ._utils import *
