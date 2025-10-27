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
.. currentmodule:: BioSimSpace.Process

Functions
=========

.. autosummary::
    :toctree: generated/

    engines
    createProcess

MD driver classes
=================

.. autosummary::
    :toctree: generated/

    Amber
    Gromacs
    Namd
    OpenMM
    Plumed
    Somd

Multi process simulation tools
==============================

.. autosummary::
    :toctree: generated/

    ProcessRunner

Examples
========

Print the list of supported molecular dynamics engines.

.. code-block:: python

   import BioSimSpace as BSS

   print(BSS.Process.engines())

Create a process to apply a minimisation protocol to a molecular system
using the `AMBER <http://ambermd.org>`_ package. Execute the process and
get the minimised molecular system.

.. code-block:: python

   import BioSimSpace as BSS

   # Load a molecular system from file.
   files = BSS.IO.expand(BSS.tutorialUrl(), ["ala.top", "ala.crd"], ".bz2")
   system = BSS.IO.readMolecules(files)

   # Create a minimisation protocol with 1000 steps.
   protocol = BSS.Protocol.Minimisation(steps=1000)

   # Create a process object to run the simulation with AMBER.
   process = BSS.Process.Amber(system, protocol)

   # Start the process in the background.
   process.start()

   # Wait for the process to finish.
   process.wait()

   # Get the minimised molecular system.
   minimised = process.getSystem()
"""

import sire as _sr

_sr.use_new_api()
del _sr

from ._amber import *
from ._gromacs import *
from ._namd import *
from ._openmm import *
from ._plumed import *
from ._process_runner import *
from ._somd import *
from ._utils import *
