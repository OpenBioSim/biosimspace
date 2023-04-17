######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2023
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
.. currentmodule:: BioSimSpace.Protocol

Functions
=========

.. autosummary::
    :toctree: generated/

    protocols
    createProtocol

Classes
=======

.. autosummary::
    :toctree: generated/

    Custom
    Equilibration
    FreeEnergy
    Metadynamics
    Minimisation
    Production
    Steering

Examples
========

Print the list of supported protocols.

.. code-block:: python

   import BioSimSpace as BSS

   print(BSS.Protocol.protocols())

Create a default minimisation protocol and print the number of steps.

.. code-block:: python

   import BioSimSpace as BSS

   protocol = BSS.Protocol.Minimisation()
   print(protocol.getSteps())

The same as above, but instead passing "Minimisation" as an argument to the
:class:`createProtocol <BioSimSpace.Protocol.createProtocol>` function. This
function should be used in any interoperable workflow
:class:`Node <BioSimSpace.Gateway.Node>` where the protocol is specified
as an input requirement by the user.

.. code-block:: python

   import BioSimSpace as BSS

   protocol = BSS.Protocol.createProtocol("minimisation")
   print(protocol.getSteps())

Create an equilibration protocol that heats the system from 0 to 300 Kelvin
while restraining the positions of any backbone atoms.

.. code-block:: python

   import BioSimSpace as BSS

   protocol = BSS.Protocol.Equilibration(
       temperature_start=0*BSS.Units.Temperature.kelvin,
       temperature_end=300*BSS.Units.Temperature.kelvin,
       restraint="backbone"
   )
"""

from ._config import *
from ._custom import *
from ._dummy import Dummy
from ._equilibration import *
from ._free_energy import *
from ._free_energy_equilibration import *
from ._free_energy_minimisation import *
from ._free_energy_mixin import *
from ._metadynamics import *
from ._minimisation import *
from ._position_restraint import *
from ._production import *
from ._steering import *
from ._utils import *
