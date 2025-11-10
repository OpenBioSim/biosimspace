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
.. currentmodule:: BioSimSpace.Node

Functions
=========

.. autosummary::
    :toctree: generated/

    help
    list
    run
    setNodeDirectory

Examples
========

Print the list of available nodes.

.. code-block:: python

   import BioSimSpace as BSS

   print(BSS.Node.list())

Get help on the "minimisation" node.

.. code-block:: python

   import BioSimSpace as BSS

   BSS.Node.help("minimisation")

Run the minimisation node.

.. code-block:: python

   import BioSimSpace as BSS

   # Generate a dictionary of input arguments.
   input = {"steps" : 1000, "files" : ["amber/ala/ala.top", "amber/ala/ala.crd"]}

   # Run the node and get the output as a dictionary.
   output = BSS.Node.run("minimisation", input)

Set a custom directory for the node library.

.. code-block:: python

   import BioSimSpace as BSS

   BSS.Node.setNodeDirectory("/path/to/node/library")
"""

import sire as _sr

_sr.use_new_api()
del _sr

from ._node import *
