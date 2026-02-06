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
.. currentmodule:: BioSimSpace.Solvent

Functions
=========

.. autosummary::
    :toctree: generated/

    solvate
    spc
    spce
    tip3p
    tip4p
    tip5p
    opc
    waterModels

Examples
========

Print the list of supported water models.

.. code-block:: python

   import BioSimSpace as BSS

   print(BSS.Solvent.waterModels())

Solvate a molecule in a 5 :class:`nanometer <BioSimSpace.Units.Length.nanometer>`
periodic box of TIP3P water with an ion concentration of 0.1 mol per litre.

.. code-block:: python

   import BioSimSpace as BSS

   # Load a system and extract the first molecule.
   molecule = BSS.IO.readMolecules("amber/ala/*")[0]

   # Solvate the molecule.
   solvated = BSS.Solvent.tip3p(
       molecule=molecule,
       box=3*[5*BSS.Units.Length.nanometer],
       ion_conc=0.1
   )

The same as above, but instead passing "TIP3P" as an argument to the
:class:`solvate <BioSimSpace.Solvent.solvate>` function. This function should
be used in any interoperable workflow :class:`Node <BioSimSpace.Gateway.Node>`
where the water model is specified as an input requirement by the user.

.. code-block:: python

   import BioSimSpace as BSS

   # Load a system and extract the first molecule.
   files = BSS.IO.expand(BSS.tutorialUrl(), ["ala.top", "ala.crd"], ".bz2")
   molecule = BSS.IO.readMolecules(files)[0]

   # Solvate the molecule.
   solvated = BSS.Solvent.solvate(
       "tip3p", molecule=molecule,
       box=3*[5*BSS.Units.Length.nanometer],
       ion_conc=0.1
   )

Solvate the molecule with a shell of at least 2 nanometers of SPC water.

.. code-block:: python

   import BioSimSpace as BSS

   # Load a system and extract the first molecule.
   files = BSS.IO.expand(BSS.tutorialUrl(), ["ala.top", "ala.crd"], ".bz2")
   molecule = BSS.IO.readMolecules(files)[0]

   # Solvate the molecule.
   solvated = BSS.Solvent.spc(
      "tip3p", molecule=molecule,
      shell=2*BSS.Units.Length.nanometer
   )

Create a 50 :class:`angstrom <BioSimSpace.Units.Length.angstrom>` periodic
box of pure SPC/E water.

.. code-block:: python

   import BioSimSpace as BSS

   water = BSS.Solvent.spce(box=3*[50*BSS.Units.Length.angstrom])
"""

import sire as _sr

_sr.use_new_api()
del _sr

from ._solvent import *
