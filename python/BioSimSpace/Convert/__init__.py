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
.. currentmodule:: BioSimSpace.Convert

Functions
=========

.. autosummary::
    :toctree: generated/

    smiles
    supportedFormats
    to
    toBioSimSpace
    toOpenMM
    toRDKit
    toSire

Examples
========

Load a system and convert to various formats.

.. code-block:: python

   import BioSimSpace as BSS

   files = BSS.IO.expand(
       BSS.tutorialUrl(),
       ["ala.crd", "ala.top"],
       ".bz2"
   )
   system = BSS.IO.readMolecules(files)

   # Convert to Sire format.
   sire_system = BSS.Convert.to(system, "sire")

   # Convert to RDKit format.
   rdkit_mols = BSS.Convert.to(system, "rdkit")

   # Convert a molecule to Sire format.
   sire_mol = BSS.Convert.to(system[0], "sire")

   # Convert a molecule to RDKit format.
   rdkit_mol = BSS.Convert.to(system[0], "rdkit")

   # Convert a residue to Sire format.
   sire_res = BSS.Convert.to(system[0].getResidues()[0], "sire")

   # Convert a residue to RDKit format. (This will return a single
   # residue RDMol.)
   rdkit_res = BSS.Convert.to(system[0].getResidues()[0], "rdkit")

   # Convert an atom to Sire format.
   sire_atom = BSS.Convert.to(system[0].getAtoms()[0], "sire")

   # Convert an atom to RDKit format. (This will return a single
   # atom RDMol.)
   sire_atom = BSS.Convert.to(system[0].getAtoms()[0], "rdkit")

   # Conversions can be chained, for example.
   bss_mols = BSS.Convert.to(
       BSS.Convert.to(
           BSS.Convert.to(system, "sire"),
           "rdkit"
       ),
       "biosimspace"
   )
"""

from ._convert import *
