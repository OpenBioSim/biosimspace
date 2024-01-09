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
.. currentmodule:: BioSimSpace.Parameters

Functions
=========

.. autosummary::
    :toctree: generated/

    __FORCE_FIELDS__
    formalCharge

Examples
========

Print the list of supported force fields.

.. code-block:: python

   import BioSimSpace as BSS

   print(BSS.Parameters.forceFields())

Parameterise a molecule using GAFF.

.. code-block:: python

   import BioSimSpace as BSS

   # Load a molecule from file.
   molecule = BSS.IO.readMolecules("molecules/benzene.pdb")

   # Initialise the parameterisation process. This will run the parameterisation
   # in the background. This is useful when you are working interactively and
   # wish to continue doing other things while the parameterisation runs.
   # (Parameterisation can be slow.)
   process = BSS.Parameters.gaff(molecule)

   # Get the parameterised molecule. This will now block until the
   # parameterisation is finished.
   molecule = process.getMolecule()

The same as above, but immediately getting the parameterised molecule from
the process object, i.e. not saving it as an intermediate variable.

.. code-block:: python

   import BioSimSpace as BSS

   # Load a molecule from file.
   molecule = BSS.IO.readMolecules("molecules/benzene.pdb")

   # Initialise the parameterisation process and block until the molecule is
   # ready to be returned.
   molecule = BSS.Parameters.gaff(molecule).getMolecule()

The same as above, but instead passing "GAFF" as an argument to the
:class:`parameterise <BioSimSpace.Parameters.parameterise>` function. This
function should be used in any interoperable workflow
:class:`Node <BioSimSpace.Gateway.Node>` where the force field is specified
as an input requirement by the user.

.. code-block:: python

   import BioSimSpace as BSS

   # Load a molecule from file.
   molecule = BSS.IO.readMolecules("molecules/benzene.pdb")

   # Initialise the parameterisation process and block until the molecule is
   # ready to be returned.
   molecule = BSS.Parameters.parameterise(molecule, "gaff").getMolecule()

Parameterise a molecule using GAFF, passing the net formal charge computed
by BioSimSpace. This is useful when your input PDB file has missing or
incorrect formal charge information.

.. code-block:: python

   import BioSimSpace as BSS

   # Load a molecule from file.
   molecule = BSS.IO.readMolecules("molecules/benzene.pdb")

   # Compute the net formal charge on the molecule.
   formal_charge = BSS.Parameters.formalCharge(molecule)

   # Initialise the parameterisation process. This will run the parameterisation
   # in the background. This is useful when you are working interactively and
   # wish to continue doing other things while the parameterisation runs.
   # (Parameterisation can be slow.)
   process = BSS.Parameters.gaff(molecule, net_charge=formal_charge)

   # Get the parameterised molecule. This will now block until the
   # parameterisation is finished.
   molecule = process.getMolecule()

Parameterise a molecule using GAFF, passing the molecule as a SMILES
string.

.. code-block:: python

   import BioSimSpace as BSS

   # Parameterise directly using a SMILES string.
   molecule = BSS.Parameters.gaff("C1=CC=CC=C1").getMolecule()

When parameterising a molecule containing a structural ion with an AMBER force
field then it is necessary to pass the name of a water model that will be used
for the ion parameters, e.g:

.. code-block:: python

   import BioSimSpace as BSS

   # Load a molecule from file.
   molecule = BSS.IO.readMolecules("molecules/4V2Y_A.pdb")

   # Initialise the parameterisation process and block until the molecule is
   # ready to be returned.
   molecule = BSS.Parameters.ff14SB(molecule, water_model="tip3p").getMolecule()

Additional parameters can be loaded by passing the ``pre_mol_commands`` option to
any compatible AMBER force field function, e.g.:

.. code-block:: python

   import BioSimSpace as BSS

   # Load a molecule from file.
   molecule = BSS.IO.readMolecules("molecules/peptide.pdb")

   # Create a list of the additional commands for the LEaP program. These will
   # be run after any default commands.
   leap_commands = ["addPath phosphate_params", "source leaprc.phosaa10"]

   # Initialise the parameterisation process and block until the molecule is
   # ready to be returned.
   molecule = BSS.Parameters.ff14SB(molecule, pre_mol_commands=leap_commands).getMolecule()

If you need to apply any additional operations on the loaded molecule, then this
can be done using the ``post_mol_commands`` option, using ``mol`` for the name of
the molecule.
"""

from ._parameters import *
from ._utils import *
