.. _ref_crystal_water:
=============
Crystal water
=============

.. toctree::
   :maxdepth: 2

Water molecules are important for the structure, function, and dynamics of
proteins. When setting up a simulation it may be important to preserve
specific water molecules, often referred to as *crystal waters*. For example,
many structures from the `Protein Data Bank <https://www.rcsb.org/>`__
contain the coordinates of water molecules that are resolvable via X-ray
crystallography. In this example, we will use BioSimSpace to parameterise
and solvate a `tyrosine kinase 2 <https://en.wikipedia.org/wiki/Tyrosine_kinase_2>`__
protein structure containing two crystal waters within its binding site.

Firstly, let's load the PDB structure:

>>> import BioSimSpace as BSS
>>> tyk2_xtal = BSS.IO.readMolecules(
...    BSS.IO.expand(
...        BSS.tutorialUrl(), "tyk2_xtal.pdb"
...    )
...)[0]

This is a single molecule that contains the oxgyen atoms of two crystal waters
at the end of the structure:

>>> for residue in tyk2_xtal.getResidues()[-2:]:
...     for atom in residue:
...         print(atom, atom.coordinates())
...
<BioSimSpace.Atom: name='O', molecule=2, index=4652> (-4.5180 A, -0.0520 A, -15.6070 A)
<BioSimSpace.Atom: name='O', molecule=2, index=4653> (-7.2040 A, -8.5540 A, -20.7520 A)

We can use BioSimSpace to search for and extract the protein and water
components of the system. In this case, the waters are part of a residue
named ``WAT``, so we can use this as our search term:

>>> tyk2 = tyk2_xtal.extract(
...    [atom.index() for atom in tyk2_xtal.search("not resname WAT").atoms()]
... )
>>> xtal_water = tyk2_xtal.extract(
...    [atom.index() for atom in tyk2_xtal.search("resname WAT").atoms()]
... )

Now we have the components, we can parameterise them individually. First the
protein:

>>> tyk2 = BSS.Parameters.ff14SB(tyk2).getMolecule()

And then the water:

>>> xtal_water = BSS.Parameters.ff14SB(
...     xtal_water,
...     water_model="tip3p",
...     ensure_compatible=False
... ).getMolecule()

In this case we need to specify the desired water model, and also that we don't
want to ensure that the parameterised molecule is compatible with the topology
of the molecule that we passed in. This is because we are only passing in
oxygen atoms, so ``tLEaP`` will add the missing hydrogens.

.. note ::

   By default, BioSimSpace ensures that the topology of the parameterised molecule(s)
   matches that of the input molecule and will raise an exception if this is not
   the case. This is because the input system is often required as a reference by
   the user, e.g. they might want to preserve a specific naming convention, chain
   identifiers, etc.

Let's check one of the parameterised crystal waters:

>>> for atom in xtal_water[0].getAtoms():
...     print(atom, atom.coordinates())
...
<BioSimSpace.Atom: name='O', molecule=9, index=0> (-4.5180 A, -0.0520 A, -15.6070 A)
<BioSimSpace.Atom: name='H1', molecule=9, index=1> (-3.5608 A, -0.0520 A, -15.6070 A)
<BioSimSpace.Atom: name='H2', molecule=9, index=2> (-4.7580 A, 0.8746 A, -15.6070 A)

Note that the coordinates of the oxygen atom are preserved.

Once both components are parameterised, we can combine them together to create
a new system:

>>> system = tyk2.toSystem() + xtal_water

The next step in our setup procedure is to solvate the system in a water box.
Here we will use a truncated octahedral box with a base length of 7 nanometers:

>>> box, angles = BSS.Box.truncatedOctahedron(7 * BSS.Units.Length.nanometer)
>>> solvated = BSS.Solvent.tip3p(system, box=box, angles=angles)

Since the protein is charged, during the solvation process ``gmx genion`` will
have been used to neutralise the solvated system by adding counter ions. We can
check this as follows:

>>> print(system.charge(), solvated.charge())
... -2.0000 |e| -1.5091e-07 |e|
>>> print(len(solvated.search("element Na")))
... 2

In this case, two sodium ions have been added to neutralise the system. In
doing so, ``gmx genion`` will have removed two *random* water molecules. In
order to ensure that crystal waters are not removed, we temporarily tag them
with a unique residue and molecule name during solvation.

To check that they are preserved we can re-solvate the system, asking to
preserve the water naming for the existing waters in the system. (By default,
they will be updated to the default ``GROMACS`` naming convention.)

>>> solvated_no_match = BSS.Solvent.tip3p(
...     system,
...     box=box,
...     angles=angles,
...     match_water=False
... )

We can check that the crystal waters are still present as follows:

>>> for residue in solvated_no_match.getWaterMolecules().toSystem().search("resname WAT").residues():
...     print(residue)
<BioSimSpace.Residue: name='WAT', molecule=8, index=0, nAtoms=3>
<BioSimSpace.Residue: name='WAT', molecule=9, index=0, nAtoms=3>

.. note ::

   When solvating, molecules in the original system will be centered within the
   solvation box, hence the coordinates of the crystal waters after solvation
   won't necessarily match those from before. By setting ``match_water=False``,
   it is possible to preserve the naming convention used for any existing
   waters in the system. However, note that some simulation engines rely on a
   specific naming convention for water molecules and a single solvent group.
