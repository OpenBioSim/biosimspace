=================
Molecular Editing
=================

In this tutorial, you will use BioSimSpace and Sire to edit force field parameters
a small molecule. Sire is a molecular simulation framework that BioSimSpace is
built upond, and provides a wide range of functionality for manipulating and editing
molecules.

--------------------------
Editing existing molecules
--------------------------

In this section you will lean how to edit the force field parameters of an existing
molecule, e.g. if you want to adjust charges, or modify bonded parameters.

First, let's import the necessary modules and use BioSimSpace to create a
a parameterised ethanol molecule from a SMILES string:

>>> import BioSimSpace as BSS
>>> import sire as sr
>>> ethanol = BSS.Parameters.gaff("CCO").getMolecule()

Forcefield parameters are stored as molecular *properties*. These can be
accessed by using the appropriate property *key* on the underlying Sire molecule,
e.g. for the charges on the atoms.

>> print(ethanol._sire_object.property("charge"))
SireMol::AtomCharges( size=9
0: -0.1361 |e|
1: 0.1264 |e|
2: -0.5998 |e|
3: 0.0423669 |e|
4: 0.0423669 |e|
5: 0.0423669 |e|
6: 0.0431999 |e|
7: 0.0431999 |e|
8: 0.396 |e|
)

Bonded terms are stored internall as expressions using a buit-in computer algebra
system:

>>> for bond in ethanol._sire_object.property("bond").potentials():
...     print(bond)
TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(3)} : 330.6 [r - 1.0969]^2 )
TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(1)} : 300.9 [r - 1.5375]^2 )
TwoAtomFunction( {CGIdx(0),Index(1)} <-> {CGIdx(0),Index(2)} : 316.7 [r - 1.4233]^2 )
TwoAtomFunction( {CGIdx(0),Index(1)} <-> {CGIdx(0),Index(7)} : 330.6 [r - 1.0969]^2 )
TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(5)} : 330.6 [r - 1.0969]^2 )
TwoAtomFunction( {CGIdx(0),Index(1)} <-> {CGIdx(0),Index(6)} : 330.6 [r - 1.0969]^2 )
TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(4)} : 330.6 [r - 1.0969]^2 )
TwoAtomFunction( {CGIdx(0),Index(2)} <-> {CGIdx(0),Index(8)} : 371.4 [r - 0.973]^2 )

Each ``TwoAtomFunction`` stores the two atoms involved in the bond, as well as the
expression used to compute the bonded energy:

>>> print(bond.atom0(), bond.atom1(), bond.function())
{CGIdx(0),Index(2)} {CGIdx(0),Index(8)} 371.4 [r - 0.973]^2

On write, these expressions are converted into the appropriate functional form
for the parser in question, e.g. for conversion to AMBER format:

>>> from sire.legacy.CAS import Symbol
>>> from sire.legacy.MM import AmberBond
>>> amber_bond = AmberBond(bond.function(), Symbol("r"))

We can now print the components of the bond function, e.g. the force constant and
the equilibrium bond length:

>>> print(amber_bond.k(), amber_bond.r0())
371.4 0.973

To go the other way, we can convert the ``AmberBond`` back into an ``Expression``:

>>> print(amber_bond.asExpression(Symbol("r")))
371.4 [r - 0.973]^2

To edit *atomic* properties, such as charges, we can use the *new* Sire Python API:

>>> cursor = ethanol._sire_object.cursor()
>>> for atom in cursor.atoms():
...     atom["charge"] = 10 * sr.units.mod_electron
>>> ethanol._sire_object = cursor.commit()
>>> print(ethanol._sire_object.property("charge"))
SireMol::AtomCharges( size=9
0: 10 |e|
1: 10 |e|
2: 10 |e|
3: 10 |e|
4: 10 |e|
5: 10 |e|
6: 10 |e|
7: 10 |e|
8: 10 |e|
)

For *molecule* properties, such as bonded parameters, we currently need to use
the legacy Sire API. Here we will create a new set of bonds by nulling the force
constant for all existing bonds:

>>> from sire.legacy.MM import TwoAtomFunctions
>>> bonds = TwoAtomFunctions(ethanol._sire_object.info())
>>> for bond in ethanol._sire_object.property("bond").potentials():
...     amber_bond = AmberBond(0, 0)
...     bonds.set(bond.atom0(), bond.atom1(), amber_bond.toExpression(Symbol("r")))
>>> cursor = ethanol._sire_object.cursor()
>>> cursor["bond"] = bonds
>>> ethanol._sire_object = cursor.commit()

Now lets print the new bonds:

>>> for bond in ethanol._sire_object.property("bond").potentials():
...     print(bond)
TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(3)} : 0 )
TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(1)} : 0 )
TwoAtomFunction( {CGIdx(0),Index(1)} <-> {CGIdx(0),Index(2)} : 0 )
TwoAtomFunction( {CGIdx(0),Index(1)} <-> {CGIdx(0),Index(7)} : 0 )
TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(5)} : 0 )
TwoAtomFunction( {CGIdx(0),Index(1)} <-> {CGIdx(0),Index(6)} : 0 )
TwoAtomFunction( {CGIdx(0),Index(0)} <-> {CGIdx(0),Index(4)} : 0 )
TwoAtomFunction( {CGIdx(0),Index(2)} <-> {CGIdx(0),Index(8)} : 0 )

.. Note::
   Here we've adjusted the parameters for the existing bonding, but we could
   also adjust which atoms are bonded too.

----------------------
Creating new molecules
----------------------

In this section you will learn how to create a new molecule from scratch. For
simplicity, we will assume that we have a reference molecule to serve as a template,
i.e. we already have a set of molecular properties to apply. In practice, you might
want to set the properties manually.

First we need to build the topology of our molecule. To do so, we will create
a new molecule and build the structure of residues and atoms from our ethanol
template. Sire works with residue-based cutting groups, so we create a new
cut-group for each residue in our template, then add atoms to them:

>> mol = sr.legacy.Mol.Molecule("ethanol")
>>> for i, res in enumerate(ethanol._sire_object.residues()):
...     new_res = mol.edit().add(sr.legacy.Mol.ResNum(i+1))
...     new_res.rename(res.name())
...     cg = new_res.molecule().add(sr.legacy.Mol.CGName(f{"i}"))
...     for j, atom in enumerate(res.atoms()):
...         new_atom = cg.add(atom.name())
...         new_atom.renumber(sr.legacy.Mol.AtomNum(j+1))
...         new_atom.reparent(sr.legacy.Mol.ResIdx(i))
...     mol = cg.molecule().commit()

We now have the basic structure of our molecule:

>>> for res in mol.residues():
...     print(res)
...     for atom in res.atoms():
...         print(f"  {atom}")
Residue( LIG:1   num_atoms=9 )
  Atom( C1:1 )
  Atom( C2:2 )
  Atom( O3:3 )
  Atom( H4:4 )
  Atom( H5:5 )
  Atom( H6:6 )
  Atom( H7:7 )
  Atom( H8:8 )
  Atom( H9:9 )

Next we can start adding the required properties, e.g. for the charge:

>>> cursor = mol.cursor()
>>> for new_atom, old_atom in zip(cursor.atoms(), ethanol._sire_object.atoms()):
...     new_atom["charge"] = old_atom["charge"]
>>> mol = cursor.commit()

.. Note::
   Here we've copied the charges from our template, but we could also set them
   manually. In this case we have added them atom-by-atom, but we could also
   copy the entire molecular property in one go.

Let's check that the charges have been added correctly:

>>> print(mol.property("charge"))
SireMol::AtomCharges( size=9
0: 10 |e|
1: 10 |e|
2: 10 |e|
3: 10 |e|
4: 10 |e|
5: 10 |e|
6: 10 |e|
7: 10 |e|
8: 10 |e|
)

Similarly, we can set molecule properties, such as `bond` or `angle`. Here we will
copy the existing `angle` property across:

>>> cursor = mol.cursor()
>>> cursor["angle"] = ethanol._sire_object.property("angle")
>>> mol = cursor.commit()

As before, this could also be done manually, e.g. by creating a new `TheeAtomFunctions`
object and adding angles one-by-one. Let's check that the angles have been added correctly:

>>> print(mol.property("angle"))
ThreeAtomFunctions( size=13
0:    C1:1-C2:2-O3:3       : 67.5 [theta - 1.92318]^2
1:    C1:1-C2:2-H7:7       : 46.4 [theta - 1.91218]^2
2:    C1:1-C2:2-H8:8       : 46.4 [theta - 1.91218]^2
3:    C2:2-C1:1-H4:4       : 46.3 [theta - 1.91637]^2
4:    C2:2-C1:1-H5:5       : 46.3 [theta - 1.91637]^2
...
8:    O3:3-C2:2-H8:8       : 50.9 [theta - 1.9244]^2
9:    H4:4-C1:1-H5:5       : 39.4 [theta - 1.87763]^2
10:    H4:4-C1:1-H6:6       : 39.4 [theta - 1.87763]^2
11:    H5:5-C1:1-H6:6       : 39.4 [theta - 1.87763]^2
12:    H7:7-C2:2-H8:8       : 39.2 [theta - 1.89298]^2
)

-----------------
Editing Dihedrals
-----------------

Dihedrals are a bit more complex to edit, since a single dihedral can contain multiple
terms, or *parts". For example, let's look at the dihedrals in our ethanol molecule:

>>> dihedrals = ethanol._sire_object.property("dihedral").potentials()
>>> for i, dihedral in enumerate(dihedrals):
...     print(i, dihedral)
0 FourAtomFunction( {CGIdx(0),Index(0)} <- {CGIdx(0),Index(1)} - {CGIdx(0),Index(2)} -> {CGIdx(0),Index(8)} : 0.25 cos(phi) + 0.16 cos(3 phi) + 0.41 )
1 FourAtomFunction( {CGIdx(0),Index(7)} <- {CGIdx(0),Index(1)} - {CGIdx(0),Index(2)} -> {CGIdx(0),Index(8)} : 0.166667 cos(3 phi) + 0.166667 )
2 FourAtomFunction( {CGIdx(0),Index(2)} <- {CGIdx(0),Index(1)} - {CGIdx(0),Index(0)} -> {CGIdx(0),Index(3)} : 0.25 cos(phi) + 0.25 )
3 FourAtomFunction( {CGIdx(0),Index(6)} <- {CGIdx(0),Index(1)} - {CGIdx(0),Index(2)} -> {CGIdx(0),Index(8)} : 0.166667 cos(3 phi) + 0.166667 )
4 FourAtomFunction( {CGIdx(0),Index(3)} <- {CGIdx(0),Index(0)} - {CGIdx(0),Index(1)} -> {CGIdx(0),Index(7)} : 0.155556 cos(3 phi) + 0.155556 )
5 FourAtomFunction( {CGIdx(0),Index(3)} <- {CGIdx(0),Index(0)} - {CGIdx(0),Index(1)} -> {CGIdx(0),Index(6)} : 0.155556 cos(3 phi) + 0.155556 )
6 FourAtomFunction( {CGIdx(0),Index(2)} <- {CGIdx(0),Index(1)} - {CGIdx(0),Index(0)} -> {CGIdx(0),Index(5)} : 0.25 cos(phi) + 0.25 )
7 FourAtomFunction( {CGIdx(0),Index(2)} <- {CGIdx(0),Index(1)} - {CGIdx(0),Index(0)} -> {CGIdx(0),Index(4)} : 0.25 cos(phi) + 0.25 )
8 FourAtomFunction( {CGIdx(0),Index(5)} <- {CGIdx(0),Index(0)} - {CGIdx(0),Index(1)} -> {CGIdx(0),Index(7)} : 0.155556 cos(3 phi) + 0.155556 )
9 FourAtomFunction( {CGIdx(0),Index(4)} <- {CGIdx(0),Index(0)} - {CGIdx(0),Index(1)} -> {CGIdx(0),Index(7)} : 0.155556 cos(3 phi) + 0.155556 )
10 FourAtomFunction( {CGIdx(0),Index(5)} <- {CGIdx(0),Index(0)} - {CGIdx(0),Index(1)} -> {CGIdx(0),Index(6)} : 0.155556 cos(3 phi) + 0.155556 )
11 FourAtomFunction( {CGIdx(0),Index(4)} <- {CGIdx(0),Index(0)} - {CGIdx(0),Index(1)} -> {CGIdx(0),Index(6)} : 0.155556 cos(3 phi) + 0.155556 )

Let's consider the first ``FourAtomFunction``, which has multiple terms, and convert it to
an ``AmberDihedral``object.

... Note::
   The containers used for bonded functions don't preserver order, so the
   dihedrals shown above may not be in the same order as you see when you run
   the code.

>>> from sire.legacy.MM import AmberDihedral
>>> from sire.legacy.CAS import Symbol
>>> amber_dihedral = AmberDihedral(dihedrals[0].function(), Symbol("phi"))

We can now print the terms in the dihedral:

>>> print(amber_dihedral.terms())
[AmberDihPart( k = 0.25, periodicity = 1, phase = 0 ), AmberDihPart( k = 0.16, periodicity = 3, phase = 0 )]

It's not currently possible to use ``AmberDihPart`` objects directly as a means
of building an ``AmberDihedral``. This is because this part of the legacy Sire
API was never intended to be used directly from Python, rather ``AmberDihedral``
objects would be created directly from expressions that are parsed from AMBER
topology files in the C++ API. However, it is easy enough to create multi-term
objects by writing custom expressions. The ``AmberDihedral`` code is written to
be robust against different AMBER-style dihedral representations from common
format. For example:

A regular AMBER-style dihedral series where all terms have positive cosine factors:

>>> from sire.legacy.CAS import Cos, Expression, Symbol
>>> Phi = Symbol("phi")
>>> f = Expression(0.3 * (1 + Cos(Phi)) + 0.8 * (1 + Cos(4 * Phi)))
>>> d = AmberDihedral(f, Phi)
>>> print("AMBER:", d)
AMBER: AmberDihedral( k[0] = 0.3, periodicity[0] = 1, phase[0] = 0, k[1] = 0.8, periodicity[1] = 4, phase[1] = 0 )
>>> assert d.toExpression(Phi) == f

An AMBER-style dihedral containing positive and negative cosine factors, which
can appear in the CHARMM force field:

>>> f = Expression(0.3 * (1 + Cos(Phi)) - 0.8 * (1 - Cos(4 * Phi)))
>>> d = AmberDihedral(f, Phi)
>>> print("CHARMM:", d)
CHARMM: AmberDihedral( k[0] = 0.3, periodicity[0] = 1, phase[0] = 0, k[1] = -0.8, periodicity[1] = 4, phase[1] = 0 )
>>> assert d.toExpression(Phi) == f

An AMBER-style dihedral containing positive and negative cosine factors, with
the negative of the form ``k [1 - Cos(Phi)]`` rather than ``-k [1 + Cos(Phi)]``.
These can appear in the GROMACS force field:

>>> f = Expression(0.3 * (1 + Cos(Phi)) + 0.8 * (1 - Cos(4 * Phi)))
>>> d = AmberDihedral(f, Phi)
>>> print("GROMACS:", d)
GROMACS: AmberDihedral( k[0] = 0.3, periodicity[0] = 1, phase[0] = 0, k[1] = 0.8, periodicity[1] = 4, phase[1] = -3.14159 )
>>> from math import isclose
>>> from sire.legacy.CAS import SymbolValue, Values
>>> val = Values(SymbolValue(Phi.ID(), 2.0))
>>> assert isclose(f.evaluate(val), d.toExpression(Phi).evaluate(val))

Finally, a three-term expression that mixes all formats:

>>> # Try a three-term expression that mixes all formats.
>>> f = Expression(
...     0.3 * (1 + Cos(Phi))
...     - 1.2 * (1 + Cos(3 * Phi))
...     + 0.8 * (1 - Cos(4 * Phi))
... )
>>> d = AmberDihedral(f, Phi)
>>> assert isclose(f.evaluate(val), d.toExpression(Phi).evaluate(val))
