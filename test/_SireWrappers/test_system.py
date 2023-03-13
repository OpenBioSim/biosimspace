import BioSimSpace as BSS

import math
import pytest

# Store the tutorial URL.
url = BSS.tutorialUrl()


@pytest.fixture(scope="session")
def system():
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules(["test/input/ala.top", "test/input/ala.crd"])


# Parameterise the function with a set of molecule indices.
@pytest.mark.parametrize("index", [0, -1])
def test_molecule_equivalence(system, index):
    # Make sure that we get the same molecule, however we extract it from
    # the system.
    assert system[index] == system.getMolecules()[index] == system.getMolecule(index)


def test_iterators(system):
    # Iterate over all molecules in the system, either directly,
    # or after calling getMolecules() and make sure they are equivalent.

    # Get the molecules from the system.
    molecules = system.getMolecules()

    # First iterate over the system object directly.
    for idx, mol in enumerate(system):
        assert mol == molecules[idx] == system.getMolecule(idx)

    # Now iterate over the molecules instead.
    for idx, mol in enumerate(molecules):
        assert mol == system[idx] == system.getMolecule(idx)

    # Now search for all molecules in the system.
    search_result = system.search("all")

    # Iterate over the search result and make sure we match all molecules.
    for idx, mol in enumerate(search_result):
        # Sire automatically converts objects to their smallest types
        # (e.g. residue, atom etc.)
        if hasattr(mol, "toMolecule"):
            mol = mol.toMolecule()

        assert mol == system[idx] == molecules[idx] == system.getMolecule(idx)


def test_atom_reindexing(system):
    # Search for all oxygen atoms in water molecules water molecules within
    # the system.
    results = system.search("waters and element oxygen")

    # There are 22 atoms in the alanine-dipeptide, then three in each water
    # molecule. The oxygen atoms come first in each water molecule, so
    # absolute indices should start at 22 and increment by 3.

    # Starting index.
    index = 22

    for atom in results:
        # Ensure the absolute index matches.
        assert system.getIndex(atom) == index

        # Increment the index. (3-point water.)
        index += 3


def test_residue_reindexing(system):
    # Search for all waters by residue name.
    results = system.search("resname WAT").residues()

    # There are 3 residues in the alanine-dipeptide, then one in each water
    # molecule. This means that residue indexing should start at 3 and
    # increment by 1.

    # Starting index.
    index = 3

    for residue in results:
        # Ensure the absolute index matches.
        assert system.getIndex(residue) == index

        index += 1


def test_molecule_reindexing(system):
    # Search for all waters by residue name.
    results = system.search("resname WAT").molecules()

    # There are 631 molecules in the system: an alanine-dipeptide, followed by
    # 630 water molecules. This means that molecule indexing should start at 1
    # and increment by 1.

    # Starting index.
    index = 1

    # Note that the waters will be returned as residue objects, since this
    # is the minimal representation, i.e. a water contains a single residue.
    # As such, we convert each result to a molecule.
    for residue in results:
        # Ensure the absolute index matches.
        assert system.getIndex(residue) == index

        index += 1


def test_contains(system):
    # Extract the first molecule.
    m = system[0]

    # Make sure the molecule is in the system.
    assert m in system

    # Make sure a copy of the molecule isn't in the system.
    assert m.copy() not in system

    # Extract the first residue of the molecule.
    r = m.getResidues()[0]

    # Extract the first atom of the residue.
    a = r[0]

    # Make sure the residue and atom are in the system.
    assert r in system
    assert a in system

    # Make sure the residue and atom are in the molecule.
    assert r in m
    assert a in m

    # Make sure that the atom is in the residue.
    assert a in r

    # Get an atom from a different residue.
    a = system[1].getResidues()[0][0]

    # Make sure the atom isn't in the residue.
    assert a not in r


def test_contains_combine(system):
    # Extract the first two molecules.
    m0 = system[0]
    m1 = system[1]

    # Combine the two molecules to make two systems, adding the
    # molecules in opposite order.
    s0 = (m0 + m1).toSystem()
    s1 = (m1 + m0).toSystem()

    # Make sure both molecules are in both systems.
    assert m0 in s0
    assert m0 in s1
    assert m1 in s0
    assert m1 in s1


def test_get_atom(system):
    # Make sure atom extraction works using absolute or relative indices,
    # i.e. absolute within the system, or relative to a molecule.
    assert system.getAtom(0) == system[0].getAtoms()[0]
    assert system.getAtom(22) == system[1].getAtoms()[0]
    assert system.getAtom(1883) == system[-10].getAtoms()[1]

    # Remove some molecules from the system.
    system2 = system.copy()
    system2.removeMolecules([system[0], system[2], system[-2]])

    # Check that the assertions still hold on the modified system,
    # making sure we map to the molecules in their new positions.
    assert system2.getAtom(0) == system2[0].getAtoms()[0]
    assert system2.getAtom(22) == system2[7].getAtoms()[1]
    assert system2.getAtom(1883) == system2[-1].getAtoms()[-1]


def test_get_residue(system):
    # Make sure residue extraction works using absolute or relative indices,
    # i.e. absolute within the system, or relative to a molecule.
    assert system.getResidue(0) == system[0].getResidues()[0]
    assert system.getResidue(3) == system[1].getResidues()[0]
    assert system.getResidue(632) == system[-1].getResidues()[-1]


def test_molecule_reordering(system):
    # Make sure that molecules are added to a container in the correct order.

    # Create a new system via a Molecules container where the last
    # molecule in the original system is placed in the first position
    # in the new one.
    new_system = (system[-1] + system[:-1]).toSystem()

    # Make sure the molecules match.
    assert new_system[0] == system[-1]


# Parameterise the function with a set of molecule indices.
@pytest.mark.parametrize(
    "restraint, expected",
    [
        ("backbone", [4, 5, 6, 8, 14, 15, 16]),
        ("heavy", [1, 4, 5, 6, 8, 10, 14, 15, 16, 18]),
        (
            "all",
            [
                0,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                13,
                14,
                15,
                16,
                17,
                18,
                19,
                20,
                21,
            ],
        ),
    ],
)
def test_restraint_atoms(system, restraint, expected):
    # Get the restraint atoms for the specicied restraint.
    atoms = system.getRestraintAtoms(restraint)

    # Make sure the indices are as expected.
    assert atoms == expected

    # Execute the restraint search at the molecular level.
    atoms = system.getRestraintAtoms(restraint, 0)

    # Make sure the indices are as expected.
    assert atoms == expected


def test_get_box(system):
    # Get the box dimensions and angles from the system.
    box, angles = system.getBox()

    # Store the expected box dimensions and angles.
    expected_box = [31.3979, 34.1000, 29.2730] * BSS.Units.Length.angstrom
    expected_angles = [90, 90, 90] * BSS.Units.Angle.degree

    # Check that the box dimensions match.
    for b0, b1 in zip(box, expected_box):
        assert math.isclose(b0.value(), b1.value(), rel_tol=1e-4)

    # Check that the angles match.
    for a0, a1 in zip(angles, expected_angles):
        assert math.isclose(a0.value(), a1.value(), rel_tol=1e-4)


def test_set_box(system):
    # Generate box dimensions and angles for a truncated octahedron.
    box, angles = BSS.Box.truncatedOctahedron(30 * BSS.Units.Length.angstrom)

    # Set the box dimensions in the system.
    system.setBox(box, angles)

    # Get the updated box dimensions and angles.
    box, angles = system.getBox()

    # Store the expected box dimensions and angles. These are different
    # to the initial values (angles at least) due to the application
    # of a lattice reduction.
    expected_box = [30, 30, 30] * BSS.Units.Length.angstrom
    expected_angles = [70.5288, 109.4712, 70.5288] * BSS.Units.Angle.degree

    # Check that the box dimensions match.
    for b0, b1 in zip(box, expected_box):
        assert math.isclose(b0.value(), b1.value(), rel_tol=1e-4)

    # Check that the angles match.
    for a0, a1 in zip(angles, expected_angles):
        assert math.isclose(a0.value(), a1.value(), rel_tol=1e-4)


def test_molecule_replace(system):
    # Make sure that molecule ordering is preserved when a molecule is
    # replaced by another with a different MolNum.

    # Extract the third molecule.
    mol0 = system[3]

    # Store the current molecule numbers.
    mol_nums0 = system._mol_nums

    # Update (replace) the third molecule with a renumbered
    # version.
    system.updateMolecule(3, mol0.copy())

    # Get the third molecule in the updated system.
    mol1 = system[3]

    # Store the updated molecule numbers.
    mol_nums1 = system._mol_nums

    # Make sure the molecules have different numbers.
    assert mol0.number() != mol1.number()

    # Make sure the molecules have the same number of atoms and residues.
    assert mol0.nAtoms() == mol1.nAtoms()
    assert mol0.nResidues() == mol1.nResidues()

    # Make sure that the third MolNum in the array is different.
    assert mol_nums0[3] != mol_nums1[3]

    # Make sure the rest match.
    for num0, num1 in zip(mol_nums0[0:2], mol_nums1[0:2]):
        assert num0 == num1
    for num0, num1 in zip(mol_nums0[4:], mol_nums1[4:]):
        assert num0 == num1


def test_isSame(system):
    # Make sure that the isSame method works correctly.

    # Make a copy of the system.
    other = system.copy()

    # Assert they are the same, making sure it is invariant to the order.
    assert system.isSame(other)
    assert other.isSame(system)

    # Translate the other system.
    other.translate(3 * [BSS.Units.Length.angstrom])

    # Assert that they are different.
    assert not system.isSame(other)
    assert not other.isSame(system)

    # Assert that they are the same, apart from their coordinates.
    assert system.isSame(other, excluded_properties=["coordinates"])
    assert other.isSame(system, excluded_properties=["coordinates"])

    # Now delete a property.
    other._sire_object.removeProperty("space")

    # Assert that they are different.
    assert not system.isSame(other, excluded_properties=["coordinates"])
    assert not other.isSame(system, excluded_properties=["coordinates"])

    # Assert that they are the same, apart from their coordinates and space.
    assert system.isSame(other, excluded_properties=["coordinates", "space"])
    assert other.isSame(system, excluded_properties=["coordinates", "space"])


def test_velocity_removal():
    # Make sure that velocities are removed when molecules are combined
    # and not all molecules have a "velocity" property.

    # Load a molecule with and without velocities.
    mol = BSS.IO.readMolecules(BSS.IO.expand(url, "methane.gro", ".bz2"))
    mol_vel = BSS.IO.readMolecules(BSS.IO.expand(url, "methane_vel.gro", ".bz2"))

    # Add together to create a new system.
    new_system = mol + mol_vel

    # Check that no molecules have a velocity property.
    assert len(new_system.search("not mol with property velocity").molecules()) == 2
