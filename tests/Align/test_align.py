import pytest
import sys

from sire.legacy.MM import InternalFF, IntraCLJFF, IntraFF
from sire.legacy.Mol import AtomIdx, Element, PartialMolecule

import BioSimSpace as BSS

# Store the tutorial URL.
url = BSS.tutorialUrl()


@pytest.fixture(scope="session")
def system0():
    return BSS.IO.readMolecules(
        [f"{url}/ligand01.prm7.bz2", f"{url}/ligand01.rst7.bz2"]
    )


@pytest.fixture(scope="session")
def system1():
    return BSS.IO.readMolecules(
        [f"{url}/ligand02.prm7.bz2", f"{url}/ligand02.rst7.bz2"]
    )


def test_flex_align(system0, system1):
    # This tests that the flex align functionality runs. We can't test
    # for consistent output, since we have occasionally observed different
    # mappings across platforms.

    # Extract the molecules.
    m0 = system0.getMolecules()[0]
    m1 = system1.getMolecules()[0]

    # Get the best mapping between the molecules that contains the prematch.
    mapping = BSS.Align.matchAtoms(
        m0, m1, timeout=BSS.Units.Time.second, scoring_function="rmsd_flex_align"
    )


# Parameterise the function with a set of valid atom pre-matches.
@pytest.mark.skipif(
    sys.platform == "win32", reason="Sire MCS currently not supported on Windows"
)
@pytest.mark.parametrize("prematch", [{3: 1}, {5: 9}, {4: 5}, {1: 0}])
def test_prematch(system0, system1, prematch):
    # Extract the molecules.
    m0 = system0.getMolecules()[0]
    m1 = system1.getMolecules()[0]

    # Get the best mapping between the molecules that contains the prematch.
    mapping = BSS.Align.matchAtoms(
        m0, m1, timeout=BSS.Units.Time.second, prematch=prematch
    )

    # Check that the prematch key:value pair is in the mapping.
    for key, value in prematch.items():
        assert mapping[key] == value


# Parameterise the function with a set of invalid atom pre-matches.
@pytest.mark.parametrize("prematch", [{-1: 1}, {50: 9}, {4: 48}, {1: -1}])
def test_invalid_prematch(system0, system1, prematch):
    # Extract the molecules.
    m0 = system0.getMolecules()[0]
    m1 = system1.getMolecules()[0]

    # Assert that the invalid prematch raises a ValueError.
    with pytest.raises(ValueError):
        mapping = BSS.Align.matchAtoms(
            m0, m1, timeout=BSS.Units.Time.second, prematch=prematch
        )


def test_merge():
    # Load the ligands.
    s0 = BSS.IO.readMolecules([f"{url}/ligand31.prm7.bz2", f"{url}/ligand31.rst7.bz2"])
    s1 = BSS.IO.readMolecules([f"{url}/ligand38.prm7.bz2", f"{url}/ligand38.rst7.bz2"])

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Get the best mapping between the molecules.
    mapping = BSS.Align.matchAtoms(m0, m1, timeout=BSS.Units.Time.second)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Create the merged molecule.
    m2 = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True)

    # Store the number of atoms in m0.
    n0 = m0._sire_object.nAtoms()

    # Test that the intramolecular energies area the same.

    # IntraCLJFF:
    #  Old interface. Uses the "intrascale" matrix. Validate that this
    #  is consistent.
    # IntraFF:
    #  New interface. Uses atom "connectivity". Validate that the bonding
    #  is consistent.

    intraclj0 = IntraCLJFF("intraclj")
    intraclj0.add(m0._sire_object)

    intraff0 = IntraFF("intraclj")
    intraff0.add(m0._sire_object)

    intraclj1 = IntraCLJFF("intraclj")
    intraclj1.add(m1._sire_object)

    intraff1 = IntraFF("intraclj")
    intraff1.add(m1._sire_object)

    intraclj2 = IntraCLJFF("intraclj")
    intraff2 = IntraFF("intraclj")

    # Create maps between property names: { "prop" : "prop0" }, { "prop" : "prop1" }
    pmap0 = {}
    pmap1 = {}
    for prop in m2._sire_object.propertyKeys():
        if prop[-1] == "0":
            pmap0[prop[:-1]] = prop
        elif prop[-1] == "1":
            pmap1[prop[:-1]] = prop

    intraclj2.add(m2._sire_object, pmap0)
    intraff2.add(m2._sire_object, pmap0)

    assert intraclj0.energy().value() == pytest.approx(intraclj2.energy().value())
    assert intraff0.energy().value() == pytest.approx(intraff2.energy().value())

    intraclj2 = IntraCLJFF("intraclj")
    intraff2 = IntraFF("intraclj")

    intraclj2.add(m2._sire_object, pmap1)
    intraff2.add(m2._sire_object, pmap1)

    assert intraclj1.energy().value() == pytest.approx(intraclj2.energy().value())
    assert intraff1.energy().value() == pytest.approx(intraff2.energy().value())

    # Test that the internal energies are consistent. This will validate that
    # bond, angle, dihedral, and improper energies are correct.

    internalff0 = InternalFF("internal")
    internalff0.setStrict(True)
    internalff0.add(m0._sire_object)

    internalff1 = InternalFF("internal")
    internalff1.setStrict(True)
    internalff1.add(m1._sire_object)

    # First extract a partial molecule using the atoms from molecule0 in
    # the merged molecule.
    selection = m2._sire_object.selection()
    selection.deselectAll()
    for atom in m0._sire_object.atoms():
        selection.select(atom.index())
    partial_mol = PartialMolecule(m2._sire_object, selection)

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(partial_mol, pmap0)

    assert internalff0.energy().value() == pytest.approx(internalff2.energy().value())

    # Extract the original molecule for the lambda=0 end state.
    amber_mol, _ = m2._extractMolecule()

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(amber_mol._sire_object)

    assert internalff0.energy().value() == pytest.approx(internalff2.energy().value())

    # Now extract a partial molecule using the atoms from molecule1 in
    # the merged molecule.
    selection = m2._sire_object.selection()
    selection.deselectAll()
    for idx in mapping.keys():
        selection.select(AtomIdx(idx))
    for idx in range(n0, m2._sire_object.nAtoms()):
        selection.select(AtomIdx(idx))
    partial_mol = PartialMolecule(m2._sire_object, selection)

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(partial_mol, pmap1)

    assert internalff1.energy().value() == pytest.approx(internalff2.energy().value())

    # Extract the original molecule for the lambda=1 end state.
    amber_mol, _ = m2._extractMolecule(is_lambda1=True)

    internalff2 = InternalFF("internal")
    internalff2.setStrict(True)
    internalff2.add(amber_mol._sire_object)

    assert internalff1.energy().value() == pytest.approx(internalff2.energy().value())


@pytest.mark.xfail(
    reason="Mapping generated with latest RDKit which requires sanitization no longer triggers the exception"
)
def test_ring_breaking_three_membered():
    # Load the ligands.
    s0 = BSS.IO.readMolecules([f"{url}/CAT-13a.prm7.bz2", f"{url}/CAT-13a.rst7.bz2"])
    s1 = BSS.IO.readMolecules([f"{url}/CAT-17g.prm7.bz2", f"{url}/CAT-17g.rst7.bz2"])

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Generate the mapping.
    mapping = BSS.Align.matchAtoms(m0, m1)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Try to merge the molecule without allowing ring breaking.
    with pytest.raises(BSS._Exceptions.IncompatibleError):
        m2 = BSS.Align.merge(m0, m1, mapping)

    # Now check that we can merge if we allow ring breaking.
    m2 = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True)


@pytest.mark.xfail(
    reason="Mapping generated with latest RDKit which requires sanitization no longer triggers the exception"
)
def test_ring_breaking_five_membered():
    # Load the ligands.
    s0 = BSS.IO.readMolecules([f"{url}/ligand31.prm7.bz2", f"{url}/ligand31.rst7.bz2"])
    s1 = BSS.IO.readMolecules([f"{url}/ligand04.prm7.bz2", f"{url}/ligand04.rst7.bz2"])

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Load the pre-defined mapping.
    mapping = BSS.Align.matchAtoms(m0, m1)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Try to merge the molecule without allowing ring breaking.
    with pytest.raises(BSS._Exceptions.IncompatibleError):
        m2 = BSS.Align.merge(m0, m1, mapping)

    # Now check that we can merge if we allow ring breaking.
    m2 = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True)


@pytest.mark.xfail(
    reason="Mapping generated with latest RDKit which requires sanitization no longer triggers the exception"
)
def test_ring_breaking_six_membered():
    # Load the ligands.
    s0 = BSS.IO.readMolecules([f"{url}/ligand31.prm7.bz2", f"{url}/ligand31.rst7.bz2"])
    s1 = BSS.IO.readMolecules([f"{url}/ligand38.prm7.bz2", f"{url}/ligand38.rst7.bz2"])

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Load the pre-defined mapping.
    mapping = BSS.Align.matchAtoms(m0, m1)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Try to merge the molecule without allowing ring breaking.
    with pytest.raises(BSS._Exceptions.IncompatibleError):
        m2 = BSS.Align.merge(m0, m1, mapping)

    # Now check that we can merge if we allow ring breaking.
    m2 = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True)


@pytest.mark.parametrize(
    "ligands",
    [
        pytest.param(
            ["CAT-13c", "CAT-17i"],
            marks=pytest.mark.xfail(
                reason="Mapping generated with latest RDKit which requires sanitization no longer triggers the exception"
            ),
        ),
        pytest.param(
            ["CAT-13e", "CAT-17g"],
            marks=pytest.mark.xfail(
                reason="Mapping generated with latest RDKit which requires sanitization no longer triggers the exception"
            ),
        ),
    ],
)
def test_ring_size_change(ligands):
    # Load the ligands.
    s0 = BSS.IO.readMolecules(
        [f"{url}/{ligands[0]}.prm7.bz2", f"{url}/{ligands[0]}.rst7.bz2"]
    )
    s1 = BSS.IO.readMolecules(
        [f"{url}/{ligands[1]}.prm7.bz2", f"{url}/{ligands[1]}.rst7.bz2"]
    )

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Generate the mapping.
    mapping = BSS.Align.matchAtoms(m0, m1)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Try to merge the molecule without allowing ring breaking.
    with pytest.raises(BSS._Exceptions.IncompatibleError):
        m2 = BSS.Align.merge(m0, m1, mapping)

    # Now check that we can merge if we allow ring breaking.
    m2 = BSS.Align.merge(
        m0, m1, mapping, allow_ring_breaking=True, allow_ring_size_change=True
    )


# Parameterise the function with a valid mapping.
@pytest.mark.parametrize(
    "ligands, mapping",
    [
        (
            ("grow1", "grow2"),
            {
                2: 21,
                4: 23,
                6: 25,
                8: 27,
                10: 18,
                1: 19,
                0: 20,
                11: 16,
                12: 17,
                13: 14,
                15: 13,
                18: 11,
                20: 9,
                22: 8,
                23: 5,
                16: 6,
                24: 3,
                26: 1,
                27: 0,
                9: 28,
                5: 24,
                3: 22,
                7: 26,
                14: 15,
                19: 12,
                21: 10,
                17: 7,
                25: 4,
            },
        ),
        (
            ("grow3", "grow4"),
            {
                1: 6,
                2: 7,
                3: 8,
                4: 9,
                5: 10,
                6: 11,
                14: 21,
                13: 20,
                12: 19,
                11: 18,
                10: 17,
            },
        ),
    ],
)
def test_grow_whole_ring(ligands, mapping):
    # Load the ligands.
    s0 = BSS.IO.readMolecules(
        [f"{url}/{ligands[0]}.prm7.bz2", f"{url}/{ligands[0]}.rst7.bz2"]
    )
    s1 = BSS.IO.readMolecules(
        [f"{url}/{ligands[1]}.prm7.bz2", f"{url}/{ligands[1]}.rst7.bz2"]
    )

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Check that we can merge without allowing ring breaking.
    m2 = BSS.Align.merge(m0, m1, mapping)


def test_hydrogen_mass_repartitioning():
    # Load the ligands.
    s0 = BSS.IO.readMolecules([f"{url}/ligand31.prm7.bz2", f"{url}/ligand31.rst7.bz2"])
    s1 = BSS.IO.readMolecules([f"{url}/ligand38.prm7.bz2", f"{url}/ligand38.rst7.bz2"])

    # Extract the molecules.
    m0 = s0.getMolecules()[0]
    m1 = s1.getMolecules()[0]

    # Get the best mapping between the molecules.
    mapping = BSS.Align.matchAtoms(m0, m1, timeout=BSS.Units.Time.second)

    # Align m0 to m1 based on the mapping.
    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    # Create the merged molecule.
    merged = BSS.Align.merge(m0, m1, mapping, allow_ring_breaking=True)

    # Create a dummy element.
    dummy = Element("Xx")

    # Get the elements in either end state.
    elements0 = merged._sire_object.property("element0").toVector()
    elements1 = merged._sire_object.property("element1").toVector()

    # Work out the initial mass of the system.
    initial_mass0 = 0
    for idx, mass in enumerate(merged._sire_object.property("mass0").toVector()):
        if elements0[idx] != dummy:
            initial_mass0 += mass.value()
    initial_mass1 = 0
    for idx, mass in enumerate(merged._sire_object.property("mass1").toVector()):
        if elements1[idx] != dummy:
            initial_mass1 += mass.value()

    # Repartition the hydrogen mass.
    merged.repartitionHydrogenMass()

    # Lists to store the mass of dummy atoms in the two end states.
    dummy_masses0 = []
    dummy_masses1 = []

    # Extract the modified end state masses.
    masses0 = merged._sire_object.property("mass0").toVector()
    masses1 = merged._sire_object.property("mass1").toVector()

    # Work out the final mass of the system.
    final_mass0 = 0
    for idx, mass in enumerate(masses0):
        if elements0[idx] != dummy:
            final_mass0 += mass.value()
        else:
            dummy_masses0.append((idx, mass))
    final_mass1 = 0
    for idx, mass in enumerate(masses1):
        if elements1[idx] != dummy:
            final_mass1 += mass.value()
        else:
            dummy_masses1.append((idx, mass))

    # Assert the the masses are approximately the same.
    assert final_mass0 == pytest.approx(initial_mass0)
    assert final_mass1 == pytest.approx(initial_mass1)

    # Assert that the dummy atom masses are the same in both end states.
    for idx, mass0 in dummy_masses0:
        assert mass0 == masses1[idx]
    for idx, mass1 in dummy_masses1:
        assert mass1 == masses0[idx]


@pytest.fixture
def proteins(request):
    p0 = BSS.IO.readMolecules(["ala_ser_size_5_capped.pdb"])[0]
    p1 = BSS.IO.readMolecules(["ala_wt_size_5_capped.pdb"])[0]
    return p0, p1


@pytest.fixture(
    params=[
        {
            0: 0,
            1: 1,
            2: 2,
            3: 3,
            4: 4,
            5: 5,
            6: 6,
            7: 7,
            8: 8,
            9: 9,
            10: 10,
            11: 11,
            12: 12,
            13: 13,
            14: 14,
            15: 15,
            16: 16,
            17: 17,
            18: 18,
            19: 19,
            20: 20,
            21: 23,
            22: 21,
            23: 22,
            25: 24,
            26: 25,
            27: 26,
            28: 27,
            29: 28,
            30: 29,
            31: 30,
            32: 31,
            33: 32,
            34: 33,
            35: 34,
            36: 35,
            37: 36,
            38: 37,
            39: 38,
            40: 39,
            41: 40,
            42: 41,
        }
    ]
)
def protein_mapping(request):
    return request.param


def test_roi_match(proteins, protein_mapping):
    p0, p1 = proteins
    mapping = BSS.Align.roiMatch(p0, p1, roi=[2])
    assert mapping == protein_mapping


def test_roi_align(proteins, protein_mapping):
    # p0 has been translated by 10 A in each direction.
    p0, p1 = proteins

    aligned_p0 = BSS.Align.roiAlign(p0, p1, roi=[2])

    # Extract sire objects for the ROI and compare their coordinates
    aligned_roi = aligned_p0.extract(
        [a.index() for a in aligned_p0.getResidues()[2].getAtoms()]
    )
    aligned_roi_coords = aligned_roi._sire_object.coordinates()

    p1_roi = p1.extract([a.index() for a in p1.getResidues()[2].getAtoms()])
    p1_roi_coords = p1_roi._sire_object.coordinates()

    for i, coord in enumerate(aligned_roi_coords):
        # assume that the coordinates are the same if they are within 0.5 A
        assert coord.value() == pytest.approx(p1_roi_coords[i].value(), abs=0.5)


def test_roi_merge(proteins, protein_mapping):
    p0, p1 = proteins

    p0 = BSS.Parameters.ff14SB(p0).getMolecule()
    p1 = BSS.Parameters.ff14SB(p1).getMolecule()

    aligned_p0 = BSS.Align.roiAlign(p0, p1, roi=[2])
    merged = BSS.Align.merge(aligned_p0, p1, protein_mapping)
    merged_system = merged.toSystem()
    assert merged_system.nPerturbableMolecules() == 1
