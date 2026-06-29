import sys

import pytest
import sire as sr
from sire.legacy.MM import InternalFF, IntraCLJFF, IntraFF
from sire.legacy.Mol import AtomIdx, Element, PartialMolecule

import BioSimSpace as BSS
from tests.conftest import has_amber, has_antechamber, has_openff, has_tleap

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


@pytest.mark.parametrize(
    "s0_files,s1_files",
    [
        (
            [f"{url}/ligand31.prm7.bz2", f"{url}/ligand31.rst7.bz2"],
            [f"{url}/ligand38.prm7.bz2", f"{url}/ligand38.rst7.bz2"],
        ),
        pytest.param(
            [f"{url}/glycam_M5.gro.bz2", f"{url}/glycam_M5.top.bz2"],
            [f"{url}/glycam_M5G0.gro.bz2", f"{url}/glycam_M5G0.top.bz2"],
            marks=pytest.mark.slow,
        ),
    ],
)
def test_merge(s0_files, s1_files):
    # Load the molecules.
    s0 = BSS.IO.readMolecules(s0_files)
    s1 = BSS.IO.readMolecules(s1_files)

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
    n0 = m0._sire_object.num_atoms()

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
    for prop in m2._sire_object.property_keys():
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
    internalff0.set_strict(True)
    internalff0.add(m0._sire_object)

    internalff1 = InternalFF("internal")
    internalff1.set_strict(True)
    internalff1.add(m1._sire_object)

    # First extract a partial molecule using the atoms from molecule0 in
    # the merged molecule.
    selection = m2._sire_object.selection()
    selection.deselect_all()
    for atom in m0._sire_object.atoms():
        selection.select(atom.index())
    partial_mol = PartialMolecule(m2._sire_object, selection)

    internalff2 = InternalFF("internal")
    internalff2.set_strict(True)
    internalff2.add(partial_mol, pmap0)

    assert internalff0.energy().value() == pytest.approx(internalff2.energy().value())

    # Extract the original molecule for the lambda=0 end state.
    amber_mol, _ = m2._extractMolecule()

    internalff2 = InternalFF("internal")
    internalff2.set_strict(True)
    internalff2.add(amber_mol._sire_object)

    assert internalff0.energy().value() == pytest.approx(internalff2.energy().value())

    # Now extract a partial molecule using the atoms from molecule1 in
    # the merged molecule.
    selection = m2._sire_object.selection()
    selection.deselect_all()
    for idx in mapping.keys():
        selection.select(AtomIdx(idx))
    for idx in range(n0, m2._sire_object.num_atoms()):
        selection.select(AtomIdx(idx))
    partial_mol = PartialMolecule(m2._sire_object, selection)

    internalff2 = InternalFF("internal")
    internalff2.set_strict(True)
    internalff2.add(partial_mol, pmap1)

    assert internalff1.energy().value() == pytest.approx(internalff2.energy().value())

    # Extract the original molecule for the lambda=1 end state.
    amber_mol, _ = m2._extractMolecule(is_lambda1=True)

    internalff2 = InternalFF("internal")
    internalff2.set_strict(True)
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
    elements0 = merged._sire_object.property("element0").to_vector()
    elements1 = merged._sire_object.property("element1").to_vector()

    # Work out the initial mass of the system.
    initial_mass0 = 0
    for idx, mass in enumerate(merged._sire_object.property("mass0").to_vector()):
        if elements0[idx] != dummy:
            initial_mass0 += mass.value()
    initial_mass1 = 0
    for idx, mass in enumerate(merged._sire_object.property("mass1").to_vector()):
        if elements1[idx] != dummy:
            initial_mass1 += mass.value()

    # Repartition the hydrogen mass.
    merged.repartitionHydrogenMass()

    # Lists to store the mass of dummy atoms in the two end states.
    dummy_masses0 = []
    dummy_masses1 = []

    # Extract the modified end state masses.
    masses0 = merged._sire_object.property("mass0").to_vector()
    masses1 = merged._sire_object.property("mass1").to_vector()

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


@pytest.fixture(
    params=[
        (
            "single_mutant",
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
            },
            [2],
        ),
        (
            "double_neighbour_mutant",
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
                21: 24,
                22: 25,
                23: 26,
                24: 27,
                25: 28,
                26: 29,
                27: 30,
                28: 34,
                29: 35,
                30: 36,
                31: 37,
                32: 38,
                33: 39,
                34: 40,
                35: 41,
                36: 42,
                37: 43,
                38: 44,
                39: 45,
                40: 46,
                41: 47,
                42: 48,
                43: 49,
                44: 50,
                45: 51,
                46: 52,
                47: 53,
                48: 54,
                49: 55,
                50: 56,
                51: 57,
                52: 58,
                53: 59,
                54: 60,
                55: 61,
            },
            [2, 3],
        ),
    ],
    ids=["single_mutant", "double_neighbour_mutant"],
)
def protein_inputs(request):
    return request.param


def test_roi_match(protein_inputs):
    proteins, protein_mapping, roi = protein_inputs
    p0 = BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), f"{proteins}_mut_peptide.pdb")
    )[0]
    p1 = BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), f"{proteins}_wt_peptide.pdb")
    )[0]
    mapping = BSS.Align.matchAtoms(p0, p1, roi=roi)
    assert mapping == protein_mapping


def test_roi_align(protein_inputs):
    # p0 has been translated by 10 A in each direction.
    proteins, protein_mapping, roi = protein_inputs
    p0 = BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), f"{proteins}_mut_peptide.pdb")
    )[0]
    p1 = BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), f"{proteins}_wt_peptide.pdb")
    )[0]

    aligned_p0 = BSS.Align.rmsdAlign(p0, p1, roi=roi)
    for res in roi:
        # Extract sire objects for the ROI and compare their coordinates
        aligned_roi = aligned_p0.extract(
            [a.index() for a in aligned_p0.getResidues()[res].getAtoms()]
        )
        aligned_roi_coords = aligned_roi._sire_object.coordinates()

        p1_roi = p1.extract([a.index() for a in p1.getResidues()[res].getAtoms()])
        p1_roi_coords = p1_roi._sire_object.coordinates()

        for i, coord in enumerate(aligned_roi_coords):
            # assume that the test passes if the coordinates are within 0.5 A
            assert coord.value() == pytest.approx(p1_roi_coords[i].value(), abs=0.5)


def test_roi_flex_align(protein_inputs):
    # p0 has been translated by 10 A in each direction.
    proteins, protein_mapping, roi = protein_inputs
    p0 = BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), f"{proteins}_mut_peptide.pdb")
    )[0]
    p1 = BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), f"{proteins}_wt_peptide.pdb")
    )[0]

    aligned_p0 = BSS.Align.flexAlign(p0, p1, roi=roi)
    for res in roi:
        # Extract sire objects for the ROI and compare their coordinates
        aligned_roi = aligned_p0.extract(
            [a.index() for a in aligned_p0.getResidues()[res].getAtoms()]
        )
        aligned_roi_coords = aligned_roi._sire_object.coordinates()

        p1_roi = p1.extract([a.index() for a in p1.getResidues()[res].getAtoms()])
        p1_roi_coords = p1_roi._sire_object.coordinates()

        for i, coord in enumerate(aligned_roi_coords):
            # assume that the test passes if the coordinates are within 0.5 A
            assert coord.value() == pytest.approx(p1_roi_coords[i].value(), abs=0.5)


def test_empty_custom_roi_mapping():
    # mut contains a proline mutation at position 15
    wt = BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), f"1choFH_apo_wt_flare_processed.pdb")
    )[0]
    mut = BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), f"1choFH_apo_mut_flare_processed.pdb")
    )[0]

    # use the custom_roi_map to specify that residue 15 in the WT protein should be
    # excluded from the ROI mapping, even though it is in the ROI list
    roi_res_idx = [a.index() for a in wt.getResidues()[15].getAtoms()]
    mapping = BSS.Align.matchAtoms(
        molecule0=wt, molecule1=mut, roi=[15], custom_roi_map={}
    )

    # check that the mapping does not contain any atoms of the region of interest of WT protein
    for atom_idx in roi_res_idx:
        assert atom_idx not in mapping.keys()

@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
def test_custom_roi_ring_break_merge():
    # wt contains a leucine at position 15
    # mut contains a proline at position 15
    wt = BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), f"1choFH_apo_wt_flare_processed.pdb")
    )[0]
    mut = BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), f"1choFH_apo_mut_flare_processed.pdb")
    )[0]

    wt = BSS.Parameters.ff14SB(wt, ensure_compatible=False).getMolecule()
    mut = BSS.Parameters.ff14SB(mut, ensure_compatible=False).getMolecule()

    # use the custom_roi_map to specify that residue 15 in the WT protein should be
    # excluded from the ROI mapping, even though it is in the ROI list
    mapping = BSS.Align.matchAtoms(
        molecule0=wt,
        molecule1=mut,
        roi=[15],
        custom_roi_map={
            204: 204,
            205: 205,
            203: 203,
            202: 202,
            211: 208,
            206: 206,
            213: 210,
            207: 207,
            214: 211,
            210: 213,
        },
    )

    aligned_wt = BSS.Align.rmsdAlign(molecule0=wt, mapping=mapping, molecule1=mut)
    merged_protein = BSS.Align.merge(
        aligned_wt, mut, mapping, allow_ring_breaking=True, roi=[15]
    )

    merged_protein_sire = merged_protein._sire_object
    pert = merged_protein_sire.perturbation()
    pert_omm = pert.to_openmm(map={"coordinates": "coordinates0"})

    changed_bonds_df = pert_omm.changed_bonds(to_pandas=True)
    n_bonds_created = (changed_bonds_df["k0"] == 0).sum()
    n_bonds_annihilated = (changed_bonds_df["k1"] == 0).sum()

    # assert that exactly one bond is being created, as mutating
    # from leucine to proline should create a new bond in the ring of proline
    assert n_bonds_created == 1
    assert n_bonds_annihilated == 0

@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
def test_custom_roi_map_invalid_outside_roi():
    wt = BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), f"1choFH_apo_wt_flare_processed.pdb")
    )[0]
    mut = BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), f"1choFH_apo_mut_flare_processed.pdb")
    )[0]

    wt = BSS.Parameters.ff14SB(wt, ensure_compatible=False).getMolecule()
    mut = BSS.Parameters.ff14SB(mut, ensure_compatible=False).getMolecule()

    # provide some invalid mapping that is outside of the ROI, which should raise an error
    with pytest.raises(ValueError):
        mapping = BSS.Align.matchAtoms(
            molecule0=wt,
            molecule1=mut,
            roi=[15],
        
            custom_roi_map={
                0: 0,
                1: 1,
                2: 2,
            },
        )


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER and to be installed.")
def test_roi_merge(protein_inputs):
    proteins, protein_mapping, roi = protein_inputs
    p0 = BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), f"{proteins}_mut_peptide.pdb")
    )[0]
    p1 = BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), f"{proteins}_wt_peptide.pdb")
    )[0]

    p0 = BSS.Parameters.ff14SB(p0).getMolecule()
    p1 = BSS.Parameters.ff14SB(p1).getMolecule()

    aligned_p0 = BSS.Align.rmsdAlign(p0, p1, roi=roi)
    merged = BSS.Align.merge(aligned_p0, p1, protein_mapping, roi=roi)
    merged_system = merged.toSystem()
    assert merged_system.nPerturbableMolecules() == 1


def test_ion_merge(system):
    from sire.legacy.IO import createSodiumIon

    # Extract a water molecule.
    water = system[-1]

    # Create a sodium ion using the water coordinates.
    ion = createSodiumIon(
        water.getAtoms()[0]._sire_object.property("coordinates"), "tip3p"
    )

    # Merge the water and ion.
    merged = BSS.Align.merge(water, BSS._SireWrappers.Molecule(ion))

    # Make sure the ion has the coordintes of the oxygen atom.
    coords0 = merged._sire_object.property("coordinates0").to_vector()[0]
    coords1 = merged._sire_object.property("coordinates1").to_vector()[0]
    water_coords = water._sire_object.property("coordinates").to_vector()[0]
    assert coords0 == coords1
    assert coords0 == water_coords


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
@pytest.mark.skipif(has_openff is False, reason="Requires OpenFF to be installed.")
@pytest.mark.parametrize(
    "ligands, mapping",
    [
        (
            ("Schindler_SYK_CHEMBL3265035", "Schindler_SYK_CHEMBL3265033"),
            {
                5: 8,
                6: 32,
                7: 31,
                8: 11,
                9: 12,
                10: 13,
                11: 14,
                12: 15,
                13: 16,
                14: 17,
                15: 18,
                16: 19,
                17: 20,
                18: 21,
                19: 22,
                20: 23,
                21: 24,
                22: 25,
                23: 26,
                24: 27,
                25: 28,
                26: 29,
                27: 30,
                28: 10,
                29: 9,
                38: 53,
                39: 52,
                40: 43,
                41: 44,
                42: 45,
                43: 46,
                44: 47,
                45: 48,
                46: 49,
                47: 50,
                48: 51,
                49: 42,
                50: 41,
                1: 1,
                2: 4,
                34: 36,
                33: 35,
                3: 3,
                4: 34,
                35: 1,
                0: 6,
                32: 37,
                31: 7,
                36: 33,
            },
        ),
        (
            ("Schindler_SYK_CHEMBL3265035", "Schindler_SYK_CHEMBL3265036"),
            {
                5: 7,
                6: 31,
                7: 30,
                8: 10,
                9: 11,
                10: 12,
                11: 13,
                12: 14,
                13: 15,
                14: 16,
                15: 17,
                16: 18,
                17: 19,
                18: 20,
                19: 21,
                20: 22,
                21: 23,
                22: 24,
                23: 25,
                24: 26,
                25: 27,
                26: 28,
                27: 29,
                28: 9,
                29: 8,
                38: 54,
                39: 53,
                40: 44,
                41: 45,
                42: 46,
                43: 47,
                44: 48,
                45: 49,
                46: 50,
                47: 51,
                48: 52,
                49: 43,
                50: 42,
                0: 4,
                1: 5,
                2: 6,
                3: 0,
                4: 1,
                30: 3,
                31: 39,
                32: 38,
                33: 40,
                34: 41,
                35: 32,
                36: 33,
                37: 34,
            },
        ),
    ],
)
@pytest.mark.skipif(has_openff is False, reason="Requires OpenFF to be installed.")
def test_ring_opening_and_size_change(ligands, mapping):
    # These perturbations involve ring formation (acyclic atoms in mol0 become
    # ring members in mol1) combined with ring size changes in the existing
    # fused bicyclic core. Check that the merge succeeds when both flags are
    # allowed.
    m0 = BSS.IO.readMolecules(f"{url}/{ligands[0]}.sdf.bz2")[0]
    m1 = BSS.IO.readMolecules(f"{url}/{ligands[1]}.sdf.bz2")[0]

    m0 = BSS.Parameters.openff_unconstrained_2_1_1(m0).getMolecule()
    m1 = BSS.Parameters.openff_unconstrained_2_1_1(m1).getMolecule()

    m0 = BSS.Align.rmsdAlign(m0, m1, mapping)

    BSS.Align.merge(
        m0, m1, mapping, allow_ring_breaking=True, allow_ring_size_change=True
    )


@pytest.mark.skipif(
    not has_antechamber or not has_tleap,
    reason="Requires antechamber and tLEaP to be installed.",
)
@pytest.mark.skipif(
    "openmm" not in sr.convert.supported_formats(),
    reason="Requires OpenMM to be installed.",
)
def test_ring_breaking_intrascale():
    """
    Test that ring-breaking merges produce correct intrascale matrices for a
    standard force field (GAFF2) with no non-default per-pair scale factors.

    The intrascale matrices are built from per-state connectivity, which
    correctly captures bonded distances across the ring-closure bond in the
    merged atom space. Since GAFF2 has no non-default per-pair scale factors,
    patchIntrascale is a no-op — verified by checking that apply_scale_factors=False
    gives the same changed bond and exception counts as the default path.
    """
    # Parameterise both molecules with GAFF2.
    cyclopentane = BSS.Parameters.gaff2("C1CCCC1").getMolecule()
    cyclohexane = BSS.Parameters.gaff2("C1CCCCC1").getMolecule()

    # Atom mapping for cyclopentane -> cyclohexane.
    mapping = {
        0: 0,
        1: 1,
        2: 2,
        3: 3,
        4: 4,
        5: 6,
        6: 7,
        7: 8,
        8: 9,
        9: 10,
        10: 11,
        11: 12,
        12: 13,
        13: 14,
        14: 15,
    }

    # Load the reference merged system produced by the old BioSimSpace merge
    # code (which used CLJNBPairs(connectivity, sf14) and was known correct).
    reference = sr.load_test_files("cyclopentane_cyclohexane.bss")
    ref_mol = reference.molecules("molecule property is_perturbable")[0]
    ref_omm = ref_mol.perturbation().to_openmm(map={"coordinates": "coordinates0"})
    ref_bonds = ref_omm.changed_bonds()
    ref_exceptions = ref_omm.changed_exceptions()

    # Merge cyclopentane -> cyclohexane and check against the reference.
    cyclopentane_aligned = BSS.Align.rmsdAlign(cyclopentane, cyclohexane, mapping)
    merged_fwd = BSS.Align.merge(
        cyclopentane_aligned,
        cyclohexane,
        mapping,
        allow_ring_size_change=True,
        allow_ring_breaking=True,
    )
    omm_fwd = merged_fwd._sire_object.perturbation().to_openmm(
        map={"coordinates": "coordinates0"}
    )
    assert len(omm_fwd.changed_bonds()) == len(ref_bonds)
    assert len(omm_fwd.changed_exceptions()) == len(ref_exceptions)

    # Merge in the reverse direction (cyclohexane -> cyclopentane) to verify
    # symmetry: the fix must hold regardless of which molecule is mol0/mol1.
    inv_mapping = {v: k for k, v in mapping.items()}
    cyclohexane_aligned = BSS.Align.rmsdAlign(cyclohexane, cyclopentane, inv_mapping)
    merged_rev = BSS.Align.merge(
        cyclohexane_aligned,
        cyclopentane,
        inv_mapping,
        allow_ring_size_change=True,
        allow_ring_breaking=True,
    )
    omm_rev = merged_rev._sire_object.perturbation().to_openmm(
        map={"coordinates": "coordinates0"}
    )
    assert len(omm_rev.changed_bonds()) == len(ref_bonds)
    assert len(omm_rev.changed_exceptions()) == len(ref_exceptions)

    # Verify patchIntrascale is a no-op for GAFF2: skipping it with
    # apply_scale_factors=False must give identical results.
    merged_nopatch = BSS.Align.merge(
        cyclopentane_aligned,
        cyclohexane,
        mapping,
        allow_ring_size_change=True,
        allow_ring_breaking=True,
        apply_scale_factors=False,
    )
    omm_nopatch = merged_nopatch._sire_object.perturbation().to_openmm(
        map={"coordinates": "coordinates0"}
    )
    assert len(omm_nopatch.changed_bonds()) == len(ref_bonds)
    assert len(omm_nopatch.changed_exceptions()) == len(ref_exceptions)


def test_ring_breaking_intrascale_m338():
    """
    Test that ring-breaking merges produce correct intrascale matrices for a
    real-world perturbation (int1 -> m338) with a standard force field (OpenFF).

    Since OpenFF has no non-default per-pair scale factors, patchIntrascale is
    a no-op and the merged intrascale matrices must exactly match those built
    directly from CLJNBPairs(conn0/conn1, sf14).
    """
    from sire.legacy import CAS as _SireCAS
    from sire.legacy import MM as _SireMM
    from sire.legacy import Mol as _SireMol

    # Atom mapping: {int1_idx: m338_idx}
    mapping = {
        21: 0,
        0: 1,
        23: 2,
        18: 3,
        1: 4,
        2: 5,
        3: 6,
        4: 7,
        5: 8,
        6: 9,
        7: 10,
        8: 11,
        9: 12,
        10: 13,
        19: 14,
        11: 15,
        12: 16,
        13: 17,
        14: 18,
        15: 19,
        16: 20,
        20: 21,
        30: 22,
        24: 23,
        22: 24,
        26: 25,
        17: 26,
        25: 27,
        27: 28,
        32: 29,
        33: 30,
        34: 31,
        35: 32,
        36: 33,
        37: 34,
        28: 35,
        38: 36,
    }

    mol0 = BSS.IO.readMolecules([f"{url}/int1.prm7", f"{url}/int1.rst7"])[0]
    mol1 = BSS.IO.readMolecules([f"{url}/m338.prm7", f"{url}/m338.rst7"])[0]

    mol0_aligned = BSS.Align.rmsdAlign(mol0, mol1, mapping)
    merged = BSS.Align.merge(
        mol0_aligned,
        mol1,
        mapping,
        allow_ring_breaking=True,
        allow_ring_size_change=True,
    )

    sire_mol = merged._sire_object

    intra0 = sire_mol.property("intrascale0")
    intra1 = sire_mol.property("intrascale1")

    # Build reference matrices directly from per-state connectivity. For OpenFF
    # these must be identical to the merge output since patchIntrascale is a no-op.
    ff = mol0._sire_object.property("forcefield")
    sf14 = _SireMM.CLJScaleFactor(
        ff.electrostatic14_scale_factor(), ff.vdw14_scale_factor()
    )

    conn0_edit = _SireMol.Connectivity(sire_mol.info()).edit()
    conn1_edit = _SireMol.Connectivity(sire_mol.info()).edit()
    for bond in sire_mol.property("bond0").potentials():
        ab = _SireMM.AmberBond(bond.function(), _SireCAS.Symbol("r"))
        if ab.k() != 0.0:
            conn0_edit.connect(bond.atom0(), bond.atom1())
    for bond in sire_mol.property("bond1").potentials():
        ab = _SireMM.AmberBond(bond.function(), _SireCAS.Symbol("r"))
        if ab.k() != 0.0:
            conn1_edit.connect(bond.atom0(), bond.atom1())

    ref_intra0 = _SireMM.CLJNBPairs(conn0_edit.commit(), sf14)
    ref_intra1 = _SireMM.CLJNBPairs(conn1_edit.commit(), sf14)

    # The two approaches must agree on every atom pair.
    n = sire_mol.num_atoms()
    for i in range(n):
        for j in range(i, n):
            idx_i = _SireMol.AtomIdx(i)
            idx_j = _SireMol.AtomIdx(j)
            assert intra0.get(idx_i, idx_j).coulomb() == pytest.approx(
                ref_intra0.get(idx_i, idx_j).coulomb()
            ), f"intra0 coulomb mismatch at ({i},{j})"
            assert intra0.get(idx_i, idx_j).lj() == pytest.approx(
                ref_intra0.get(idx_i, idx_j).lj()
            ), f"intra0 lj mismatch at ({i},{j})"
            assert intra1.get(idx_i, idx_j).coulomb() == pytest.approx(
                ref_intra1.get(idx_i, idx_j).coulomb()
            ), f"intra1 coulomb mismatch at ({i},{j})"
            assert intra1.get(idx_i, idx_j).lj() == pytest.approx(
                ref_intra1.get(idx_i, idx_j).lj()
            ), f"intra1 lj mismatch at ({i},{j})"


@pytest.mark.skipif(
    not has_antechamber or not has_tleap,
    reason="Requires antechamber and tLEaP to be installed.",
)
@pytest.mark.skipif(
    not has_openff,
    reason="Requires OpenFF to be installed.",
)
def test_ring_breaking_cross_bond_cleanup():
    """
    Test that bonded terms spanning a ring-breaking bond are removed from the
    end-state properties where that bond is absent (SYK 5035→5033).

    In this perturbation a ring opens, leaving one bond present at λ=0 but
    absent at λ=1. Any angle, dihedral or improper whose geometry depends on
    that bond must be removed from the λ=1 properties; retaining them would
    constrain atoms toward a bonded geometry that no longer exists and cause
    large repulsion at the nonbonded/bonded lambda boundary.
    """

    # MCS mapping: {5033_idx: 5035_idx} — mol0=5033 (ring present), mol1=5035
    # (ring absent), so the ring bond appears in connectivity0 but not
    # connectivity1, giving ring_breaking = {(1, 7)} in the merged molecule.
    mapping = {
        6: 0,
        5: 1,
        4: 2,
        3: 3,
        34: 4,
        8: 5,
        32: 6,
        31: 7,
        11: 8,
        12: 9,
        13: 10,
        14: 11,
        15: 12,
        16: 13,
        17: 14,
        18: 15,
        19: 16,
        20: 17,
        21: 18,
        22: 19,
        23: 20,
        24: 21,
        25: 22,
        26: 23,
        27: 24,
        28: 25,
        29: 26,
        30: 27,
        10: 28,
        9: 29,
        7: 30,
        38: 31,
        37: 32,
        35: 33,
        36: 34,
        1: 35,
        33: 36,
        53: 38,
        52: 39,
        43: 40,
        44: 41,
        45: 42,
        46: 43,
        47: 44,
        48: 45,
        49: 46,
        50: 47,
        51: 48,
        42: 49,
        41: 50,
    }

    mol0 = BSS.Parameters.openff_unconstrained_2_2_1(
        BSS.IO.readMolecules(f"{url}/5033.sdf")[0]
    ).getMolecule()
    mol1 = BSS.Parameters.openff_unconstrained_2_2_1(
        BSS.IO.readMolecules(f"{url}/5035.sdf")[0]
    ).getMolecule()

    mol0_aligned = BSS.Align.rmsdAlign(mol0, mol1, mapping)
    merged = BSS.Align.merge(
        mol0_aligned,
        mol1,
        mapping,
        allow_ring_breaking=True,
    )

    sire_mol = merged._sire_object
    mol_info = sire_mol.info()

    conn0 = sire_mol.property("connectivity0")
    conn1 = sire_mol.property("connectivity1")

    bonds0 = {
        (
            min(b.atom0().value(), b.atom1().value()),
            max(b.atom0().value(), b.atom1().value()),
        )
        for b in conn0.get_bonds()
    }
    bonds1 = {
        (
            min(b.atom0().value(), b.atom1().value()),
            max(b.atom0().value(), b.atom1().value()),
        )
        for b in conn1.get_bonds()
    }

    # Bonds present only at λ=1 must not appear in angle0/dihedral0/improper0.
    ring_making = bonds1 - bonds0
    # Bonds present only at λ=0 must not appear in angle1/dihedral1/improper1.
    ring_breaking = bonds0 - bonds1

    # This perturbation must have at least one ring-breaking bond.
    assert ring_breaking, "Expected ring-breaking bonds in SYK 5035→5033"

    for changing, suffix in [(ring_making, "0"), (ring_breaking, "1")]:
        if not changing:
            continue

        for p in sire_mol.property(f"angle{suffix}").potentials():
            i = mol_info.atom_idx(p.atom0()).value()
            j = mol_info.atom_idx(p.atom1()).value()
            k = mol_info.atom_idx(p.atom2()).value()
            assert (min(i, j), max(i, j)) not in changing, (
                f"angle{suffix} ({i},{j},{k}) spans absent bond "
                f"({min(i, j)},{max(i, j)})"
            )
            assert (min(j, k), max(j, k)) not in changing, (
                f"angle{suffix} ({i},{j},{k}) spans absent bond "
                f"({min(j, k)},{max(j, k)})"
            )

        for p in sire_mol.property(f"dihedral{suffix}").potentials():
            j = mol_info.atom_idx(p.atom1()).value()
            k = mol_info.atom_idx(p.atom2()).value()
            assert (
                min(j, k),
                max(j, k),
            ) not in changing, (
                f"dihedral{suffix} central bond ({j},{k}) spans absent bond"
            )

        for p in sire_mol.property(f"improper{suffix}").potentials():
            atoms = {
                mol_info.atom_idx(p.atom0()).value(),
                mol_info.atom_idx(p.atom1()).value(),
                mol_info.atom_idx(p.atom2()).value(),
                mol_info.atom_idx(p.atom3()).value(),
            }
            for a, b in changing:
                assert not (
                    a in atoms and b in atoms
                ), f"improper{suffix} spans absent bond ({a},{b})"

    # Check that the ring-breaking and ring-making bond properties are set.
    def _read_pairs(prop_name):
        if not sire_mol.has_property(prop_name):
            return set()
        flat = list(sire_mol.property(prop_name).to_list())
        return {(flat[i], flat[i + 1]) for i in range(0, len(flat), 2)}

    stored_breaking = _read_pairs("ring_breaking_bonds")
    stored_making = _read_pairs("ring_making_bonds")
    assert (
        stored_breaking == ring_breaking
    ), f"ring_breaking_bonds property mismatch: {stored_breaking} != {ring_breaking}"
    assert (
        stored_making == ring_making
    ), f"ring_making_bonds property mismatch: {stored_making} != {ring_making}"
