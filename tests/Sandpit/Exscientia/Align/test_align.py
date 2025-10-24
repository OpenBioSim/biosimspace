import pytest
import sys

from sire.legacy.MM import InternalFF, IntraCLJFF, IntraFF
from sire.legacy.Mol import AtomIdx, Element, PartialMolecule
from sire.legacy.Maths import Vector

import BioSimSpace.Sandpit.Exscientia as BSS

from tests.Sandpit.Exscientia.conftest import url, has_antechamber, has_openff


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


@pytest.fixture()
def propane():
    return BSS.Parameters.parameterise(
        "CCC", "openff_unconstrained-2.0.0"
    ).getMolecule()


@pytest.fixture()
def butane():
    return BSS.Parameters.parameterise(
        "CCCC", "openff_unconstrained-2.0.0"
    ).getMolecule()


@pytest.fixture()
def propane_butane(propane, butane):
    mapping = BSS.Align.matchAtoms(propane, butane)
    # We make sure we have a dummy atom in both endstates
    mapping.pop(3)
    return BSS.Align.merge(propane, butane, mapping)


@pytest.mark.skipif(
    has_antechamber is False or has_openff is False,
    reason="Requires AmberTools/antechamber and OpenFF to be installed.",
)
def test_propane_butane_1_4(propane_butane, tmp_path_factory):
    # This is a regression test for a bug introduced in a recent commit
    # I don't know of a more rational way of checking it so I am doing it the old-fashioned way
    tmpdir = tmp_path_factory.mktemp("TestMerge")
    BSS.IO.saveMolecules(str(tmpdir / "propane_butane"), propane_butane, "grotop")
    contents = open(tmpdir / "propane_butane.top").readlines()
    n_nb = 0
    reached_nb = False

    # Get the number of 1-4 interactions
    for line in contents:
        if line[0] == ";":
            continue
        if reached_nb:
            if line.strip():
                n_nb += 1
            else:
                break
        if "[ pairs ]" in line:
            reached_nb = True

    assert n_nb == 33


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

    # Test that the intramolecular energies are the same.

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


@pytest.fixture(scope="module")
def roi_mol0():
    files = BSS.IO.expand(url, ["wild.prmtop", "wild.inpcrd"], ".bz2")
    return BSS.IO.readMolecules(files)[0]


@pytest.fixture(scope="module")
def roi_mol1():
    files = BSS.IO.expand(url, ["mutated.prmtop", "mutated.inpcrd"], ".bz2")
    return BSS.IO.readMolecules(files)[0]


@pytest.fixture(scope="module")
def roi_merged_mol(roi_mol0, roi_mol1):
    # Extract the residue name and index
    def get_res_info(mol):
        mol_res_name = []
        mol_res_idx = []
        for bss_res in mol.getResidues():
            mol_res_name.append(bss_res.name())
            mol_res_idx.append([x.index() for x in bss_res.getAtoms()])
        return mol_res_name, mol_res_idx

    def recover_mapping(mut1_idx, mut2_idx, mapping):
        new_mapping = {}
        for k, v in mapping.items():
            new_mapping[mut1_idx[k]] = mut2_idx[v]
        return new_mapping

    def update_coordinate(mol, coord_dict):
        edit_mol = mol._sire_object.edit()
        for idx in range(0, mol.nAtoms()):
            if idx in coord_dict:
                vec = coord_dict[idx]
                vec = Vector(
                    vec.x().angstroms().value(),
                    vec.y().angstroms().value(),
                    vec.z().angstroms().value(),
                )
                edit_mol = (
                    edit_mol.atom(AtomIdx(idx))
                    .set_property("coordinates", vec)
                    .molecule()
                )
        mol._sire_object = edit_mol.commit()

    mol0_res_name, mol0_res_idx = get_res_info(roi_mol0)
    mol1_res_name, mol1_res_idx = get_res_info(roi_mol1)

    # Get the best mapping between the molecules.
    mapping = {}
    for i in range(len(mol0_res_name)):
        # The mutated residue would have different name
        if mol0_res_name[i] == mol1_res_name[i]:
            for k, v in zip(mol0_res_idx[i], mol1_res_idx[i]):
                mapping[k] = v
        else:
            # Assume the atom order is kept during conversion
            mut0 = roi_mol0.getResidues()[i].toMolecule()
            mut1 = roi_mol1.getResidues()[i].toMolecule()
            mut0_idx = [x.index() for x in roi_mol0.getResidues()[i].getAtoms()]
            mut1_idx = [x.index() for x in roi_mol1.getResidues()[i].getAtoms()]

            mapping_mut = BSS.Align.matchAtoms(mut0, mut1)

            # We should also translate the moelcule according to the mapping
            # or there will be odd dummy atom in the merged topology
            mut0 = BSS.Align.rmsdAlign(mut0, mut1, mapping=mapping_mut)
            aligned_coordinates = mut0.coordinates()
            coord_dict = dict(zip(mut0_idx, aligned_coordinates))

            mapping_recovered = recover_mapping(mut0_idx, mut1_idx, mapping_mut)
    mapping.update(mapping_recovered)
    update_coordinate(roi_mol0, coord_dict)

    # Create the merged molecule.
    merged_mol = BSS.Align.merge(
        roi_mol0, roi_mol1, mapping=mapping, roi=[mut0_idx, mut1_idx]
    )

    return merged_mol


# IntraCLJFF:
#  Old interface. Uses the "intrascale" matrix. Validate that this
#  is consistent.
# IntraFF:
#  New interface. Uses atom "connectivity". Validate that the bonding
#  is consistent.
@pytest.fixture(scope="module")
def roi_intraclj0(roi_mol0):
    res = IntraCLJFF("intraclj")
    res.add(roi_mol0._sire_object)
    return res


@pytest.fixture(scope="module")
def roi_intraclj1(roi_mol1):
    res = IntraCLJFF("intraclj")
    res.add(roi_mol1._sire_object)
    return res


@pytest.fixture(scope="module")
def roi_intraff0(roi_mol0):
    res = IntraFF("intraclj")
    res.add(roi_mol0._sire_object)
    return res


@pytest.fixture(scope="module")
def roi_intraff1(roi_mol1):
    res = IntraFF("intraclj")
    res.add(roi_mol1._sire_object)
    return res


@pytest.fixture(scope="module")
def roi_internal0(roi_mol0):
    res = InternalFF("internal")
    res.setStrict(True)
    res.add(roi_mol0._sire_object)
    return res


@pytest.fixture(scope="module")
def roi_internal1(roi_mol1):
    res = InternalFF("internal")
    res.setStrict(True)
    res.add(roi_mol1._sire_object)
    return res


@pytest.fixture(scope="module")
def roi_pmap0(roi_merged_mol):
    res = {}
    for prop in roi_merged_mol._sire_object.property_keys():
        if prop[-1] == "0":
            res[prop[:-1]] = prop
    return res


@pytest.fixture(scope="module")
def roi_pmap1(roi_merged_mol):
    res = {}
    for prop in roi_merged_mol._sire_object.property_keys():
        if prop[-1] == "1":
            res[prop[:-1]] = prop
    return res


def test_roi_nonbonded0(roi_merged_mol, roi_intraclj0, roi_intraff0, roi_pmap0):
    # Test that the nonbonded energies are the same in mol0.
    intraclj_merged = IntraCLJFF("intraclj")
    intraff_merged = IntraFF("intraclj")
    intraclj_merged.add(roi_merged_mol._sire_object, roi_pmap0)
    intraff_merged.add(roi_merged_mol._sire_object, roi_pmap0)
    assert roi_intraclj0.energy().value() == pytest.approx(
        intraclj_merged.energy().value()
    )
    assert roi_intraff0.energy().value() == pytest.approx(
        intraff_merged.energy().value()
    )


def test_roi_nonbonded1(roi_merged_mol, roi_intraclj1, roi_intraff1, roi_pmap1):
    # Test that the nonbonded energies are the same in mol1.
    intraclj_merged = IntraCLJFF("intraclj")
    intraff_merged = IntraFF("intraclj")
    intraclj_merged.add(roi_merged_mol._sire_object, roi_pmap1)
    intraff_merged.add(roi_merged_mol._sire_object, roi_pmap1)
    assert roi_intraclj1.energy().value() == pytest.approx(
        intraclj_merged.energy().value()
    )
    assert roi_intraff1.energy().value() == pytest.approx(
        intraff_merged.energy().value()
    )


def test_roi_bonded0(roi_mol0, roi_merged_mol, roi_pmap0, roi_internal0):
    # Test that the internal energies are consistent. This will validate that
    # bond, angle, dihedral, and improper energies are correct.
    # In this test, we extract the original molecule for the lambda=0 end state.
    amber_mol, _ = roi_merged_mol._extractMolecule()
    roi_internal_merged = InternalFF("internal")
    roi_internal_merged.setStrict(True)
    roi_internal_merged.add(amber_mol._sire_object)
    assert roi_internal0.energy().value() == pytest.approx(
        roi_internal_merged.energy().value()
    )


def test_roi_bonded1(roi_mol1, roi_merged_mol, roi_pmap1, roi_internal1):
    # Test that the internal energies are consistent. This will validate that
    # bond, angle, dihedral, and improper energies are correct.
    # In this test, we extract the original molecule for the lambda=1 end state.
    amber_mol, _ = roi_merged_mol._extractMolecule(is_lambda1=True)
    roi_internal_merged = InternalFF("internal")
    roi_internal_merged.setStrict(True)
    roi_internal_merged.add(amber_mol._sire_object)
    assert roi_internal1.energy().value() == pytest.approx(
        roi_internal_merged.energy().value()
    )


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


def test_ion_merge():
    from sire.legacy.IO import createSodiumIon
    from tests.conftest import root_fp

    # Extract a water molecule from the system.
    water = BSS.IO.readMolecules(
        [f"{root_fp}/input/ala.crd", f"{root_fp}/input/ala.top"]
    )[-1]

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
