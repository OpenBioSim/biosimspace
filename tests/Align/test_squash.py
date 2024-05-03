import os
import numpy as np
import pickle
import pytest

import sire

from sire.maths import Vector

import BioSimSpace as BSS

# Make sure AMBER is installed.
if BSS._amber_home is not None:
    exe = "%s/bin/sander" % BSS._amber_home
    if os.path.isfile(exe):
        has_amber = True
    else:
        has_amber = False
else:
    has_amber = False


@pytest.fixture(scope="session")
def perturbed_system():
    # N_atoms are: 12, 15, 18, 21, 24, 27 and 30.
    mol_smiles = [
        "c1ccccc1",
        "c1ccccc1C",
        "c1ccccc1CC",
        "c1ccccc1CCC",
        "c1ccccc1CCCC",
        "c1ccccc1CCCCC",
        "c1ccccc1CCCCCC",
    ]
    mols = [BSS.Parameters.gaff(smi).getMolecule() for smi in mol_smiles]
    pert_mols = [
        mols[0],
        BSS.Align.merge(mols[1], mols[2]),
        mols[3],
        mols[4],
        BSS.Align.merge(mols[5], mols[6]),
    ]
    system = BSS._SireWrappers.System(pert_mols)
    return system


@pytest.fixture(scope="session")
def dual_topology_system():
    mol_smiles = ["c1ccccc1", "c1ccccc1C"]
    mols = [BSS.Parameters.gaff(smi).getMolecule() for smi in mol_smiles]
    pertmol = BSS.Align.merge(mols[0], mols[1], mapping={0: 0})
    c = pertmol._sire_object.cursor()
    # Translate all atoms so that we have different coordinates between both endstates
    for atom in c.atoms():
        atom["coordinates1"] = atom["coordinates0"] + Vector(1, 1, 1)
    pertmol = BSS._SireWrappers.Molecule(c.commit())
    system = pertmol.toSystem()
    return system


@pytest.fixture
def perturbed_tripeptide():
    return pickle.load(open(f"tests/input/merged_tripeptide.pickle", "rb"))


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
@pytest.mark.parametrize(
    "explicit,expected_n_atoms",
    [
        (False, [12, 21, 24, 15, 18, 27, 30]),
        (True, [12, 21, 24, 18, 18, 30, 30]),
    ],
)
def test_squash(perturbed_system, explicit, expected_n_atoms):
    squashed_system, mapping = BSS.Align._squash._squash(
        perturbed_system, explicit_dummies=explicit
    )
    assert len(squashed_system) == 7
    n_atoms = [mol.nAtoms() for mol in squashed_system]
    assert squashed_system[-2].getResidues()[0].name() == "LIG"
    assert squashed_system[-1].getResidues()[0].name() == "LIG"
    # First we must have the unperturbed molecules, and then the perturbed ones.
    assert n_atoms == expected_n_atoms
    python_mapping = {k.value(): v.value() for k, v in mapping.items()}
    assert python_mapping == {0: 0, 2: 1, 3: 2}


@pytest.mark.parametrize("explicit", [False, True])
def test_squash_multires(perturbed_tripeptide, explicit):
    squashed_system, mapping = BSS.Align._squash._squash(
        perturbed_tripeptide, explicit_dummies=explicit
    )
    assert len(squashed_system) == 1
    assert len(squashed_system[0].getResidues()) == 4


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
@pytest.mark.parametrize("is_lambda1", [False, True])
def test_squashed_molecule_mapping(perturbed_system, is_lambda1):
    res = BSS.Align._squash._squashed_molecule_mapping(
        perturbed_system, is_lambda1=is_lambda1
    )
    if not is_lambda1:
        expected = {0: 0, 2: 1, 3: 2, 1: 3, 4: 5}
    else:
        expected = {0: 0, 2: 1, 3: 2, 1: 4, 4: 6}
    assert res == expected


@pytest.mark.parametrize("is_lambda1", [False, True])
def test_squashed_atom_mapping_implicit(perturbed_tripeptide, is_lambda1):
    res = BSS.Align._squash._squashed_atom_mapping(
        perturbed_tripeptide, is_lambda1=is_lambda1, explicit_dummies=False
    )
    if not is_lambda1:
        merged_indices = list(range(16)) + list(range(16, 30)) + list(range(43, 51))
        squashed_indices = list(range(16)) + list(range(16, 30)) + list(range(30, 38))
    else:
        merged_indices = (
            list(range(16))
            + list(range(16, 21))
            + list(range(23, 26))
            + list(range(30, 43))
            + list(range(43, 51))
        )
        squashed_indices = list(range(16)) + list(range(38, 59)) + list(range(30, 38))
    expected = dict(zip(merged_indices, squashed_indices))
    assert res == expected


@pytest.mark.parametrize("is_lambda1", [False, True])
def test_squashed_atom_mapping_explicit(perturbed_tripeptide, is_lambda1):
    res = BSS.Align._squash._squashed_atom_mapping(
        perturbed_tripeptide, is_lambda1=is_lambda1, explicit_dummies=True
    )
    merged_indices = list(range(51))
    if not is_lambda1:
        squashed_indices = list(range(16)) + list(range(16, 43)) + list(range(43, 51))
    else:
        squashed_indices = list(range(16)) + list(range(51, 78)) + list(range(43, 51))
    expected = dict(zip(merged_indices, squashed_indices))
    assert res == expected


@pytest.mark.skipif(has_amber is False, reason="Requires AMBER to be installed.")
@pytest.mark.parametrize("explicit", [False, True])
def test_unsquash(dual_topology_system, explicit):
    squashed_system, mapping = BSS.Align._squash._squash(
        dual_topology_system, explicit_dummies=explicit
    )
    new_perturbed_system = BSS.Align._squash._unsquash(
        dual_topology_system, squashed_system, mapping, explicit_dummies=explicit
    )
    assert [
        mol0.nAtoms() == mol1.nAtoms()
        for mol0, mol1 in zip(dual_topology_system, new_perturbed_system)
    ]
    assert [
        mol0.isPerturbable() == mol1.isPerturbable()
        for mol0, mol1 in zip(dual_topology_system, new_perturbed_system)
    ]
    if explicit:
        # Check that we have loaded the correct coordinates
        coords0_before = sire.io.get_coords_array(
            dual_topology_system[0]._sire_object, map={"coordinates": "coordinates0"}
        )
        coords1_before = sire.io.get_coords_array(
            dual_topology_system[0]._sire_object, map={"coordinates": "coordinates1"}
        )
        coords0_after = sire.io.get_coords_array(
            new_perturbed_system[0]._sire_object, map={"coordinates": "coordinates0"}
        )
        coords1_after = sire.io.get_coords_array(
            new_perturbed_system[0]._sire_object, map={"coordinates": "coordinates1"}
        )

        # The coordinates at the first endstate should be completely preserved
        # Because in this case they are either common core, or separate dummies at lambda = 0
        assert np.allclose(coords0_before, coords0_after)

        # The coordinates at the first endstate should be partially preserved
        # The common core must have the same coordinates as lambda = 0
        # Here this is just a single atom in the beginning
        # The extra atoms which are dummies at lambda = 0 have separate coordinates here
        assert np.allclose(coords0_before[:1, :], coords1_after[:1, :])
        assert np.allclose(coords1_before[1:, :], coords1_after[1:, :])


@pytest.mark.parametrize("explicit", [False, True])
def test_unsquash_multires(perturbed_tripeptide, explicit):
    squashed_system, mapping = BSS.Align._squash._squash(
        perturbed_tripeptide, explicit_dummies=explicit
    )
    new_perturbed_system = BSS.Align._squash._unsquash(
        perturbed_tripeptide, squashed_system, mapping, explicit_dummies=explicit
    )
    assert [
        mol0.nAtoms() == mol1.nAtoms()
        for mol0, mol1 in zip(perturbed_tripeptide, new_perturbed_system)
    ]
    assert [
        mol0.isPerturbable() == mol1.isPerturbable()
        for mol0, mol1 in zip(perturbed_tripeptide, new_perturbed_system)
    ]
