import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from tests.Sandpit.Exscientia.conftest import has_gromacs


@pytest.fixture(scope="module")
def solvated_system():
    """
    A BSS-solvated alanine dipeptide system with no existing ions.

    Uses BSS.Solvent.tip3p so water is named SOL (as required by genion).
    """
    if not has_gromacs:
        pytest.skip("Requires GROMACS to be installed")
    system = BSS.IO.readMolecules(["tests/input/ala.top", "tests/input/ala.crd"])
    box, angles = BSS.Box.cubic(4 * BSS.Units.Length.nanometer)
    return BSS.Solvent.tip3p(molecule=system[0], box=box, angles=angles)


@pytest.fixture(scope="module")
def solvated_system_with_ions():
    """
    A BSS-solvated alanine dipeptide system that already contains Na+/Cl- ions.

    Uses BSS.Solvent.tip3p so water is named SOL (as required by genion).
    """
    if not has_gromacs:
        pytest.skip("Requires GROMACS to be installed")
    system = BSS.IO.readMolecules(["tests/input/ala.top", "tests/input/ala.crd"])
    box, angles = BSS.Box.cubic(4 * BSS.Units.Length.nanometer)
    return BSS.Solvent.tip3p(molecule=system[0], box=box, angles=angles, ion_conc=0.15)


def test_ions():
    """Test that ions() returns the expected list of supported ion names."""
    ion_list = BSS.Solvent.ions()
    assert isinstance(ion_list, list)
    assert len(ion_list) > 0
    # Check a representative selection of expected ions are present.
    for ion in ["br", "ca", "cl", "cs", "f", "k", "li", "mg", "na", "rb", "zn"]:
        assert ion in ion_list


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed")
def test_add_ions_no_existing_ions(solvated_system):
    """Test adding Mg2+ ions to a solvated system that has no existing ions."""
    num_ions = 2
    num_water_before = solvated_system.nWaterMolecules()

    result = BSS.Solvent.addIons(solvated_system, "mg", num_ions=num_ions)

    # Check that the correct number of MG ions were added.
    try:
        mg_molecules = result.search("resname MG").molecules()
    except Exception:
        pytest.fail("No MG ions found in the result system.")

    assert len(mg_molecules) == num_ions

    # Each ion replaces one water molecule.
    assert result.nWaterMolecules() == num_water_before - num_ions


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed")
def test_add_ions_neutral(solvated_system):
    """
    Test that is_neutral=True causes genion to add counter-ions to neutralise
    the system charge.

    Adding 2 Mg2+ ions introduces a +4 charge; with is_neutral=True genion
    should add 4 Cl- counter-ions, resulting in a net system charge of zero.
    """
    num_ions = 2

    result = BSS.Solvent.addIons(
        solvated_system, "mg", num_ions=num_ions, is_neutral=True
    )

    # Check that MG ions were added.
    try:
        mg_molecules = result.search("resname MG").molecules()
    except Exception:
        pytest.fail("No MG ions found in the result system.")

    assert len(mg_molecules) == num_ions

    # 2 Mg2+ introduce a +4 charge; genion should add 4 Cl- to neutralise.
    try:
        cl_molecules = result.search("resname CL").molecules()
    except Exception:
        pytest.fail("No CL counter-ions found in the result system.")

    assert len(cl_molecules) == num_ions * 2

    # The total system charge should be zero.
    assert round(result.charge().value()) == 0


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed")
def test_add_ions_with_existing_ions(solvated_system_with_ions):
    """
    Test adding Mg2+ ions to a solvated system that already has Na+/Cl- ions.

    Verifies that existing ions are preserved after adding the new ion species.
    """
    num_ions = 2
    num_water_before = solvated_system_with_ions.nWaterMolecules()

    # Record the existing ion counts.
    try:
        num_na_before = len(solvated_system_with_ions.search("resname NA").molecules())
    except Exception:
        num_na_before = 0
    try:
        num_cl_before = len(solvated_system_with_ions.search("resname CL").molecules())
    except Exception:
        num_cl_before = 0

    result = BSS.Solvent.addIons(solvated_system_with_ions, "mg", num_ions=num_ions)

    # Check that the correct number of MG ions were added.
    try:
        mg_molecules = result.search("resname MG").molecules()
    except Exception:
        pytest.fail("No MG ions found in the result system.")

    assert len(mg_molecules) == num_ions

    # Each ion replaces one water molecule.
    assert result.nWaterMolecules() == num_water_before - num_ions

    # Existing Na+ and Cl- ions should be unchanged.
    try:
        num_na_after = len(result.search("resname NA").molecules())
    except Exception:
        num_na_after = 0
    try:
        num_cl_after = len(result.search("resname CL").molecules())
    except Exception:
        num_cl_after = 0

    assert num_na_after == num_na_before
    assert num_cl_after == num_cl_before
