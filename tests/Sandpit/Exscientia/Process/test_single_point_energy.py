import pytest

import BioSimSpace.Sandpit.Exscientia as BSS

from tests.Sandpit.Exscientia.conftest import url, has_amber, has_gromacs, has_pyarrow


@pytest.fixture(scope="session")
def system():
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules(
        [f"{url}/ubiquitin.prm7.bz2", f"{url}/ubiquitin.rst7.bz2"]
    )


@pytest.mark.skipif(
    has_amber is False or has_gromacs is False or has_pyarrow is False,
    reason="Requires that AMBER, GROMACS, and pyarrow are installed.",
)
def test_amber_gromacs(system):
    """Single point energy comparison between AMBER and GROMACS."""

    # Create a single-step minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=1)

    # Create a process to run with AMBER.
    process_amb = BSS.Process.Amber(system, protocol)

    # Create a process to run with GROMACS.
    process_gmx = BSS.Process.Gromacs(system, protocol, extra_options={"nsteps": 0})

    # Run the AMBER process and wait for it to finish.
    process_amb.start()
    process_amb.wait()

    # Run the GROMACS process and wait for it to finish.
    process_gmx.start()
    process_gmx.wait()

    # Compare bond energies. (In kJ / mol)
    nrg_amb = process_amb.getBondEnergy().kj_per_mol().value()
    nrg_gmx = process_gmx.getBondEnergy().kj_per_mol().value()
    assert nrg_amb == pytest.approx(nrg_gmx, rel=1e-2)

    # Compare angle energies. (In kJ / mol)
    nrg_amb = process_amb.getAngleEnergy().kj_per_mol().value()
    nrg_gmx = process_gmx.getAngleEnergy().kj_per_mol().value()
    assert nrg_amb == pytest.approx(nrg_gmx, rel=1e-2)

    # Compare dihedral energies. (In kJ / mol)
    nrg_amb = process_amb.getDihedralEnergy().kj_per_mol().value()
    nrg_gmx = process_gmx.getDihedralEnergy().kj_per_mol().value()
    assert nrg_amb == pytest.approx(nrg_gmx, rel=1e-2)


@pytest.mark.skipif(
    has_amber is False or has_gromacs is False or has_pyarrow is False,
    reason="Requires that AMBER, GROMACS, and pyarrow are installed.",
)
def test_amber_gromacs_triclinic(system):
    """Single point energy comparison between AMBER and GROMACS in a triclinic box."""

    # Swap the space for a triclinic cell (truncated octahedron).
    from sire.legacy.Vol import TriclinicBox

    triclinic_box = TriclinicBox.truncated_octahedron(50)
    system._sire_object.set_property("space", triclinic_box)

    # Create a single-step minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=1)

    # Create a process to run with AMBER.
    process_amb = BSS.Process.Amber(system, protocol)

    # Create a process to run with GROMACS.
    process_gmx = BSS.Process.Gromacs(system, protocol, extra_options={"nsteps": 0})

    # Run the AMBER process and wait for it to finish.
    process_amb.start()
    process_amb.wait()

    # Run the GROMACS process and wait for it to finish.
    process_gmx.start()
    process_gmx.wait()

    # Compare bond energies. (In kJ / mol)
    nrg_amb = process_amb.getBondEnergy().kj_per_mol().value()
    nrg_gmx = process_gmx.getBondEnergy().kj_per_mol().value()
    assert nrg_amb == pytest.approx(nrg_gmx, rel=1e-2)

    # Compare angle energies. (In kJ / mol)
    nrg_amb = process_amb.getAngleEnergy().kj_per_mol().value()
    nrg_gmx = process_gmx.getAngleEnergy().kj_per_mol().value()
    assert nrg_amb == pytest.approx(nrg_gmx, rel=1e-2)

    # Compare dihedral energies. (In kJ / mol)
    nrg_amb = process_amb.getDihedralEnergy().kj_per_mol().value()
    nrg_gmx = process_gmx.getDihedralEnergy().kj_per_mol().value()
    assert nrg_amb == pytest.approx(nrg_gmx, rel=1e-2)
