import pytest

import BioSimSpace as BSS

from tests.conftest import url, has_namd

# Store the allowed restraints.
restraints = BSS.Protocol._position_restraint_mixin._PositionRestraintMixin.restraints()


@pytest.fixture(scope="module")
def namd_system():
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules(
        [
            "tests/input/alanin.psf",
            f"tests/input/alanin.pdb",
            f"tests/input/alanin.params",
        ]
    )


@pytest.mark.skipif(has_namd is False, reason="Requires NAMD to be installed.")
@pytest.mark.parametrize("restraint", restraints)
def test_minimise(namd_system, restraint):
    """Test a minimisation protocol."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100, restraint=restraint)

    # Run the process, check that it finished without error, and returns a system.
    run_process(namd_system, protocol)


@pytest.mark.skipif(has_namd is False, reason="Requires NAMD to be installed.")
@pytest.mark.parametrize("restraint", restraints)
def test_equilibrate(namd_system, restraint):
    """Test an equilibration protocol."""

    # Create a short equilibration protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"), restraint=restraint
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(namd_system, protocol)


@pytest.mark.skipif(has_namd is False, reason="Requires NAMD to be installed.")
def test_heat(namd_system):
    """Test a heating protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"),
        temperature_start=BSS.Types.Temperature(0, "kelvin"),
        temperature_end=BSS.Types.Temperature(300, "kelvin"),
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(namd_system, protocol)


@pytest.mark.skipif(has_namd is False, reason="Requires NAMD to be installed.")
def test_cool(namd_system):
    """Test a cooling protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"),
        temperature_start=BSS.Types.Temperature(300, "kelvin"),
        temperature_end=BSS.Types.Temperature(0, "kelvin"),
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(namd_system, protocol)


@pytest.mark.skipif(has_namd is False, reason="Requires NAMD to be installed.")
@pytest.mark.parametrize("restraint", restraints)
def test_production(namd_system, restraint):
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(
        runtime=BSS.Types.Time(0.001, "nanoseconds"), restraint=restraint
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(namd_system, protocol)


@pytest.mark.skipif(has_namd is False, reason="Requires NAMD to be installed.")
def test_perturbable_restraint(perturbable_system):
    """Test a free energy perturbation protocol."""

    # Create a short minimisation prototocol with a restraint.
    protocol = BSS.Protocol.Minimisation(steps=100, restraint="heavy")

    # Run the process, check that it finished without error, and returns a system.
    run_process(perturbable_system, protocol)


def run_process(namd_system, protocol):
    """Helper function to run various simulation protocols."""

    # Initialise the NAMD process.
    process = BSS.Process.Namd(namd_system, protocol, name="test")

    # Start the NAMD simulation.
    process.start()

    # Wait for the process to end.
    process.wait()

    # Make sure the process didn't error.
    assert not process.isError()

    # Make sure that we get a molecular system back.
    assert process.getSystem() is not None
