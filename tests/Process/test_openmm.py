import pytest

import BioSimSpace as BSS

from tests.conftest import url, has_amber, has_gromacs, has_openff

# Store the allowed restraints.
restraints = BSS.Protocol._position_restraint_mixin._PositionRestraintMixin.restraints()


@pytest.fixture(scope="session")
def system():
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules(["tests/input/ala.top", "tests/input/ala.crd"])


@pytest.mark.parametrize("restraint", restraints)
def test_minimise(system, restraint):
    """Test a minimisation protocol."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100, restraint=restraint)

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.parametrize("restraint", restraints)
def test_equilibrate(system, restraint):
    """Test an equilibration protocol."""

    # Create a short equilibration protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"), restraint=restraint
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


def test_heat(system):
    """Test a heating protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"),
        temperature_start=BSS.Types.Temperature(0, "kelvin"),
        temperature_end=BSS.Types.Temperature(300, "kelvin"),
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


def test_cool(system):
    """Test a cooling protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"),
        temperature_start=BSS.Types.Temperature(300, "kelvin"),
        temperature_end=BSS.Types.Temperature(0, "kelvin"),
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.parametrize("restraint", restraints)
def test_production(system, restraint):
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(
        runtime=BSS.Types.Time(0.001, "nanoseconds"), restraint=restraint
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.skipif(
    has_amber is False or has_gromacs is False or has_openff is False,
    reason="Requires AMBER, GROMACS, and OpenFF to be installed",
)
def test_rhombic_dodecahedron():
    """Test that OpenMM can load and run rhombic dodecahedral triclinic spaces."""

    # Create a methane molecule.
    mol = BSS.Parameters.openff_unconstrained_2_0_0("C").getMolecule()

    # Generate box dimensions and angles for a hexagonal rhombic dodecahedron.
    box, angles = BSS.Box.rhombicDodecahedronHexagon(5 * BSS.Units.Length.nanometer)

    # Create a solvated system.
    solvated = BSS.Solvent.tip3p(mol, box=box, angles=angles)

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Run the process, check that it finished without error, and returns a system.
    run_process(solvated, protocol)


def test_parmed_triclinic():
    """Test the workaround for fixing ParmEd triclinic lattice reduction."""

    # Load the test system.
    system = BSS.IO.readMolecules(
        BSS.IO.expand(url, ["parmed_issue.rst7", "parmed_issue.prm7"])
    )

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


def run_process(system, protocol):
    """Helper function to run various simulation protocols."""

    # Initialise the OpenMM  process.
    process = BSS.Process.OpenMM(system, protocol, name="test")

    # Start the OpenMM simulation.
    process.start()

    # Wait for the process to end.
    process.wait()

    # Make sure the process didn't error.
    assert not process.isError()

    # Make sure that we get a molecular system back.
    assert process.getSystem() is not None
