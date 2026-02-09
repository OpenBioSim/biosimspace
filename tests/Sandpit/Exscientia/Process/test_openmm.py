import math

import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from tests.conftest import root_fp
from tests.Sandpit.Exscientia.conftest import has_amber, has_gromacs, has_openff, url


@pytest.fixture(scope="session")
def system():
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules(
        [f"{root_fp}/input/ala.top", f"{root_fp}/input/ala.crd"]
    )


@pytest.mark.parametrize("restraint", ["backbone", "heavy", "all", "none"])
@pytest.mark.parametrize(
    "protocol",
    [
        BSS.Protocol.Minimisation,
        BSS.Protocol.Equilibration,
        BSS.Protocol.Production,
    ],
)
def test_restrain_atoms(system, restraint, protocol):
    process = BSS.Process.OpenMM(
        system, protocol(restraint=restraint, force_constant=1000)
    )

    expected_lines = [  # (with restraint, without restraint)
        (
            f"restrained_atoms = {system.getRestraintAtoms(restraint) if restraint != 'none' else []}",
            None,
        ),
        (
            "restraint.addBond(i, j, 0 * nanometers, 418400.0 * kilojoules_per_mole / nanometer**2)",
            None,
        ),
        (
            "simulation.context.setPositions(positions)",
            "simulation.context.setPositions(prm.positions)",
        ),
    ]

    with open(f"{process._work_dir}/{process._name}_script.py") as fp:
        content = fp.read()

    for with_restraint_line, without_restraint_line in expected_lines:
        if restraint == "none":
            if without_restraint_line is not None:
                assert without_restraint_line in content
            if with_restraint_line is not None:
                assert with_restraint_line not in content
        else:
            if with_restraint_line is not None:
                assert with_restraint_line in content
            if without_restraint_line is not None:
                assert without_restraint_line not in content


@pytest.mark.parametrize("restraint", ["backbone", "heavy", "all", "none"])
def test_minimise(system, restraint):
    """Test a minimisation protocol."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(
        steps=100, restraint=restraint, force_constant=1000
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol, restraint=restraint, tolerance=0.1)


@pytest.mark.parametrize("restraint", ["backbone", "heavy", "all", "none"])
def test_equilibrate(system, restraint):
    """Test an equilibration protocol."""

    # Create a short equilibration protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"),
        restraint=restraint,
        force_constant=1000,
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol, restraint=restraint, tolerance=0.2)


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


@pytest.mark.parametrize("restraint", ["backbone", "heavy", "all", "none"])
def test_production(system, restraint):
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(
        runtime=BSS.Types.Time(0.001, "nanoseconds"),
        restraint=restraint,
        force_constant=1000,
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol, restraint=restraint, tolerance=1.1)


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


def run_process(system, protocol, restraint="none", tolerance=0.05):
    """Helper function to run various simulation protocols.

    Args:
        system: A BSS test system instance
        protocol: A BSS protocol instance

    Returns:
        A BSS system

    Raises:
        AssertionError:
    """

    # Need to adjust protocol defaults to get trajectory frames
    # in short test runs.
    if not isinstance(protocol, BSS.Protocol.Minimisation):
        protocol.setReportInterval(100)
        protocol.setRestartInterval(100)

    # Initialise the OpenMM  process.
    process = BSS.Process.OpenMM(system, protocol, name="test")

    # Start the OpenMM simulation.
    process.start()

    # Wait for the process to end.
    process.wait()

    # Make sure the process didn't error.
    assert not process.isError()

    # Make sure that we get a molecular system back.
    new_system = process.getSystem()
    assert new_system is not None

    if restraint == "none":
        return

    # Check if restrained atoms stayed in place
    restraint_atoms = system.getRestraintAtoms(restraint)
    did_not_move = [
        get_dist_atoms(system.getAtom(x), new_system.getAtom(x)) < tolerance
        for x in restraint_atoms
    ]

    assert all(did_not_move)


def get_dist_atoms(a, b):
    """Return the distance between two `BioSimSpace.Atom` instances"""
    v = a.coordinates() - b.coordinates()
    return math.sqrt(
        sum(
            [
                v.x().angstroms().value() ** 2,
                v.y().angstroms().value() ** 2,
                v.z().angstroms().value() ** 2,
            ]
        )
    )
