import BioSimSpace as BSS

import pytest

# Make sure GROMSCS is installed.
has_gromacs = BSS._gmx_exe is not None

# Store the tutorial URL.
url = BSS.tutorialUrl()


@pytest.fixture(scope="session")
def system():
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules(["test/input/ala.top", "test/input/ala.crd"])


@pytest.fixture(scope="session")
def perturbable_system():
    """Re-use the same perturbable system for each test."""
    return BSS.IO.readPerturbableSystem(
        f"{url}/complex_vac0.prm7.bz2",
        f"{url}/complex_vac0.rst7.bz2",
        f"{url}/complex_vac1.prm7.bz2",
        f"{url}/complex_vac1.rst7.bz2",
    )


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_minimise(system):
    """Test a minimisation protocol."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
@pytest.mark.parametrize("restraint", ["backbone", "heavy", "all", "none"])
def test_equilibrate(system, restraint):
    """Test an equilibration protocol."""

    # Create a short equilibration protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"), restraint=restraint
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
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


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_cool(system):
    """Test a cooling protocol."""

    # Create a short heating protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"),
        temperature_start=BSS.Types.Temperature(300, "kelvin"),
        temperature_end=BSS.Types.Temperature(0, "kelvin"),
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol, extra_options={"verlet-buffer-tolerance": "2e-07"})


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_production(system):
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_vacuum_water(system):
    """Regression test for ensuring the water topology is swapped for vacuum simulations."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Create a new system using the first two molecules of 'system'.
    # This will be an alanine-dipeptide and water in vacuum.
    new_system = (system[0] + system[1]).toSystem()

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
@pytest.mark.parametrize("restraint", ["backbone", "heavy"])
def test_restraints(perturbable_system, restraint):
    """Regression test for correct injection of restraint file into GROMACS topology."""

    # Create an equilibration protocol with backbone restraints.
    protocol = BSS.Protocol.Equilibration(restraint=restraint)

    # Create the simulation process.
    process = BSS.Process.Gromacs(perturbable_system, protocol)


def run_process(system, protocol, **kwargs):
    """Helper function to run various simulation protocols."""

    # Initialise the GROMACS process.
    process = BSS.Process.Gromacs(system, protocol, name="test", **kwargs)

    # Only run on a single MPI rank.
    process.setArg("-ntmpi", 1)

    # Start the GROMACS simulation.
    process.start()

    # Wait for the process to end.
    process.wait()

    # Make sure the process didn't error.
    assert not process.isError()

    # Make sure that we get a molecular system back.
    assert process.getSystem() is not None
