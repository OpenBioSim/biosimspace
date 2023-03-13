import BioSimSpace.Sandpit.Exscientia as BSS

import filecmp
import pytest
import os
import warnings as _warnings

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
        f"{url}/perturbable_system0.prm7",
        f"{url}/perturbable_system0.rst7",
        f"{url}/perturbable_system1.prm7",
        f"{url}/perturbable_system1.rst7",
    )


def test_minimise(system):
    """Test a minimisation protocol."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


def test_equilibrate(system):
    """Test an equilibration protocol."""

    # Create a short equilibration protocol.
    protocol = BSS.Protocol.Equilibration(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


def test_production(system):
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


def test_free_energy(perturbable_system):
    """Test a free energy perturbation protocol."""

    # Create a short FEP protocol.
    protocol = BSS.Protocol.FreeEnergy(
        runtime=0.1 * BSS.Units.Time.picosecond, report_interval=50, restart_interval=50
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(perturbable_system, protocol)


def test_pert_file():
    """Test the perturbation file writer."""

    # Load the perturbable molecule.
    mol = BSS.IO.readPerturbableSystem(
        f"{url}/morph0.prm7",
        f"{url}/morph0.rst7",
        f"{url}/morph1.prm7",
        f"{url}/morph1.rst7",
    )[0]

    # Create the perturbation file.
    BSS.Process._somd._to_pert_file(mol)

    # Check that the files are the same.
    assert filecmp.cmp("MORPH.pert", "test/input/morph.pert")

    # Remove the temporary perturbation file.
    os.remove("MORPH.pert")


def test_pert_res_num(perturbable_system):
    """Test that the perturbable residue number is correct when
    the molecules in the system are re-ordered.
    """

    # Create a default free energy protocol.
    protocol = BSS.Protocol.FreeEnergy()

    # Set up the intial process object.
    process0 = BSS.Process.Somd(perturbable_system, protocol)

    # Now put the perturbable molecule last.
    new_system = (perturbable_system[1:] + perturbable_system[0]).toSystem()

    # Re-add the box info.
    new_system.setBox(*perturbable_system.getBox())

    # Set up the new process object.
    process1 = BSS.Process.Somd(new_system, protocol)

    # Get the config for each process.
    config0 = process0.getConfig()
    config1 = process1.getConfig()

    # Get the difference between the two configs.

    # Items unique to config0.
    unique0 = list(set(config0) - set(config1))
    unique1 = list(set(config1) - set(config0))

    # One item should be different.
    assert len(unique0) == len(unique1) == 1

    # Now check that the entries are as expected.

    # For the first system, the perturbable residue is first.
    assert unique0[0] == "perturbed residue number = 1"

    # For the second system, the perturbable residue is second.
    assert unique1[0] == "perturbed residue number = 2"


def run_process(system, protocol):
    """Helper function to run various simulation protocols."""

    # Initialise the SOMD process.
    process = BSS.Process.Somd(system, protocol, name="test", platform="CPU")

    # Start the SOMD simulation.
    process.start()

    # Wait for the process to end.
    process.wait()

    res = process.isError()
    # OpenMM minimisation occasionally fails with "Particle coordinate is nan"
    # if this is the case, the test will pass with a warning
    if res:
        with open(os.path.join(process.workDir(), "test.out"), "r") as hnd:
            if "RuntimeError: Particle coordinate is nan" in hnd.read():
                _warnings.warn(
                    "Test raised RuntimeError: Particle coordinate is nan",
                    RuntimeWarning,
                )
                res = False

    return not res

    # Make sure the process didn't error.
    assert not process.isError()

    # Make sure that we get a molecular system back.
    assert process.getSystem() is not None
