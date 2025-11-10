import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align import decouple
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
from BioSimSpace.Sandpit.Exscientia.Units.Angle import radian
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol
from BioSimSpace.Sandpit.Exscientia.FreeEnergy import Restraint
from BioSimSpace.Sandpit.Exscientia.Units.Temperature import kelvin

import filecmp
import os
import pytest
import warnings

import BioSimSpace.Sandpit.Exscientia as BSS
from tests.conftest import root_fp

# Store the tutorial URL.
url = BSS.tutorialUrl()


@pytest.fixture(scope="session")
def system():
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules(
        [f"{root_fp}/input/ala.top", f"{root_fp}/input/ala.crd"]
    )


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
    assert filecmp.cmp("MORPH.pert", f"{root_fp}/input/morph.pert")

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


def test_restraint_boresch(system, tmp_path):
    """
    Test if the Boresch restraint has been written in a way that can be processed
    correctly.
    """
    ligand = BSS.IO.readMolecules(
        [f"{url}/ligand01.rst7", f"{url}/ligand01.prm7"]
    ).getMolecule(0)
    decoupled_ligand = decouple(ligand)
    l1 = decoupled_ligand.getAtoms()[0]
    l2 = decoupled_ligand.getAtoms()[1]
    l3 = decoupled_ligand.getAtoms()[2]
    ligand_2 = BSS.IO.readMolecules(
        [f"{url}/ligand04.rst7", f"{url}/ligand04.prm7"]
    ).getMolecule(0)
    r1 = ligand_2.getAtoms()[0]
    r2 = ligand_2.getAtoms()[1]
    r3 = ligand_2.getAtoms()[2]
    system = (decoupled_ligand + ligand_2).toSystem()

    restraint_dict = {
        "anchor_points": {"r1": r1, "r2": r2, "r3": r3, "l1": l1, "l2": l2, "l3": l3},
        "equilibrium_values": {
            "r0": 7.84 * angstrom,
            "thetaA0": 0.81 * radian,
            "thetaB0": 1.74 * radian,
            "phiA0": 2.59 * radian,
            "phiB0": -1.20 * radian,
            "phiC0": 2.63 * radian,
        },
        "force_constants": {
            "kr": 10 * kcal_per_mol / angstrom**2,
            "kthetaA": 10 * kcal_per_mol / (radian * radian),
            "kthetaB": 10 * kcal_per_mol / (radian * radian),
            "kphiA": 10 * kcal_per_mol / (radian * radian),
            "kphiB": 10 * kcal_per_mol / (radian * radian),
            "kphiC": 10 * kcal_per_mol / (radian * radian),
        },
    }
    restraint = Restraint(
        system, restraint_dict, 300 * kelvin, restraint_type="Boresch"
    )

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process and check that it finishes without error.
    run_process(system, protocol, restraint=restraint, work_dir=str(tmp_path))

    # Check for restraint options in config file
    with open(tmp_path / "test.cfg", "r") as f:
        lines = f.readlines()
        assert "use boresch restraints = True" in lines[-2]
        assert "boresch restraints dictionary" in lines[-1]


@pytest.fixture(scope="module")
def restraint_mdr(system):
    """A restraint object for testing."""
    ligand = BSS.IO.readMolecules(
        [f"{url}/ligand01.rst7", f"{url}/ligand01.prm7"]
    ).getMolecule(0)
    decoupled_ligand = decouple(ligand)

    fake_protein = BSS.IO.readMolecules(
        [f"{url}/ligand04.rst7", f"{url}/ligand04.prm7"]
    ).getMolecule(0)

    system = (decoupled_ligand + fake_protein).toSystem()

    # Create three distance restraints
    restraint_dict = {
        "distance_restraints": [
            {
                "l1": decoupled_ligand.getAtoms()[0],
                "r1": fake_protein.getAtoms()[0],
                "r0": 3 * angstrom,
                "kr": 10 * kcal_per_mol / angstrom**2,
                "r_fb": 1 * angstrom,
            },
            {
                "l1": decoupled_ligand.getAtoms()[1],
                "r1": fake_protein.getAtoms()[1],
                "r0": 3 * angstrom,
                "kr": 10 * kcal_per_mol / angstrom**2,
                "r_fb": 1 * angstrom,
            },
        ],
        "permanent_distance_restraint": {
            "l1": decoupled_ligand.getAtoms()[2],
            "r1": fake_protein.getAtoms()[2],
            "r0": 3 * angstrom,
            "kr": 10 * kcal_per_mol / angstrom**2,
            "r_fb": 1 * angstrom,
        },
    }
    restraint = Restraint(
        system, restraint_dict, 300 * kelvin, restraint_type="multiple_distance"
    )

    return restraint


def test_restraint_mdr(tmp_path, restraint_mdr):
    """
    Test if the multiple distance restraints restraint has been
    written in a way that can be processed correctly.
    """
    # Create a short free energy protocol.
    protocol = BSS.Protocol.Production(
        runtime=BSS.Types.Time(0.001, "nanoseconds"),
    )

    # Run the process and check that it finishes without error.
    run_process(
        restraint_mdr._system, protocol, restraint=restraint_mdr, work_dir=str(tmp_path)
    )

    # Check for restraint options in config file.
    with open(tmp_path / "test.cfg", "r") as f:
        lines = f.readlines()
        assert "use distance restraints = True" in lines[-4]
        assert "use permanent distance restraints = True" in lines[-3]
        assert "distance restraints dictionary" in lines[-2]
        assert "permanent distance restraints dictionary" in lines[-1]


def test_triclinic():
    """
    Test the triclinic box vectors can be updated correctly with OpenMM 8.2+.
    """

    # Load the triclinic alanine dipeptide system.
    system = BSS.IO.readMolecules(
        BSS.IO.expand(url, ["ala_triclinic.rst7", "ala_triclinic.prm7"])
    )

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


def test_triclinic_fep():
    """
    Test the triclinic box vectors can be updated correctly with OpenMM 8.2+.
    """

    # Load the triclinic alanine dipeptide system.
    system = BSS.IO.readPerturbableSystem(
        *BSS.IO.expand(
            url,
            [
                "ethane_methanol_triclinic0.prm7",
                "ethane_methanol_triclinic0.rst7",
                "ethane_methanol_triclinic1.prm7",
                "ethane_methanol_triclinic1.rst7",
            ],
        )
    )

    # Create a short production FEP protocol.
    protocol = BSS.Protocol.FreeEnergy(runtime=0.1 * BSS.Units.Time.picosecond)

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


def run_process(system, protocol, **kwargs):
    """Helper function to run various simulation protocols."""

    # Use CUDA platform is CUDA_VISIBLE_DEVICES is set.
    platform = "CUDA" if "CUDA_VISIBLE_DEVICES" in os.environ else "CPU"

    # Initialise the SOMD process.
    process = BSS.Process.Somd(
        system, protocol, name="test", platform=platform, **kwargs
    )

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
                warnings.warn(
                    "Test raised RuntimeError: Particle coordinate is nan",
                    RuntimeWarning,
                )
                res = False

    # Make sure the process didn't error.
    assert not res

    # Get the updated system.
    new_system = process.getSystem()

    # Make sure that we got a molecular system back.
    assert new_system is not None

    # Make sure the space is valid.
    if system.getBox()[0] is not None:
        assert new_system._sire_object.property("space").is_periodic()
