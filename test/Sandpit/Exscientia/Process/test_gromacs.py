import numpy as np
import shutil
import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align import decouple
from BioSimSpace.Sandpit.Exscientia.FreeEnergy import Restraint
from BioSimSpace.Sandpit.Exscientia.Units.Angle import radian
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol, kj_per_mol
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
from BioSimSpace.Sandpit.Exscientia.Units.Pressure import bar
from BioSimSpace.Sandpit.Exscientia.Units.Temperature import kelvin
from BioSimSpace.Sandpit.Exscientia.Units.Time import picosecond
from BioSimSpace.Sandpit.Exscientia.Units.Volume import nanometer3

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


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_write_restraint(system, tmp_path):
    """Test if the restraint has been written in a way that could be processed
    correctly.
    """
    ligand = ligand = BSS.IO.readMolecules(
        [f"{url}/ligand01.prm7.bz2", f"{url}/ligand01.rst7.bz2"]
    ).getMolecule(0)
    decoupled_ligand = decouple(ligand)
    l1 = decoupled_ligand.getAtoms()[0]
    l2 = decoupled_ligand.getAtoms()[1]
    l3 = decoupled_ligand.getAtoms()[2]
    ligand_2 = BSS.IO.readMolecules(
        [f"{url}/ligand04.prm7.bz2", f"{url}/ligand04.rst7.bz2"]
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
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(0.0001, "nanoseconds"))

    # Run the process and check that it finishes without error.
    run_process(system, protocol, restraint=restraint, work_dir=str(tmp_path))
    with open(tmp_path / "test.top", "r") as f:
        assert "intermolecular_interactions" in f.read()


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


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
class TestGetRecord:
    @staticmethod
    @pytest.fixture()
    def setup(system):
        protocol = BSS.Protocol.Production(
            runtime=BSS.Types.Time(60, "picosecond"),
            timestep=BSS.Types.Time(4, "femtosecond"),
            report_interval=200,
        )
        process = BSS.Process.Gromacs(system, protocol)
        shutil.copyfile(
            "test/Sandpit/Exscientia/output/gromacs.edr",
            process.workDir() + "/gromacs.edr",
        )
        process._update_energy_dict()
        return process

    @pytest.mark.parametrize(
        "func,time_series,value,unit",
        [
            ("getStep", True, 0, 1),
            ("getStep", False, 15000, 1),
            ("getTime", True, 0.0, picosecond),
            ("getTime", False, 60.0, picosecond),
            ("getBondEnergy", True, 191.762558, kj_per_mol),
            ("getBondEnergy", False, 214.545654, kj_per_mol),
            ("getAngleEnergy", True, 908.064758, kj_per_mol),
            ("getAngleEnergy", False, 925.678772, kj_per_mol),
            ("getProperEnergy", True, 521.629150, kj_per_mol),
            ("getProperEnergy", False, 507.826538, kj_per_mol),
            ("getImproperEnergy", True, 0.000000, kj_per_mol),
            ("getImproperEnergy", False, 0.000000, kj_per_mol),
            ("getLennardJones14", True, 149.906967, kj_per_mol),
            ("getLennardJones14", False, 163.615982, kj_per_mol),
            ("getLennardJonesSR", True, 8568.126953, kj_per_mol),
            ("getLennardJonesSR", False, 8971.750000, kj_per_mol),
            ("getCoulomb14", True, -11420.338867, kj_per_mol),
            ("getCoulomb14", False, -11493.004883, kj_per_mol),
            ("getCoulombSR", True, -63410.824219, kj_per_mol),
            ("getCoulombSR", False, -64019.128906, kj_per_mol),
            ("getCoulombReciprocal", True, 586.427307, kj_per_mol),
            ("getCoulombReciprocal", False, 251.100098, kj_per_mol),
            ("getDispersionCorrection", True, -575.343933, kj_per_mol),
            ("getDispersionCorrection", False, -577.368042, kj_per_mol),
            ("getPotentialEnergy", True, -64480.589844, kj_per_mol),
            ("getPotentialEnergy", False, -65054.984375, kj_per_mol),
            ("getKineticEnergy", True, 11556.777344, kj_per_mol),
            ("getKineticEnergy", False, 11280.177734, kj_per_mol),
            ("getTotalEnergy", True, -52923.812500, kj_per_mol),
            ("getTotalEnergy", False, -53774.804688, kj_per_mol),
            ("getConservedEnergy", True, -52937.847656, kj_per_mol),
            ("getConservedEnergy", False, -51628.320312, kj_per_mol),
            ("getTemperature", True, 306.766907, kelvin),
            ("getTemperature", False, 299.424744, kelvin),
            ("getPressure", True, 119.490417, bar),
            ("getPressure", False, 756.571045, bar),
            ("getPressureDC", True, -214.083145, bar),
            ("getPressureDC", False, -215.590363, bar),
            ("getVolume", True, 44.679958, nanometer3),
            ("getVolume", False, 44.523510, nanometer3),
        ],
    )
    def test_get(self, setup, func, time_series, value, unit):
        energy = getattr(setup, func)(time_series, block=False)
        if time_series:
            assert len(energy) == 76
            np.testing.assert_almost_equal(energy[0] / unit, value, decimal=3)
        else:
            np.testing.assert_almost_equal(energy / unit, value, decimal=3)
