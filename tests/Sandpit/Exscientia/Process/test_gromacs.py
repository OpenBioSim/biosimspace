import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
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
from tests.Sandpit.Exscientia.conftest import (
    has_alchemtest,
    has_amber,
    has_gromacs,
    has_openff,
    has_pyarrow,
    url,
)
from tests.conftest import root_fp


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
        f"{url}/complex_vac0.prm7.bz2",
        f"{url}/complex_vac0.rst7.bz2",
        f"{url}/complex_vac1.prm7.bz2",
        f"{url}/complex_vac1.rst7.bz2",
    )


@pytest.mark.skipif(
    has_gromacs is False or has_pyarrow is False,
    reason="Requires GROMACS and pyarrow to be installed.",
)
def test_minimise(system):
    """Test a minimisation protocol."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.skipif(
    has_gromacs is False or has_pyarrow is False,
    reason="Requires GROMACS and pyarrow to be installed.",
)
@pytest.mark.parametrize("restraint", ["backbone", "heavy", "all", "none"])
def test_equilibrate(system, restraint):
    """Test an equilibration protocol."""

    # Create a short equilibration protocol.
    protocol = BSS.Protocol.Equilibration(
        runtime=BSS.Types.Time(0.001, "nanoseconds"), restraint=restraint
    )

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.skipif(
    has_gromacs is False or has_pyarrow is False,
    reason="Requires GROMACS and pyarrow to be installed.",
)
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


@pytest.mark.skipif(
    has_gromacs is False or has_pyarrow is False,
    reason="Requires GROMACS and pyarrow to be installed.",
)
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


@pytest.mark.skipif(
    has_gromacs is False or has_pyarrow is False,
    reason="Requires GROMACS and pyarrow to be installed.",
)
def test_production(system):
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(0.001, "nanoseconds"))

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.skipif(
    has_gromacs is False or has_pyarrow is False,
    reason="Requires GROMACS and pyarrow to be installed.",
)
def test_vacuum_water(system):
    """Regression test for ensuring the water topology is swapped for vacuum simulations."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Create a new system using the first two molecules of 'system'.
    # This will be an alanine-dipeptide and water in vacuum.
    new_system = (system[0] + system[1]).toSystem()

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.skipif(
    has_gromacs is False or has_pyarrow is False,
    reason="Requires GROMACS and pyarrow to be installed.",
)
@pytest.mark.parametrize("restraint", ["backbone", "heavy"])
def test_restraints(perturbable_system, restraint):
    """Regression test for correct injection of restraint file into GROMACS topology."""

    # Create an equilibration protocol with backbone restraints.
    protocol = BSS.Protocol.Equilibration(restraint=restraint)

    # Create the simulation process.
    process = BSS.Process.Gromacs(perturbable_system, protocol)


@pytest.mark.skipif(
    has_gromacs is False or has_pyarrow is False,
    reason="Requires GROMACS and pyarrow to be installed.",
)
class TestRestraint:
    @pytest.fixture(scope="class")
    def setup(self):
        ligand = BSS.IO.readMolecules(
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
            "anchor_points": {
                "r1": r1,
                "r2": r2,
                "r3": r3,
                "l1": l1,
                "l2": l2,
                "l3": l3,
            },
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
        return system, restraint

    def test_regular_protocol(self, setup, tmp_path_factory):
        """Test if the restraint has been written in a way that could be processed
        correctly.
        """
        tmp_path = tmp_path_factory.mktemp("out")
        system, restraint = setup
        # Create a short production protocol.
        protocol = BSS.Protocol.FreeEnergy(
            runtime=BSS.Types.Time(0.0001, "nanoseconds"), perturbation_type="full"
        )

        # Run the process and check that it finishes without error.
        run_process(system, protocol, restraint=restraint, work_dir=str(tmp_path))
        with open(tmp_path / "test.top", "r") as f:
            assert "intermolecular_interactions" in f.read()

    def test_restraint_lambda(self, setup, tmp_path_factory):
        """Test if the restraint has been written correctly when restraint lambda is evoked."""
        tmp_path = tmp_path_factory.mktemp("out")
        system, restraint = setup
        # Create a short production protocol.
        protocol = BSS.Protocol.FreeEnergy(
            runtime=BSS.Types.Time(0.0001, "nanoseconds"),
            lam=pd.Series(data={"bonded": 0.0, "restraint": 0.0}),
            lam_vals=pd.DataFrame(data={"bonded": [0.0, 1.0], "restraint": [0.0, 1.0]}),
        )

        # Run the process and check that it finishes without error.
        run_process(system, protocol, restraint=restraint, work_dir=str(tmp_path))
        with open(tmp_path / "test.top", "r") as f:
            assert "dihedral_restraints" in f.read()
        with open(tmp_path / "test.mdp", "r") as f:
            assert "restraint-lambdas" in f.read()


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


@pytest.mark.skipif(
    has_gromacs is False or has_pyarrow is False or has_alchemtest is False,
    reason="Requires GROMACS, alchemtest, and pyarrow to be installed.",
)
class TestGetRecord:
    @staticmethod
    @pytest.fixture()
    def setup(perturbable_system):
        from alchemtest.gmx import load_ABFE

        protocol = BSS.Protocol.FreeEnergy(
            runtime=BSS.Types.Time(60, "picosecond"),
            timestep=BSS.Types.Time(4, "femtosecond"),
            report_interval=200,
        )
        process = BSS.Process.Gromacs(perturbable_system, protocol)
        shutil.copyfile(
            f"{root_fp}/Sandpit/Exscientia/output/gromacs.edr",
            process.workDir() + "/gromacs.edr",
        )
        shutil.copyfile(
            load_ABFE().data["ligand"][0],
            process.workDir() + "/gromacs.xvg",
        )
        process.saveMetric()
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

    def test_metric_parquet_exist(self, setup):
        assert Path(f"{setup.workDir()}/metric.parquet").exists()

    def test_metric_parquet(self, setup):
        df = pd.read_parquet(f"{setup.workDir()}/metric.parquet")
        assert np.isclose(df["PotentialEnergy (kJ/mol)"][0.0], -64480.589844)
        assert np.isclose(df["Volume (nm^3)"][0.0], 44.679958)
        assert np.isclose(df["Pressure (bar)"][0.0], 119.490417)
        assert np.isclose(df["Temperature (kelvin)"][0.0], 306.766907)

    def test_dhdl_parquet_exist(self, setup):
        assert Path(f"{setup.workDir()}/dHdl.parquet").exists()

    def test_dhdl_parquet(self, setup):
        df = pd.read_parquet(f"{setup.workDir()}/dHdl.parquet")
        assert df.shape == (1001, 2)

    def test_u_nk_parquet_exist(self, setup):
        assert Path(f"{setup.workDir()}/u_nk.parquet").exists()

    def test_u_nk_parquet(self, setup):
        df = pd.read_parquet(f"{setup.workDir()}/u_nk.parquet")
        assert df.shape == (1001, 20)

    def test_error_alchemlyb_extract(self, perturbable_system, monkeypatch):
        def extract(*args):
            raise ValueError("alchemlyb.parsing.gmx.extract failed.")

        monkeypatch.setattr("alchemlyb.parsing.gmx.extract", extract)
        # Create a process using any system and the protocol.
        process = BSS.Process.Gromacs(
            perturbable_system,
            BSS.Protocol.FreeEnergy(temperature=298 * BSS.Units.Temperature.kelvin),
        )
        process.saveMetric()
        with open(process.workDir() + "/gromacs.err", "r") as f:
            text = f.read()
            assert "Exception Information" in text

    @staticmethod
    @pytest.mark.parametrize(
        ("filename", "u_nk", "dHdl"),
        [
            ("metric.parquet", None, None),
            (None, "u_nk.parquet", None),
            (None, None, "dHdl.parquet"),
        ],
    )
    def test_selective(perturbable_system, filename, u_nk, dHdl):
        from alchemtest.gmx import load_ABFE

        protocol = BSS.Protocol.FreeEnergy(
            runtime=BSS.Types.Time(60, "picosecond"),
            timestep=BSS.Types.Time(4, "femtosecond"),
            report_interval=200,
        )
        process = BSS.Process.Gromacs(perturbable_system, protocol)
        shutil.copyfile(
            f"{root_fp}/Sandpit/Exscientia/output/gromacs.edr",
            process.workDir() + "/gromacs.edr",
        )
        shutil.copyfile(
            load_ABFE().data["ligand"][0],
            process.workDir() + "/gromacs.xvg",
        )

        process.saveMetric(filename, u_nk, dHdl)
        if filename is not None:
            assert (Path(process.workDir()) / filename).exists()
        else:
            assert not (Path(process.workDir()) / "metric.parquet").exists()
        if u_nk is not None:
            assert (Path(process.workDir()) / u_nk).exists()
        else:
            assert not (Path(process.workDir()) / "u_nk.parquet").exists()
        if dHdl is not None:
            assert (Path(process.workDir()) / dHdl).exists()
        else:
            assert not (Path(process.workDir()) / "dHdl.parquet").exists()


def test_error_saveMetric(perturbable_system):
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.Gromacs(perturbable_system, protocol)
    process.saveMetric()
    with open(f"{process.workDir()}/{process._name}.err", "r") as f:
        assert "Exception Information during saveMetric():" in f.read()


@pytest.mark.skipif(
    has_amber is False
    or has_gromacs is False
    or has_openff is False
    or has_pyarrow is False,
    reason="Requires AMBER, GROMACS, OpenFF, and pyarrow to be installed.",
)
def test_vacuum_com():
    """
    Test to ensure that the center of mass is moved to the origin following
    vacuum simulation. Because of the need to fake vacuum simulations by
    using an extremeley large box, molecular coordinates can exceed the
    format limit used by certain molecular file formats, e.g. AMBER RST7.
    As such, we translate the center of mass of the system to the origin.
    """

    # Create a test molecule.
    mol = BSS.Parameters.openff_unconstrained_2_0_0("CC").getMolecule()

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100)

    # Create a process object.
    process = BSS.Process.Gromacs(mol.toSystem(), protocol)

    # Start the process and wait for it to finish.
    process.start()
    process.wait()

    # Make sure it worked.
    assert not process.isError()

    # Get the minimised system.
    system = process.getSystem()

    # Make sure it worked.
    assert system is not None

    # Get the center of mass.
    com = system._getCenterOfMass()

    # Make sure each component is close to zero.
    for x in com:
        assert math.isclose(x.value(), 0, abs_tol=1e-6)
