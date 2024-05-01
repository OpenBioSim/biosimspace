import math
import numpy as np
import shutil
import pytest

import BioSimSpace as BSS

from BioSimSpace.Units.Angle import radian
from BioSimSpace.Units.Energy import kcal_per_mol, kj_per_mol
from BioSimSpace.Units.Length import angstrom
from BioSimSpace.Units.Pressure import bar
from BioSimSpace.Units.Temperature import kelvin
from BioSimSpace.Units.Time import picosecond
from BioSimSpace.Units.Volume import nanometer3

from tests.conftest import url, has_amber, has_gromacs, has_openff

# Store the allowed restraints.
restraints = BSS.Protocol._position_restraint_mixin._PositionRestraintMixin.restraints()


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
@pytest.mark.parametrize("restraint", restraints)
def test_minimise(system, restraint):
    """Test a minimisation protocol."""

    # Create a short minimisation protocol.
    protocol = BSS.Protocol.Minimisation(steps=100, restraint=restraint)

    # Run the process, check that it finished without error, and returns a system.
    run_process(system, protocol)


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
@pytest.mark.parametrize("restraint", restraints)
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
@pytest.mark.parametrize("restraint", restraints)
def test_production(system, restraint):
    """Test a production protocol."""

    # Create a short production protocol.
    protocol = BSS.Protocol.Production(
        runtime=BSS.Types.Time(0.001, "nanoseconds"), restraint=restraint
    )

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
            "tests/output/gromacs.edr",
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


@pytest.mark.skipif(
    has_amber is False or has_gromacs is False or has_openff is False,
    reason="Requires AMBER, GROMACS, and OpenFF to be installed.",
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
