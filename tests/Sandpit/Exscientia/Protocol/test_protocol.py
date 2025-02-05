import pytest

import BioSimSpace as BSS

# Unit tests for equivalence of protocol settings when instantiated
# using strings of unit-based types.


@pytest.fixture(scope="session")
def system():
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules(["tests/input/ala.top", "tests/input/ala.crd"])


def test_equilibration():
    # Instantiate from types.
    p0 = BSS.Protocol.Equilibration(
        timestep=1 * BSS.Units.Time.femtosecond,
        runtime=100 * BSS.Units.Time.picosecond,
        temperature_start=0 * BSS.Units.Temperature.kelvin,
        temperature_end=298 * BSS.Units.Temperature.kelvin,
        pressure=1 * BSS.Units.Pressure.atm,
        restraint="backbone",
        force_constant=100
        * BSS.Units.Energy.kj_per_mol
        / BSS.Units.Length.nanometer**2,
    )

    # Instantiate from strings.
    p1 = BSS.Protocol.Equilibration(
        timestep="1fs",
        runtime="100ps",
        temperature_start="0K",
        temperature_end="298K",
        pressure="1atm",
        restraint="backbone",
        force_constant="100 kJ / (mol nm^2)",
    )

    # Make sure the protocols are equivalent.
    assert p0 == p1


def test_production():
    # Instantiate from types.
    p0 = BSS.Protocol.Production(
        timestep=1 * BSS.Units.Time.femtosecond,
        runtime=100 * BSS.Units.Time.picosecond,
        temperature=298 * BSS.Units.Temperature.kelvin,
        pressure=1 * BSS.Units.Pressure.atm,
        restraint="backbone",
        force_constant=100
        * BSS.Units.Energy.kj_per_mol
        / BSS.Units.Length.nanometer**2,
    )

    # Instantiate from strings.
    p1 = BSS.Protocol.Production(
        timestep="1fs",
        runtime="100ps",
        temperature="298K",
        pressure="1atm",
        restraint="backbone",
        force_constant="100 kJ / (mol nm^2)",
    )

    # Make sure the protocols are equivalent.
    assert p0 == p1


def test_free_energy():
    # Instantiate from types.
    p0 = BSS.Protocol.FreeEnergy(
        timestep=1 * BSS.Units.Time.femtosecond,
        runtime=100 * BSS.Units.Time.picosecond,
        temperature=298 * BSS.Units.Temperature.kelvin,
        pressure=1 * BSS.Units.Pressure.atm,
        restraint="backbone",
        force_constant=100
        * BSS.Units.Energy.kj_per_mol
        / BSS.Units.Length.nanometer**2,
    )

    # Instantiate from strings.
    p1 = BSS.Protocol.FreeEnergy(
        timestep="1fs",
        runtime="100ps",
        temperature="298K",
        pressure="1atm",
        restraint="backbone",
        force_constant="100 kJ / (mol nm^2)",
    )

    # Make sure the protocols are equivalent.
    assert p0 == p1


def test_metadynamics():
    # Create a collective variable so we can instantiate a metadynamics protocol.
    cv = BSS.Metadynamics.CollectiveVariable.Torsion(atoms=[0, 1, 2, 3])

    # Instantiate from types.
    p0 = BSS.Protocol.Metadynamics(
        collective_variable=cv,
        timestep=1 * BSS.Units.Time.femtosecond,
        runtime=100 * BSS.Units.Time.picosecond,
        temperature=298 * BSS.Units.Temperature.kelvin,
        pressure=1 * BSS.Units.Pressure.atm,
    )

    # Instantiate from strings.
    p1 = BSS.Protocol.Metadynamics(
        collective_variable=cv,
        timestep="1fs",
        runtime="100ps",
        temperature="298K",
        pressure="1atm",
    )

    # Make sure the protocols are equivalent.
    assert p0 == p1


def test_steering(system):
    # Create a reference containing the first and third molecule from the system.
    reference = (system[0] + system[2]).toSystem()

    # Create the collective variable. Here we align on molecule index 1
    # and all atoms not belonging to to the ALA residue in molecule index 0.
    # We compute the RSMD using the atoms belonging to the ALA residue.
    cv = BSS.Metadynamics.CollectiveVariable.RMSD(
        system,
        reference,
        "(molidx 0 and not resname ALA) or molidx 1",
        "resname ALA",
        reference_mapping={0: 0, 1: 2},
    )

    # Create a schedule.
    start = 0 * BSS.Units.Time.nanosecond
    apply_force = 4 * BSS.Units.Time.picosecond
    steer = 150 * BSS.Units.Time.nanosecond
    relax = 152 * BSS.Units.Time.nanosecond

    # Create restraints.
    nm = BSS.Units.Length.nanometer
    restraint0 = BSS.Metadynamics.Restraint(cv.getInitialValue(), 0)
    restraint1 = BSS.Metadynamics.Restraint(cv.getInitialValue(), 3500)
    restraint2 = BSS.Metadynamics.Restraint(0 * nm, 3500)
    restraint3 = BSS.Metadynamics.Restraint(0 * nm, 0)

    # Instantiate from types.
    p0 = BSS.Protocol.Steering(
        cv,
        [start, apply_force, steer, relax],
        [restraint0, restraint1, restraint2, restraint3],
        timestep=2 * BSS.Units.Time.femtosecond,
        runtime=200 * BSS.Units.Time.nanosecond,
        temperature=298 * BSS.Units.Temperature.kelvin,
        pressure=1 * BSS.Units.Pressure.atm,
    )

    # Instantiate from strings.
    p1 = BSS.Protocol.Steering(
        cv,
        [start, apply_force, steer, relax],
        [restraint0, restraint1, restraint2, restraint3],
        timestep="2fs",
        runtime="200ns",
        temperature="298K",
        pressure="1atm",
    )

    # Make sure the protocols are equivalent.
    assert p0 == p1
