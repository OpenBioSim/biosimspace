import pytest
import socket

import BioSimSpace as BSS


@pytest.mark.skipif(
    socket.gethostname() != "porridge",
    reason="Local test requiring PLUMED patched GROMACS.",
)
def test_metadynamics(system):
    # Search for the first molecule containing ALA.
    molecule = system.search("resname ALA").molecules()[0]

    # Store the torsion indices.
    phi_idx = [4, 6, 8, 14]
    psi_idx = [6, 8, 14, 16]

    # Create the collective variables.
    phi = BSS.Metadynamics.CollectiveVariable.Torsion(atoms=phi_idx)
    psi = BSS.Metadynamics.CollectiveVariable.Torsion(atoms=psi_idx)

    # Create the metadynamics protocol.
    protocol = BSS.Protocol.Metadynamics(
        collective_variable=[phi, psi], runtime=100 * BSS.Units.Time.picosecond
    )

    # Run the metadynamics simulation.
    process = BSS.Metadynamics.run(molecule.toSystem(), protocol, gpu_support=True)

    # Wait for the process to finish.
    process.wait()

    # Check if the process has finished successfully.
    assert not process.isError()

    free_nrg = process.getFreeEnergy(kt=BSS.Units.Energy.kt)

    # Check if the free energy is not None.
    assert free_nrg is not None


@pytest.mark.skipif(
    socket.gethostname() != "porridge",
    reason="Local test requiring PLUMED patched GROMACS.",
)
def test_steering(system):
    # Create the collective variable. Here we align on molecule index 1
    # and all atoms not belonging to to the ALA residue in molecule index 0.
    # We compute the RSMD using the atoms belonging to the ALA residue.
    cv = BSS.Metadynamics.CollectiveVariable.RMSD(
        system, system, "(molidx 0 and not resname ALA) or molidx 1", "resname ALA"
    )

    # Add some stages.
    start = 0 * BSS.Units.Time.nanosecond
    apply_force = 4 * BSS.Units.Time.picosecond
    steer = 50 * BSS.Units.Time.picosecond
    relax = 100 * BSS.Units.Time.picosecond

    # Create some restraints.
    nm = BSS.Units.Length.nanometer
    restraint_1 = BSS.Metadynamics.Restraint(cv.getInitialValue(), 0)
    restraint_2 = BSS.Metadynamics.Restraint(cv.getInitialValue(), 3500)
    restraint_3 = BSS.Metadynamics.Restraint(0 * nm, 3500)
    restraint_4 = BSS.Metadynamics.Restraint(0 * nm, 0)

    # Create the steering protocol.
    protocol = BSS.Protocol.Steering(
        cv,
        [start, apply_force, steer, relax],
        [restraint_1, restraint_2, restraint_3, restraint_4],
        runtime=100 * BSS.Units.Time.picosecond,
    )

    # Create the steering process.
    process = BSS.Process.Gromacs(system, protocol, extra_args={"-ntmpi": 1})

    # Start the process and wait for it to finish.
    process.start()
    process.wait()

    # Check if the process has finished successfully.
    assert not process.isError()


@pytest.mark.skipif(
    socket.gethostname() != "porridge",
    reason="Local test requiring PLUMED.",
)
def test_funnel_metadynamics():
    # Load the protein-ligand system.
    system = BSS.IO.readMolecules(
        BSS.IO.expand(BSS.tutorialUrl(), ["funnel_system.rst7", "funnel_system.prm7"])
    )

    # Get the p0 and p1 points for defining the funnel.
    p0, p1 = BSS.Metadynamics.CollectiveVariable.makeFunnel(system)

    # Expected p0 and p1 points.
    expected_p0 = [1017, 1031, 1050, 1186, 1205, 1219, 1238, 2585, 2607, 2623]
    expected_p1 = [
        519,
        534,
        553,
        572,
        583,
        597,
        608,
        619,
        631,
        641,
        1238,
        1254,
        1280,
        1287,
        1306,
        1313,
        1454,
        1473,
        1480,
        1863,
        1879,
        1886,
        1906,
        2081,
        2116,
        2564,
        2571,
        2585,
        2607,
    ]

    # Make sure the p0 and p1 points are as expected.
    assert p0 == expected_p0
    assert p1 == expected_p1

    # Set the upper bound for the funnel collective variable.
    upper_bound = BSS.Metadynamics.Bound(value=3.5 * BSS.Units.Length.nanometer)

    # Create the funnel collective variable.
    cv = BSS.Metadynamics.CollectiveVariable.Funnel(p0, p1, upper_bound=upper_bound)

    # Create the metadynamics protocol.
    protocol = BSS.Protocol.Metadynamics(
        cv,
        runtime=1 * BSS.Units.Time.picosecond,
        hill_height=1.5 * BSS.Units.Energy.kj_per_mol,
        hill_frequency=500,
        restart_interval=1000,
        bias_factor=10,
    )

    # Create the metadynamics process.
    process = BSS.Process.OpenMM(system, protocol)

    # Start the process and wait for it to finish.
    process.start()
    process.wait()

    # Check if the process has finished successfully.
    assert not process.isError()
