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
