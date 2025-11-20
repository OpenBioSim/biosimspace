import os
import pytest
import tempfile

import BioSimSpace as BSS
from BioSimSpace._SireWrappers import System, ReplicaSystem

from tests.Sandpit.Exscientia.conftest import url


@pytest.fixture(scope="module")
def system():
    """Solvated alanine dipeptide system."""
    return BSS.IO.readMolecules(["tests/input/ala.top", "tests/input/ala.crd"])


@pytest.fixture(scope="module")
def perturbable_system():
    """A vacuum perturbable system."""
    return BSS.IO.readPerturbableSystem(
        f"{url}/perturbable_system0.prm7",
        f"{url}/perturbable_system0.rst7",
        f"{url}/perturbable_system1.prm7",
        f"{url}/perturbable_system1.rst7",
    )


@pytest.fixture(scope="module")
def replica_system(system):
    """A replica system for testing."""
    return ReplicaSystem(system, num_replicas=10)


@pytest.fixture(scope="module")
def perturbable_replica_system(perturbable_system):
    """A perturbable replica system for testing."""
    return ReplicaSystem(perturbable_system, num_replicas=10)


@pytest.fixture(scope="module")
def squashed_perturbable_replica_system(perturbable_system):
    """A squashed perturbable replica system for testing."""
    return ReplicaSystem(perturbable_system, num_replicas=10, is_squashed=True)


@pytest.mark.parametrize("rs", ["replica_system", "perturbable_replica_system"])
def test_num_replicas(rs, request):
    """Test the number of replicas in the replica system."""
    replica_system = request.getfixturevalue(rs)
    assert replica_system.nReplicas() == 10


@pytest.mark.parametrize("rs", ["replica_system", "perturbable_replica_system"])
def test_get_replica(rs, request):
    """Test retrieving a replica system."""
    replica_system = request.getfixturevalue(rs)

    replica = replica_system.getReplica(0)
    assert isinstance(replica, System)

    # Make sure negative indexing works.
    replica_neg = replica_system.getReplica(-1)
    assert isinstance(replica_neg, System)

    # Make sure __getitem__ works.
    replica_item = replica_system[-1]

    # Make sure bounds checking works.
    with pytest.raises(IndexError):
        replica_system.getReplica(10)
    with pytest.raises(IndexError):
        replica_system.getReplica(-11)


@pytest.mark.parametrize("rs", ["replica_system", "perturbable_replica_system"])
def test_stream(rs, request):
    """Test streaming the replica system to a file."""
    replica_system = request.getfixturevalue(rs)

    with tempfile.TemporaryDirectory() as tmpdir:
        stream, trajectory = replica_system.save(f"{tmpdir}/replica_system")

        # Check that files were created.
        assert os.path.exists(stream)
        assert os.path.exists(trajectory)

        # Check that we can load the replica system back.
        rs = ReplicaSystem.load(stream, trajectory)

        # Check that the number of replicas matches.
        assert rs.nReplicas() == replica_system.nReplicas()


@pytest.mark.parametrize("rs", ["replica_system", "perturbable_replica_system"])
def test_save_load_replicas(rs, request):
    """Test saving/loading individual replicas to/from files."""
    replica_system = request.getfixturevalue(rs)

    with tempfile.TemporaryDirectory() as tmpdir:
        # Create the list of filenames, one for each replica.
        filenames = [
            f"{tmpdir}/replica_{i:03d}.gro" for i in range(replica_system.nReplicas())
        ]

        # Save the replicas.
        replica_system.saveReplicas(filenames)

        # Check that files were created.
        for filename in filenames:
            assert os.path.exists(filename)

        # Make sure an exception is raised if the number of filenames does not
        # match the number of replicas.
        with pytest.raises(ValueError):
            replica_system.saveReplicas(filenames[:-1])

        # Make sure an exception is raised the file names don't all have the same
        # extension.
        bad_filenames = filenames.copy()
        bad_filenames[0] = f"{tmpdir}/replica_000.rst7"
        with pytest.raises(IOError):
            replica_system.saveReplicas(bad_filenames)

        # Try to load the replicas back.
        rs = ReplicaSystem.loadReplicas(replica_system, filenames)

        # Check that the number of replicas matches.
        assert rs.nReplicas() == replica_system.nReplicas()

        # Only load back a subset of replicas.
        rs = ReplicaSystem.loadReplicas(replica_system, filenames[:5])

        # Check that the number of replicas matches.
        assert rs.nReplicas() == 5


def test_squashed_representation(squashed_perturbable_replica_system):
    """Test that the internal squashed trajectory representation works."""
    # Check that the number of atoms in the internal system is less
    # than the number of atoms in the internal trajectory.
    assert (
        squashed_perturbable_replica_system._new_sire_object.num_atoms()
        < squashed_perturbable_replica_system._trajectory[0].current().num_atoms()
    )

    # Check that we can retrieve a replica system.
    replica = squashed_perturbable_replica_system.getReplica(0)
    assert isinstance(replica, System)

    # Make sure that the replica has the correct number of atoms.
    assert (
        replica.nAtoms()
        == squashed_perturbable_replica_system._new_sire_object.num_atoms()
    )

    # Test saving and loading the squashed replica system.
    with tempfile.TemporaryDirectory(delete=False) as tmpdir:
        stream, trajectory = squashed_perturbable_replica_system.save(
            f"{tmpdir}/squashed_replica_system"
        )

        # Check that files were created.
        assert os.path.exists(stream)
        assert os.path.exists(trajectory)

        # Check that we can load the replica system back.
        rs = ReplicaSystem.load(stream, trajectory)

        # Check that the number of replicas matches.
        assert rs.nReplicas() == squashed_perturbable_replica_system.nReplicas()

        # Make sure that the trajectory is still squashed.
        assert rs._is_squashed is True
        assert rs._new_sire_object.num_atoms() < rs._trajectory[0].current().num_atoms()
