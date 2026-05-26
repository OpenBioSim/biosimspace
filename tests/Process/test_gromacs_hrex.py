"""Tests for BioSimSpace.Process.GromacsHREX."""

import os
import textwrap

import pytest

import BioSimSpace as BSS
from tests.conftest import has_gromacs

# ---------------------------------------------------------------------------
# Input validation tests: no GROMACS required; errors are raised before the
# exe lookup in Gromacs.__init__().
# ---------------------------------------------------------------------------


def test_non_fep_protocol(perturbable_system):
    """Non-FEP protocol raises TypeError."""
    with pytest.raises(TypeError, match="_FreeEnergyMixin"):
        BSS.Process.GromacsHREX(
            perturbable_system, BSS.Protocol.Minimisation(steps=100)
        )


def test_invalid_use_mpi(perturbable_system):
    """Non-bool use_mpi raises TypeError."""
    protocol = BSS.Protocol.FreeEnergy()
    with pytest.raises(TypeError, match="use_mpi"):
        BSS.Process.GromacsHREX(perturbable_system, protocol, use_mpi="yes")


def test_invalid_repex_frequency_zero(perturbable_system):
    """repex_frequency=0 raises ValueError."""
    protocol = BSS.Protocol.FreeEnergy()
    with pytest.raises(ValueError, match="repex_frequency"):
        BSS.Process.GromacsHREX(perturbable_system, protocol, repex_frequency=0)


def test_invalid_repex_frequency_float(perturbable_system):
    """Float repex_frequency raises ValueError."""
    protocol = BSS.Protocol.FreeEnergy()
    with pytest.raises(ValueError, match="repex_frequency"):
        BSS.Process.GromacsHREX(perturbable_system, protocol, repex_frequency=1.5)


def test_replica_count_mismatch(perturbable_system):
    """ReplicaSystem with wrong replica count raises ValueError."""
    from BioSimSpace._SireWrappers import ReplicaSystem

    protocol = BSS.Protocol.FreeEnergy()
    n_lam = len(protocol.getLambdaValues())
    wrong = ReplicaSystem(perturbable_system, num_replicas=n_lam + 2)
    with pytest.raises(ValueError, match="number of replicas"):
        BSS.Process.GromacsHREX(wrong, protocol)


def test_system_type_error():
    """Passing a non-System, non-ReplicaSystem raises TypeError."""
    protocol = BSS.Protocol.FreeEnergy()
    with pytest.raises(TypeError, match="'system' must be"):
        BSS.Process.GromacsHREX("not-a-system", protocol)


# ---------------------------------------------------------------------------
# Tests that require GROMACS to be installed.
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed.")
def test_system_promotion(perturbable_system):
    """Plain System input is promoted to a ReplicaSystem with the right replica count."""
    from BioSimSpace._SireWrappers import ReplicaSystem

    protocol = BSS.Protocol.FreeEnergy()
    n_lam = len(protocol.getLambdaValues())

    process = BSS.Process.GromacsHREX(perturbable_system, protocol)

    assert isinstance(process._replica_system, ReplicaSystem)
    assert process._replica_system.nReplicas() == n_lam


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed.")
def test_replica_system_passthrough(perturbable_system):
    """A ReplicaSystem with the correct replica count is accepted without re-wrapping."""
    from BioSimSpace._SireWrappers import ReplicaSystem

    protocol = BSS.Protocol.FreeEnergy()
    n_lam = len(protocol.getLambdaValues())
    rs = ReplicaSystem(perturbable_system, num_replicas=n_lam)

    process = BSS.Process.GromacsHREX(rs, protocol)

    assert isinstance(process._replica_system, ReplicaSystem)
    assert process._replica_system.nReplicas() == n_lam


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed.")
def test_file_layout(perturbable_system, tmp_path):
    """After setup, the shared topology and all per-lambda input files must exist."""
    protocol = BSS.Protocol.FreeEnergy()
    lam_vals = protocol.getLambdaValues()

    process = BSS.Process.GromacsHREX(
        perturbable_system, protocol, work_dir=str(tmp_path)
    )

    work_dir = str(process._work_dir)

    # Shared topology.
    assert os.path.isfile(os.path.join(work_dir, "gromacs.top"))

    # Per-lambda files.
    for lam in lam_vals:
        lam_dir = os.path.join(work_dir, f"lambda_{lam:5.4f}")
        assert os.path.isdir(lam_dir), f"Missing lambda dir: {lam_dir}"
        assert os.path.isfile(os.path.join(lam_dir, "gromacs.gro"))
        assert os.path.isfile(os.path.join(lam_dir, "gromacs.mdp"))
        assert os.path.isfile(os.path.join(lam_dir, "gromacs.tpr"))


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed.")
def test_topology_written_once(perturbable_system, tmp_path):
    """Only one .top file should exist at the top level, not inside lambda dirs."""
    protocol = BSS.Protocol.FreeEnergy()

    process = BSS.Process.GromacsHREX(
        perturbable_system, protocol, work_dir=str(tmp_path)
    )

    work_dir = str(process._work_dir)

    top_files = []
    for root, dirs, files in os.walk(work_dir):
        for fname in files:
            if fname.endswith(".top"):
                top_files.append(os.path.join(root, fname))

    assert len(top_files) == 1
    assert os.path.dirname(top_files[0]) == work_dir


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed.")
def test_lambda_dirs_are_absolute(perturbable_system, tmp_path):
    """_lambda_dirs should all be absolute paths."""
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.GromacsHREX(
        perturbable_system, protocol, work_dir=str(tmp_path)
    )

    for d in process._lambda_dirs:
        assert os.path.isabs(d), f"Lambda dir is not absolute: {d}"


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed.")
def test_thread_mpi_command(perturbable_system, tmp_path):
    """With use_mpi=False, the launch command should use -ntmpi and not mpirun."""
    protocol = BSS.Protocol.FreeEnergy()
    n_lam = len(protocol.getLambdaValues())

    process = BSS.Process.GromacsHREX(
        perturbable_system, protocol, use_mpi=False, work_dir=str(tmp_path)
    )

    # Trigger command construction without actually running.
    process.start()
    process.kill()

    assert f"-ntmpi {n_lam}" in process._command
    assert "mpirun" not in process._command
    assert f"-replex {process._repex_frequency}" in process._command


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed.")
def test_replex_flag_in_command(perturbable_system, tmp_path):
    """The -replex flag with the correct frequency must appear in the command."""
    freq = 500
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.GromacsHREX(
        perturbable_system, protocol, repex_frequency=freq, work_dir=str(tmp_path)
    )

    process.start()
    process.kill()

    assert f"-replex {freq}" in process._command


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed.")
def test_get_system_no_output(perturbable_system, tmp_path):
    """getSystem() returns None when no output coordinate files exist yet."""
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.GromacsHREX(
        perturbable_system, protocol, work_dir=str(tmp_path)
    )

    result = process.getSystem(block=False)
    assert result is None


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed.")
def test_get_exchange_statistics_no_log(perturbable_system, tmp_path):
    """getExchangeStatistics() returns None when no log file exists."""
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.GromacsHREX(
        perturbable_system, protocol, work_dir=str(tmp_path)
    )

    assert process.getExchangeStatistics() is None


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed.")
def test_get_exchange_statistics_parsed(perturbable_system, tmp_path):
    """getExchangeStatistics() correctly parses a synthetic GROMACS log."""
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.GromacsHREX(
        perturbable_system, protocol, work_dir=str(tmp_path)
    )

    # Write a synthetic log with exchange statistics to the first lambda dir.
    log_content = textwrap.dedent("""\
        Replica exchange statistics
        Repl  0 <-> 1  : accepted  45 out of 100
        Repl  1 <-> 2  : accepted  30 out of 100
        """)
    log_file = os.path.join(process._lambda_dirs[0], "gromacs.log")
    with open(log_file, "w") as f:
        f.write(log_content)

    df = process.getExchangeStatistics()

    assert df is not None
    assert list(df.columns) == [
        "replica_i",
        "replica_j",
        "lambda_i",
        "lambda_j",
        "n_attempts",
        "n_accepted",
        "acceptance_rate",
    ]
    assert len(df) == 2
    assert df.iloc[0]["replica_i"] == 0
    assert df.iloc[0]["replica_j"] == 1
    assert df.iloc[0]["n_accepted"] == 45
    assert df.iloc[0]["n_attempts"] == 100
    assert df.iloc[0]["acceptance_rate"] == pytest.approx(0.45)
    assert df.iloc[1]["n_accepted"] == 30
