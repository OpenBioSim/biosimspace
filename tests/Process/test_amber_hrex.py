"""Tests for BioSimSpace.Process.AmberHREX."""

import os
import socket
import textwrap

import pytest

import BioSimSpace as BSS

# Full path to the pmemd.MPI binary on hosts where these tests run.
_PMEMD_MPI = "/home/lester/.conda/envs/pmemd/bin/pmemd.MPI"


# ---------------------------------------------------------------------------
# Input validation tests: errors raised before the exe lookup in
# Amber.__init__(), so no AMBER installation is required.
# ---------------------------------------------------------------------------


def test_non_fep_protocol(perturbable_system):
    """Non-FEP protocol raises TypeError."""
    with pytest.raises(TypeError, match="_FreeEnergyMixin"):
        BSS.Process.AmberHREX(perturbable_system, BSS.Protocol.Minimisation(steps=100))


def test_invalid_repex_frequency_zero(perturbable_system):
    """repex_frequency=0 raises ValueError."""
    protocol = BSS.Protocol.FreeEnergy()
    with pytest.raises(ValueError, match="repex_frequency"):
        BSS.Process.AmberHREX(perturbable_system, protocol, repex_frequency=0)


def test_invalid_repex_frequency_float(perturbable_system):
    """Float repex_frequency raises ValueError."""
    protocol = BSS.Protocol.FreeEnergy()
    with pytest.raises(ValueError, match="repex_frequency"):
        BSS.Process.AmberHREX(perturbable_system, protocol, repex_frequency=1.5)


def test_invalid_is_gpu(perturbable_system):
    """Non-bool is_gpu raises TypeError."""
    protocol = BSS.Protocol.FreeEnergy()
    with pytest.raises(TypeError, match="is_gpu"):
        BSS.Process.AmberHREX(perturbable_system, protocol, is_gpu="yes")


def test_replica_count_mismatch(perturbable_system):
    """ReplicaSystem with wrong replica count raises ValueError."""
    from BioSimSpace._SireWrappers import ReplicaSystem

    protocol = BSS.Protocol.FreeEnergy()
    n_lam = len(protocol.getLambdaValues())
    wrong = ReplicaSystem(perturbable_system, num_replicas=n_lam + 2)
    with pytest.raises(ValueError, match="number of replicas"):
        BSS.Process.AmberHREX(wrong, protocol)


def test_system_type_error():
    """Passing a non-System, non-ReplicaSystem raises TypeError."""
    protocol = BSS.Protocol.FreeEnergy()
    with pytest.raises(TypeError, match="'system' must be"):
        BSS.Process.AmberHREX("not-a-system", protocol)


# ---------------------------------------------------------------------------
# Tests that require AMBER to be installed.
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    socket.gethostname() != "hulk",
    reason="Local test requiring pmemd.MPI installation.",
)
def test_system_promotion(perturbable_system):
    """Plain System input is promoted to a ReplicaSystem with the right replica count."""
    from BioSimSpace._SireWrappers import ReplicaSystem

    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.AmberHREX(perturbable_system, protocol, exe=_PMEMD_MPI)

    n_lam = len(protocol.getLambdaValues())

    assert isinstance(process._replica_system, ReplicaSystem)
    assert process._replica_system.nReplicas() == n_lam


@pytest.mark.skipif(
    socket.gethostname() != "hulk",
    reason="Local test requiring pmemd.MPI installation.",
)
def test_replica_system_passthrough(perturbable_system):
    """A ReplicaSystem with the correct replica count is accepted."""
    from BioSimSpace._SireWrappers import ReplicaSystem

    protocol = BSS.Protocol.FreeEnergy()
    n_lam = len(protocol.getLambdaValues())
    rs = ReplicaSystem(perturbable_system, num_replicas=n_lam)
    process = BSS.Process.AmberHREX(rs, protocol, exe=_PMEMD_MPI)

    assert isinstance(process._replica_system, ReplicaSystem)
    assert process._replica_system.nReplicas() == n_lam


@pytest.mark.skipif(
    socket.gethostname() != "hulk",
    reason="Local test requiring pmemd.MPI installation.",
)
def test_file_layout(perturbable_system, tmp_path):
    """After setup, the shared topology, groupfile, and per-lambda inputs must exist."""
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.AmberHREX(
        perturbable_system, protocol, exe=_PMEMD_MPI, work_dir=str(tmp_path)
    )

    lam_vals = protocol.getLambdaValues()

    work_dir = str(process._work_dir)

    # Shared topology.
    assert os.path.isfile(os.path.join(work_dir, "amber.prm7"))

    # Groupfile.
    assert os.path.isfile(os.path.join(work_dir, "groupfile"))

    # Per-lambda files.
    for lam in lam_vals:
        lam_dir = os.path.join(work_dir, f"lambda_{lam:5.4f}")
        assert os.path.isdir(lam_dir), f"Missing lambda dir: {lam_dir}"
        assert os.path.isfile(os.path.join(lam_dir, "amber.rst7"))
        assert os.path.isfile(os.path.join(lam_dir, "amber.cfg"))


@pytest.mark.skipif(
    socket.gethostname() != "hulk",
    reason="Local test requiring pmemd.MPI installation.",
)
def test_topology_written_once(perturbable_system, tmp_path):
    """Only one .prm7 file should exist at the top level."""
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.AmberHREX(
        perturbable_system, protocol, exe=_PMEMD_MPI, work_dir=str(tmp_path)
    )

    work_dir = str(process._work_dir)

    prm7_files = []
    for root, dirs, files in os.walk(work_dir):
        for fname in files:
            if fname.endswith(".prm7"):
                prm7_files.append(os.path.join(root, fname))

    assert len(prm7_files) == 1
    assert os.path.dirname(prm7_files[0]) == work_dir


@pytest.mark.skipif(
    socket.gethostname() != "hulk",
    reason="Local test requiring pmemd.MPI installation.",
)
def test_groupfile_content(perturbable_system, tmp_path):
    """Each line in the groupfile references a per-lambda cfg, prm7, rst7, and out file."""
    protocol = BSS.Protocol.FreeEnergy()
    lam_vals = protocol.getLambdaValues()
    process = BSS.Process.AmberHREX(
        perturbable_system, protocol, exe=_PMEMD_MPI, work_dir=str(tmp_path)
    )

    groupfile = os.path.join(str(process._work_dir), "groupfile")
    with open(groupfile) as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    assert len(lines) == len(lam_vals)

    for line, lam in zip(lines, lam_vals):
        lam_name = f"lambda_{lam:5.4f}"
        assert f"-i {lam_name}/amber.cfg" in line
        assert "-p amber.prm7" in line
        assert f"-c {lam_name}/amber.rst7" in line
        assert f"-o {lam_name}/amber.out" in line


@pytest.mark.skipif(
    socket.gethostname() != "hulk",
    reason="Local test requiring pmemd.MPI installation.",
)
def test_groupfile_paths_relative(perturbable_system, tmp_path):
    """Groupfile paths must be relative (not absolute) so pmemd.MPI can resolve them."""
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.AmberHREX(
        perturbable_system, protocol, exe=_PMEMD_MPI, work_dir=str(tmp_path)
    )

    groupfile = os.path.join(str(process._work_dir), "groupfile")
    with open(groupfile) as f:
        content = f.read()

    # No line should start with '/' after the flag prefix.
    for line in content.splitlines():
        for token in line.split():
            if token not in ("-O", "-i", "-p", "-c", "-o", "-r", "-x", "-inf"):
                assert not os.path.isabs(token), (
                    f"Groupfile contains absolute path: {token}"
                )


@pytest.mark.skipif(
    socket.gethostname() != "hulk",
    reason="Local test requiring pmemd.MPI installation.",
)
def test_nstlim_in_cfg(perturbable_system, tmp_path):
    """Each amber.cfg must contain nstlim equal to repex_frequency."""
    freq = 250
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.AmberHREX(
        perturbable_system,
        protocol,
        exe=_PMEMD_MPI,
        repex_frequency=freq,
        work_dir=str(tmp_path),
    )

    lam_vals = protocol.getLambdaValues()

    work_dir = str(process._work_dir)
    for lam in lam_vals:
        cfg = os.path.join(work_dir, f"lambda_{lam:5.4f}", "amber.cfg")
        content = open(cfg).read()
        assert f"nstlim={freq}" in content or f"nstlim = {freq}" in content


@pytest.mark.skipif(
    socket.gethostname() != "hulk",
    reason="Local test requiring pmemd.MPI installation.",
)
def test_get_system_no_output(perturbable_system, tmp_path):
    """getSystem() returns None when no output coordinate files exist yet."""
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.AmberHREX(
        perturbable_system, protocol, exe=_PMEMD_MPI, work_dir=str(tmp_path)
    )

    result = process.getSystem(block=False)
    assert result is None


@pytest.mark.skipif(
    socket.gethostname() != "hulk",
    reason="Local test requiring pmemd.MPI installation.",
)
def test_get_exchange_statistics_no_remlog(perturbable_system, tmp_path):
    """getExchangeStatistics() returns None when no remlog file exists."""
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.AmberHREX(
        perturbable_system, protocol, exe=_PMEMD_MPI, work_dir=str(tmp_path)
    )

    assert process.getExchangeStatistics() is None


@pytest.mark.skipif(
    socket.gethostname() != "hulk",
    reason="Local test requiring pmemd.MPI installation.",
)
def test_get_exchange_statistics_parsed(perturbable_system, tmp_path):
    """getExchangeStatistics() correctly parses a synthetic AMBER remlog."""
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.AmberHREX(
        perturbable_system, protocol, exe=_PMEMD_MPI, work_dir=str(tmp_path)
    )

    # Write a synthetic remlog with exchange attempt records.
    log_content = textwrap.dedent("""\
        # Replica exchange log
        exchange:  0 (0.00000) <->  1 (0.10000) : success
        exchange:  1 (0.10000) <->  2 (0.20000) : failed
        exchange:  0 (0.00000) <->  1 (0.10000) : success
        exchange:  1 (0.10000) <->  2 (0.20000) : success
        """)
    remlog = os.path.join(str(process._work_dir), "remlog")
    with open(remlog, "w") as f:
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

    # Pair (0,1): 2 attempts, 2 successes.
    row01 = df[df["replica_i"] == 0].iloc[0]
    assert row01["n_attempts"] == 2
    assert row01["n_accepted"] == 2
    assert row01["acceptance_rate"] == pytest.approx(1.0)

    # Pair (1,2): 2 attempts, 1 success.
    row12 = df[df["replica_i"] == 1].iloc[0]
    assert row12["n_attempts"] == 2
    assert row12["n_accepted"] == 1
    assert row12["acceptance_rate"] == pytest.approx(0.5)
