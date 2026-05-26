"""Tests for FreeEnergy.Relative with replica_exchange=True."""

import socket

import pytest

import BioSimSpace as BSS
from tests.conftest import has_gromacs

_PMEMD_MPI = "/home/lester/.conda/envs/pmemd/bin/pmemd.MPI"


# ---------------------------------------------------------------------------
# Input validation tests: no engine required.
# ---------------------------------------------------------------------------


def test_repex_invalid_type(perturbable_system):
    """Non-bool repex raises TypeError."""
    with pytest.raises(TypeError, match="'repex'"):
        BSS.FreeEnergy.Relative(perturbable_system, engine="GROMACS", repex="yes")


def test_repex_frequency_invalid_zero(perturbable_system):
    """repex_frequency=0 raises ValueError."""
    with pytest.raises(ValueError, match="'repex_frequency'"):
        BSS.FreeEnergy.Relative(
            perturbable_system, engine="GROMACS", repex=True, repex_frequency=0
        )


def test_repex_frequency_invalid_float(perturbable_system):
    """Float repex_frequency raises ValueError."""
    with pytest.raises(ValueError, match="'repex_frequency'"):
        BSS.FreeEnergy.Relative(
            perturbable_system, engine="GROMACS", repex=True, repex_frequency=1.5
        )


def test_repex_unsupported_engine(perturbable_system):
    """repex=True with SOMD raises NotImplementedError."""
    with pytest.raises(NotImplementedError, match="SOMD"):
        BSS.FreeEnergy.Relative(perturbable_system, engine="SOMD", repex=True)


def test_get_exchange_statistics_no_repex(perturbable_system, tmp_path):
    """getExchangeStatistics() returns None when repex=False."""
    relative = BSS.FreeEnergy.Relative(
        perturbable_system, engine="GROMACS", setup_only=True, work_dir=str(tmp_path)
    )
    assert relative.getExchangeStatistics() is None


# ---------------------------------------------------------------------------
# GROMACS repex tests: require GROMACS.
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS.")
def test_repex_creates_gromacs_hrex(perturbable_system, tmp_path):
    """repex=True with GROMACS creates a GromacsHREX process."""
    relative = BSS.FreeEnergy.Relative(
        perturbable_system,
        engine="GROMACS",
        repex=True,
        work_dir=str(tmp_path),
    )
    assert isinstance(relative._process, BSS.Process.GromacsHREX)


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS.")
def test_repex_frequency_propagated_gromacs(perturbable_system, tmp_path):
    """repex_frequency is propagated to the GromacsHREX process."""
    freq = 500
    relative = BSS.FreeEnergy.Relative(
        perturbable_system,
        engine="GROMACS",
        repex=True,
        repex_frequency=freq,
        work_dir=str(tmp_path),
    )
    assert relative._process._repex_frequency == freq


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS.")
def test_use_mpi_propagated_gromacs(perturbable_system, tmp_path):
    """use_mpi kwarg is forwarded to GromacsHREX."""
    relative = BSS.FreeEnergy.Relative(
        perturbable_system,
        engine="GROMACS",
        repex=True,
        work_dir=str(tmp_path),
        use_mpi=False,
    )
    assert relative._process._use_mpi is False


# ---------------------------------------------------------------------------
# AMBER repex tests: require pmemd.MPI on hulk.
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    socket.gethostname() != "hulk",
    reason="Local test requiring pmemd.MPI installation.",
)
def test_repex_creates_amber_hrex(perturbable_system, tmp_path):
    """repex=True with AMBER creates an AmberHREX process."""
    relative = BSS.FreeEnergy.Relative(
        perturbable_system,
        engine="AMBER",
        repex=True,
        work_dir=str(tmp_path),
        exe=_PMEMD_MPI,
    )
    assert isinstance(relative._process, BSS.Process.AmberHREX)


@pytest.mark.skipif(
    socket.gethostname() != "hulk",
    reason="Local test requiring pmemd.MPI installation.",
)
def test_repex_frequency_propagated_amber(perturbable_system, tmp_path):
    """repex_frequency is propagated to AmberHREX and appears in amber.cfg."""
    import os

    freq = 250
    protocol = BSS.Protocol.FreeEnergy()
    lam_vals = protocol.getLambdaValues()

    relative = BSS.FreeEnergy.Relative(
        perturbable_system,
        protocol=protocol,
        engine="AMBER",
        repex=True,
        repex_frequency=freq,
        work_dir=str(tmp_path),
        exe=_PMEMD_MPI,
    )

    assert relative._process._repex_frequency == freq

    for lam in lam_vals:
        cfg = os.path.join(str(tmp_path), f"lambda_{lam:5.4f}", "amber.cfg")
        content = open(cfg).read()
        assert f"nstlim={freq}" in content or f"nstlim = {freq}" in content
