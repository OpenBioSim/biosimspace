"""Tests for FreeEnergy.AlchemicalFreeEnergy with repex=True."""

import socket

import pytest

from BioSimSpace.Sandpit import Exscientia as BSS
from tests.Sandpit.Exscientia.conftest import has_gromacs

_PMEMD_MPI = "/home/lester/.conda/envs/pmemd/bin/pmemd.MPI"
_url = BSS.tutorialUrl()


@pytest.fixture()
def perturbable_system():
    """A vacuum perturbable system (same files as main conftest)."""
    return BSS.IO.readPerturbableSystem(
        f"{_url}/perturbable_system0.prm7",
        f"{_url}/perturbable_system0.rst7",
        f"{_url}/perturbable_system1.prm7",
        f"{_url}/perturbable_system1.rst7",
    )


@pytest.fixture(scope="module")
def abfe_system_and_restraint():
    """A decoupled ligand + protein system with a multiple-distance restraint."""
    from BioSimSpace.Sandpit.Exscientia.FreeEnergy import Restraint
    from BioSimSpace.Sandpit.Exscientia.Align._decouple import decouple
    from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol
    from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
    from BioSimSpace.Sandpit.Exscientia.Units.Temperature import kelvin

    ligand = BSS.IO.readMolecules(
        [f"{_url}/ligand01.prm7.bz2", f"{_url}/ligand01.rst7.bz2"]
    ).getMolecule(0)
    decoupled_ligand = decouple(ligand)

    protein = BSS.IO.readMolecules(
        [f"{_url}/1jr5.crd.bz2", f"{_url}/1jr5.top.bz2"]
    ).getMolecule(0)

    system = (protein + decoupled_ligand).toSystem()

    restraint_dict = {
        "distance_restraints": [
            {
                "l1": decoupled_ligand.getAtoms()[0],
                "r1": protein.getAtoms()[0],
                "r0": 3 * angstrom,
                "kr": 10 * kcal_per_mol / angstrom**2,
                "r_fb": 1 * angstrom,
            },
            {
                "l1": decoupled_ligand.getAtoms()[1],
                "r1": protein.getAtoms()[1],
                "r0": 3 * angstrom,
                "kr": 10 * kcal_per_mol / angstrom**2,
                "r_fb": 1 * angstrom,
            },
        ],
        "permanent_distance_restraint": {
            "l1": decoupled_ligand.getAtoms()[2],
            "r1": protein.getAtoms()[2],
            "r0": 3 * angstrom,
            "kr": 10 * kcal_per_mol / angstrom**2,
            "r_fb": 1 * angstrom,
        },
    }

    restraint = Restraint(
        system, restraint_dict, 300 * kelvin, restraint_type="multiple_distance"
    )
    return system, restraint


# ---------------------------------------------------------------------------
# Input validation tests: no engine required.
# ---------------------------------------------------------------------------


def test_repex_invalid_type(perturbable_system):
    """Non-bool repex raises TypeError."""
    with pytest.raises(TypeError, match="'repex'"):
        BSS.FreeEnergy.AlchemicalFreeEnergy(
            perturbable_system, engine="GROMACS", repex="yes"
        )


def test_repex_frequency_invalid_zero(perturbable_system):
    """repex_frequency=0 raises ValueError."""
    with pytest.raises(ValueError, match="'repex_frequency'"):
        BSS.FreeEnergy.AlchemicalFreeEnergy(
            perturbable_system, engine="GROMACS", repex=True, repex_frequency=0
        )


def test_repex_frequency_invalid_float(perturbable_system):
    """Float repex_frequency raises ValueError."""
    with pytest.raises(ValueError, match="'repex_frequency'"):
        BSS.FreeEnergy.AlchemicalFreeEnergy(
            perturbable_system, engine="GROMACS", repex=True, repex_frequency=1.5
        )


def test_repex_unsupported_engine(perturbable_system):
    """repex=True with SOMD raises NotImplementedError."""
    with pytest.raises(NotImplementedError, match="SOMD"):
        BSS.FreeEnergy.AlchemicalFreeEnergy(
            perturbable_system, engine="SOMD", repex=True
        )


def test_get_exchange_statistics_no_repex(perturbable_system, tmp_path):
    """getExchangeStatistics() returns None when repex=False."""
    relative = BSS.FreeEnergy.AlchemicalFreeEnergy(
        perturbable_system, engine="GROMACS", setup_only=True, work_dir=str(tmp_path)
    )
    assert relative.getExchangeStatistics() is None


# ---------------------------------------------------------------------------
# GROMACS repex tests: require GROMACS.
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS.")
def test_repex_creates_gromacs_hrex(perturbable_system, tmp_path):
    """repex=True with GROMACS creates a GromacsHREX process."""
    relative = BSS.FreeEnergy.AlchemicalFreeEnergy(
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
    relative = BSS.FreeEnergy.AlchemicalFreeEnergy(
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
    relative = BSS.FreeEnergy.AlchemicalFreeEnergy(
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
    relative = BSS.FreeEnergy.AlchemicalFreeEnergy(
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

    relative = BSS.FreeEnergy.AlchemicalFreeEnergy(
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


# ---------------------------------------------------------------------------
# ABFE + repex restraint validation tests: no engine required.
# ---------------------------------------------------------------------------


def test_repex_release_restraint_no_restraint(abfe_system_and_restraint):
    """repex=True with release_restraint but no restraint raises ValueError."""
    system, _ = abfe_system_and_restraint
    protocol = BSS.Protocol.FreeEnergy(perturbation_type="release_restraint")
    with pytest.raises(ValueError, match="restraint"):
        BSS.FreeEnergy.AlchemicalFreeEnergy(
            system, protocol=protocol, engine="GROMACS", repex=True
        )


def test_repex_release_restraint_wrong_object_type(abfe_system_and_restraint):
    """repex=True with a non-Restraint object raises TypeError."""
    system, _ = abfe_system_and_restraint
    protocol = BSS.Protocol.FreeEnergy(perturbation_type="release_restraint")
    with pytest.raises(TypeError, match="'restraint'"):
        BSS.FreeEnergy.AlchemicalFreeEnergy(
            system,
            protocol=protocol,
            engine="GROMACS",
            repex=True,
            restraint="not-a-restraint",
        )


# ---------------------------------------------------------------------------
# ABFE + repex GROMACS integration tests: require GROMACS.
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS.")
def test_repex_restraint_forwarded_to_process(abfe_system_and_restraint, tmp_path):
    """The restraint passed to AlchemicalFreeEnergy is forwarded to GromacsHREX."""
    from BioSimSpace.Sandpit.Exscientia.FreeEnergy import Restraint

    system, restraint = abfe_system_and_restraint
    protocol = BSS.Protocol.FreeEnergy(
        perturbation_type="release_restraint",
        lam_vals=[0.0, 0.5, 1.0],
    )
    relative = BSS.FreeEnergy.AlchemicalFreeEnergy(
        system,
        protocol=protocol,
        engine="GROMACS",
        repex=True,
        restraint=restraint,
        work_dir=str(tmp_path),
    )

    assert isinstance(relative._process, BSS.Process.GromacsHREX)
    assert isinstance(relative._process._restraint, Restraint)


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS.")
def test_repex_restraint_topology_written(abfe_system_and_restraint, tmp_path):
    """The ABFE restraint is written into the shared GROMACS topology."""
    import os

    system, restraint = abfe_system_and_restraint
    protocol = BSS.Protocol.FreeEnergy(
        perturbation_type="release_restraint",
        lam_vals=[0.0, 0.5, 1.0],
    )
    relative = BSS.FreeEnergy.AlchemicalFreeEnergy(
        system,
        protocol=protocol,
        engine="GROMACS",
        repex=True,
        restraint=restraint,
        work_dir=str(tmp_path),
    )

    top_file = os.path.join(str(tmp_path), "gromacs.top")
    with open(top_file) as f:
        top_content = f.read()

    assert "intermolecular_interactions" in top_content
