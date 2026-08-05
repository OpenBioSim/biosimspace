"""Tests for BioSimSpace.Sandpit.Exscientia.Process.GromacsHREX."""

import os
import textwrap

import pytest

from BioSimSpace.Sandpit import Exscientia as BSS
from tests.Sandpit.Exscientia.conftest import has_gromacs

url = BSS.tutorialUrl()


# ---------------------------------------------------------------------------
# Shared fixtures.
# ---------------------------------------------------------------------------


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
def abfe_system_and_restraint():
    """A decoupled ligand + protein system with a multiple-distance restraint."""
    from BioSimSpace.Sandpit.Exscientia.FreeEnergy import Restraint
    from BioSimSpace.Sandpit.Exscientia.Align._decouple import decouple
    from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol
    from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
    from BioSimSpace.Sandpit.Exscientia.Units.Temperature import kelvin

    ligand = BSS.IO.readMolecules(
        [f"{url}/ligand01.prm7.bz2", f"{url}/ligand01.rst7.bz2"]
    ).getMolecule(0)
    decoupled_ligand = decouple(ligand)

    protein = BSS.IO.readMolecules(
        [f"{url}/1jr5.crd.bz2", f"{url}/1jr5.top.bz2"]
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


@pytest.fixture(scope="module")
def boresch_restraint(abfe_system_and_restraint):
    """A Boresch restraint (wrong type for release_restraint)."""
    from BioSimSpace.Sandpit.Exscientia.FreeEnergy import Restraint
    from BioSimSpace.Sandpit.Exscientia.Units.Angle import degree, radian
    from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol
    from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
    from BioSimSpace.Sandpit.Exscientia.Units.Temperature import kelvin

    system, _ = abfe_system_and_restraint

    ligand = system.getMolecule(-1)
    protein = system.getMolecule(0)

    restraint_dict = {
        "anchor_points": {
            "r1": protein.getAtoms()[0],
            "r2": protein.getAtoms()[1],
            "r3": protein.getAtoms()[2],
            "l1": ligand.getAtoms()[0],
            "l2": ligand.getAtoms()[1],
            "l3": ligand.getAtoms()[2],
        },
        "equilibrium_values": {
            "r0": 5.08 * angstrom,
            "thetaA0": 64.051 * degree,
            "thetaB0": 39.618 * degree,
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

    return Restraint(system, restraint_dict, 298 * kelvin, restraint_type="Boresch")


# ---------------------------------------------------------------------------
# Input validation tests: no GROMACS required.
# ---------------------------------------------------------------------------


def test_non_fep_protocol(perturbable_system):
    """Non-FEP protocol raises TypeError."""
    with pytest.raises(TypeError, match="_FreeEnergyMixin"):
        BSS.Process.GromacsHREX(
            perturbable_system, BSS.Protocol.Minimisation(steps=100)
        )


def test_invalid_oversubscribe(perturbable_system):
    """Non-bool oversubscribe raises TypeError."""
    protocol = BSS.Protocol.FreeEnergy()
    with pytest.raises(TypeError, match="oversubscribe"):
        BSS.Process.GromacsHREX(perturbable_system, protocol, oversubscribe="yes")


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
    from BioSimSpace.Sandpit.Exscientia._SireWrappers import ReplicaSystem

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


def test_unsupported_perturbation_type(perturbable_system):
    """Unsupported perturbation type raises NotImplementedError."""
    protocol = BSS.Protocol.FreeEnergy(perturbation_type="discharge_soft")
    with pytest.raises(NotImplementedError):
        BSS.Process.GromacsHREX(perturbable_system, protocol)


# ---------------------------------------------------------------------------
# Restraint validation tests: no GROMACS required.
# ---------------------------------------------------------------------------


def test_release_restraint_requires_restraint(abfe_system_and_restraint):
    """release_restraint without a restraint object raises ValueError."""
    system, _ = abfe_system_and_restraint
    protocol = BSS.Protocol.FreeEnergy(perturbation_type="release_restraint")
    with pytest.raises(ValueError, match="restraint"):
        BSS.Process.GromacsHREX(system, protocol)


def test_release_restraint_wrong_object_type(abfe_system_and_restraint):
    """Passing a non-Restraint as restraint raises TypeError."""
    system, _ = abfe_system_and_restraint
    protocol = BSS.Protocol.FreeEnergy(perturbation_type="release_restraint")
    with pytest.raises(TypeError, match="'restraint'"):
        BSS.Process.GromacsHREX(system, protocol, restraint="not-a-restraint")


def test_release_restraint_non_multiple_distance(
    abfe_system_and_restraint, boresch_restraint
):
    """A Boresch restraint (not multiple_distance) raises ValueError."""
    system, _ = abfe_system_and_restraint
    protocol = BSS.Protocol.FreeEnergy(perturbation_type="release_restraint")
    with pytest.raises(ValueError, match="multiple distance"):
        BSS.Process.GromacsHREX(system, protocol, restraint=boresch_restraint)


# ---------------------------------------------------------------------------
# Tests that require GROMACS to be installed.
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed.")
def test_system_promotion(perturbable_system):
    """Plain System input is promoted to a ReplicaSystem with the right replica count."""
    from BioSimSpace.Sandpit.Exscientia._SireWrappers import ReplicaSystem

    protocol = BSS.Protocol.FreeEnergy()
    n_lam = len(protocol.getLambdaValues())

    process = BSS.Process.GromacsHREX(perturbable_system, protocol)

    assert isinstance(process._replica_system, ReplicaSystem)
    assert process._replica_system.nReplicas() == n_lam


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed.")
def test_replica_system_passthrough(perturbable_system):
    """A ReplicaSystem with the correct replica count is accepted without re-wrapping."""
    from BioSimSpace.Sandpit.Exscientia._SireWrappers import ReplicaSystem

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

    assert os.path.isfile(os.path.join(work_dir, "gromacs.top"))

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
def test_oversubscribe_in_command(perturbable_system, tmp_path):
    """oversubscribe=True adds --oversubscribe to the mpirun command."""
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.GromacsHREX(
        perturbable_system, protocol, oversubscribe=True, work_dir=str(tmp_path)
    )
    process.start()
    process.kill()
    assert "--oversubscribe" in process._command


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed.")
def test_no_oversubscribe_by_default(perturbable_system, tmp_path):
    """oversubscribe defaults to False — --oversubscribe must not appear."""
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.GromacsHREX(
        perturbable_system, protocol, work_dir=str(tmp_path)
    )
    process.start()
    process.kill()
    assert "--oversubscribe" not in process._command


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed.")
def test_mpi_exe_in_command(perturbable_system, tmp_path):
    """mpi_exe is used as the MPI launcher instead of mpirun."""
    protocol = BSS.Protocol.FreeEnergy()
    process = BSS.Process.GromacsHREX(
        perturbable_system, protocol, mpi_exe="srun", work_dir=str(tmp_path)
    )
    process.start()
    process.kill()
    assert process._command.startswith("srun")


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


@pytest.mark.skipif(not has_gromacs, reason="Requires GROMACS to be installed.")
def test_abfe_restraint_stored(abfe_system_and_restraint, tmp_path):
    """GromacsHREX stores the restraint and writes restraint topology."""
    from BioSimSpace.Sandpit.Exscientia.FreeEnergy import Restraint

    system, restraint = abfe_system_and_restraint
    protocol = BSS.Protocol.FreeEnergy(
        perturbation_type="release_restraint",
        lam_vals=[0.0, 0.5, 1.0],
    )

    process = BSS.Process.GromacsHREX(
        system, protocol, restraint=restraint, work_dir=str(tmp_path)
    )

    assert isinstance(process._restraint, Restraint)

    top_file = os.path.join(str(process._work_dir), "gromacs.top")
    with open(top_file) as f:
        top_content = f.read()

    assert "intermolecular_interactions" in top_content
