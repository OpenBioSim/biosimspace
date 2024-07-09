import itertools
import pandas as pd
import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Process._process import Process

from tests.Sandpit.Exscientia.conftest import (
    has_amber,
    has_antechamber,
    has_gromacs,
    has_openff,
)


@pytest.fixture
def system():
    ff = "openff_unconstrained-2.0.0"
    mol0 = BSS.Parameters.parameterise("c1ccccc1C", ff).getMolecule()
    mol1 = BSS.Parameters.parameterise("c1ccccc1", ff).getMolecule()
    return BSS.Align.merge(mol0, mol1).toSystem()


@pytest.fixture
def system_vel(system):
    mol = system.getMolecule(0)
    mol_sire = mol._sire_object

    # Edit the molecule
    mol_edit = mol_sire.edit()

    mol_edit.setProperty("velocity", 1)

    # Update the Sire molecule object of the new molecule.
    mol._sire_object = mol_edit.commit()
    return mol.toSystem()


@pytest.fixture
def system_vel0(system):
    mol = system.getMolecule(0)
    mol_sire = mol._sire_object

    # Edit the molecule
    mol_edit = mol_sire.edit()

    mol_edit.setProperty("velocity0", 1)

    # Update the Sire molecule object of the new molecule.
    mol._sire_object = mol_edit.commit()
    return mol.toSystem()


@pytest.fixture
def free_energy():
    return {
        "lam": pd.Series(data={"fep": 0.5}),
        "lam_vals": pd.DataFrame(
            data={"fep": [0.5, 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0]}
        ),
    }


@pytest.fixture()
def equilibration(production):
    return {
        **production,
        "temperature_start": BSS.Types.Temperature(300, "kelvin"),
        "temperature_end": BSS.Types.Temperature(310, "kelvin"),
    }


@pytest.fixture()
def production():
    return {"runtime": BSS.Types.Time(0.001, "nanoseconds")}


@pytest.fixture(
    params=itertools.product(
        [True, False], [True, False], ["Equilibration", "Production"]
    )
)
def protocol(request, free_energy, equilibration, production):
    FE, restart, name = request.param
    if FE:
        if name == "Equilibration":
            return BSS.Protocol.FreeEnergyEquilibration(
                **equilibration, **free_energy, restart=restart
            )
        else:
            return BSS.Protocol.FreeEnergy(**production, **free_energy, restart=restart)
    else:
        if name == "Equilibration":
            return BSS.Protocol.Equilibration(**equilibration, restart=restart)
        else:
            return BSS.Protocol.Production(**production, restart=restart)


@pytest.mark.skipif(
    has_gromacs is False or has_openff is False,
    reason="Requires GROMACS and OpenFF to be installed.",
)
def test_gromacs(protocol, system, tmp_path):
    BSS.Process.Gromacs(system, protocol, work_dir=str(tmp_path), ignore_warnings=True)
    with open(tmp_path / "gromacs.mdp", "r") as f:
        cfg = f.read()
    if protocol.isRestart():
        assert "gen-vel" not in cfg
    else:
        assert "gen-vel" in cfg
        assert "gen-temp" in cfg


@pytest.mark.skipif(
    has_amber is False or has_openff is False,
    reason="Requires AMBER and OpenFF to be installed.",
)
def test_amber(protocol, system, tmp_path):
    if not isinstance(protocol, BSS.Protocol._FreeEnergyMixin):
        pytest.skip("AMBER position restraint only works for free energy protocol")
    BSS.Process.Amber(system, protocol, work_dir=str(tmp_path))
    with open(tmp_path / "amber.cfg", "r") as f:
        cfg = f.read()
        if protocol.isRestart():
            assert "irest=1" in cfg
            assert "ntx=5" in cfg
        else:
            assert "irest=0" in cfg
            assert "ntx=1" in cfg


@pytest.mark.skipif(
    has_antechamber is False or has_openff is False,
    reason="Requires AmberTools/antechamber and OpenFF to be installed.",
)
@pytest.mark.parametrize("restart", [True, False])
@pytest.mark.parametrize("name", ["system", "system_vel", "system_vel0"])
@pytest.mark.parametrize(
    "protocol", [BSS.Protocol.Production, BSS.Protocol.Equilibration]
)
def test_process(protocol, name, request, restart, production):
    """Ensure that when there is no velocity in the system. The restart is False."""
    _protocol = protocol(**production, restart=restart)
    process = Process(request.getfixturevalue(name), _protocol)
    if name == "system":
        assert process._protocol.isRestart() is False
    else:
        assert process._protocol.isRestart() is restart
