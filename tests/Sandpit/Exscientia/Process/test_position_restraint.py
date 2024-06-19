import itertools
import os
from difflib import unified_diff

import pandas as pd
import pytest
import sire as sr
from sire.legacy import Units as SireUnits
from sire.legacy.IO import AmberRst

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align._alch_ion import _mark_alchemical_ion
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kj_per_mol
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
from BioSimSpace.Sandpit.Exscientia._SireWrappers import Molecule
from tests.conftest import root_fp


@pytest.fixture
def system():
    ff = "openff_unconstrained-2.0.0"
    mol0 = BSS.Parameters.parameterise("c1ccccc1C", ff).getMolecule()
    mol1 = BSS.Parameters.parameterise("c1ccccc1", ff).getMolecule()
    return BSS.Align.merge(mol0, mol1).toSystem()


@pytest.fixture(scope="session")
def alchemical_ion_system():
    """A large protein system for re-use."""
    system = BSS.IO.readMolecules(
        [f"{root_fp}/input/ala.top", f"{root_fp}/input/ala.crd"]
    )
    solvated = BSS.Solvent.tip3p(
        system, box=[4 * BSS.Units.Length.nanometer] * 3, ion_conc=0.15
    )
    ion = solvated.getMolecule(-1)
    pert_ion = BSS.Align.merge(ion, ion, mapping={0: 0})
    for lambda_ in [0, 1]:
        atomtype = pert_ion._sire_object.property(f"atomtype{lambda_}")
        pert_ion._sire_object = (
            pert_ion._sire_object.edit()
            .setProperty(f"ambertype{lambda_}", atomtype)
            .molecule()
        )
    pert_ion._sire_object = (
        pert_ion.getAtoms()[0]
        ._sire_object.edit()
        .setProperty("charge1", 0 * SireUnits.mod_electron)
        .molecule()
    )

    alchemcial_ion = _mark_alchemical_ion(pert_ion)
    solvated.updateMolecule(solvated.getIndex(ion), alchemcial_ion)
    return solvated


@pytest.fixture(scope="session")
def alchemical_ion_system_psores(alchemical_ion_system):
    # Create a reference system with different coordinate
    system = alchemical_ion_system.copy()
    mol = system.getMolecule(0)
    sire_mol = mol._sire_object
    atoms = sire_mol.cursor().atoms()
    atoms[0]["coordinates"] = sr.maths.Vector(0, 0, 0)
    new_mol = atoms.commit()
    system.updateMolecule(0, Molecule(new_mol))
    return system


@pytest.fixture
def ref_system(system):
    mol = system[0]
    mol.translate(3 * [BSS.Units.Length.nanometer * 1])
    return mol.toSystem()


@pytest.fixture
def restraint():
    return {
        "restraint": "heavy",
        "force_constant": 1000 * kj_per_mol / angstrom**2,
    }


@pytest.fixture
def free_energy():
    return {
        "lam": pd.Series(data={"fep": 0.5}),
        "lam_vals": pd.DataFrame(
            data={"fep": [0.5, 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0]}
        ),
    }


@pytest.fixture()
def minimisation():
    return {"steps": 100}


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
        [True, False], ["Minimisation", "Equilibration", "Production"]
    )
)
def protocol(request, restraint, free_energy, minimisation, equilibration, production):
    FE, name = request.param
    if FE:
        if name == "Minimisation":
            return BSS.Protocol.FreeEnergyMinimisation(
                **minimisation, **restraint, **free_energy
            )
        elif name == "Equilibration":
            return BSS.Protocol.FreeEnergyEquilibration(
                **equilibration, **restraint, **free_energy
            )
        else:
            return BSS.Protocol.FreeEnergy(**production, **restraint, **free_energy)
    else:
        if name == "Minimisation":
            return BSS.Protocol.Minimisation(**minimisation, **restraint)
        elif name == "Equilibration":
            return BSS.Protocol.Equilibration(**equilibration, **restraint)
        else:
            return BSS.Protocol.Production(**production, **restraint)


def test_gromacs(protocol, system, ref_system, tmp_path):
    proc = BSS.Process.Gromacs(
        system,
        protocol,
        reference_system=ref_system,
        work_dir=str(tmp_path),
        ignore_warnings=True,
    )

    # We have the restraints
    assert (tmp_path / "posre_0001.itp").is_file()
    with open(tmp_path / "gromacs.top", "r") as f:
        assert "posre_0001.itp" in f.read()

    # We have generated a separate restraint reference
    assert os.path.exists(proc._ref_file)
    contents_ref, contents_gro = (
        open(proc._ref_file).readlines(),
        open(proc._gro_file).readlines(),
    )
    diff = list(unified_diff(contents_ref, contents_gro))
    assert len(diff)


def test_amber(protocol, system, ref_system, tmp_path):
    if not isinstance(protocol, BSS.Protocol._FreeEnergyMixin):
        pytest.skip("AMBER position restraint only works for free energy protocol")
    proc = BSS.Process.Amber(
        system, protocol, reference_system=ref_system, work_dir=str(tmp_path)
    )

    # We have the restraints
    with open(tmp_path / "amber.cfg", "r") as f:
        cfg = f.read()
        assert "restraint_wt" in cfg
        assert "restraintmask" in cfg

    # We have generated a separate restraint reference
    assert os.path.exists(proc._ref_file)

    ref = AmberRst(proc._ref_file).getFrame(0)
    rst = AmberRst(proc._rst_file).getFrame(0)

    assert ref != rst

    # We are pointing the reference to the correct file
    assert f"{proc._work_dir}/{proc.getArgs()['-ref']}" == proc._ref_file


@pytest.mark.parametrize(
    "restraint",
    ["backbone", "heavy", "all", "none"],
)
def test_gromacs_alchemical_ion(
    alchemical_ion_system, restraint, alchemical_ion_system_psores
):
    protocol = BSS.Protocol.FreeEnergy(restraint=restraint)
    process = BSS.Process.Gromacs(
        alchemical_ion_system,
        protocol,
        name="test",
        reference_system=alchemical_ion_system_psores,
        ignore_warnings=True,
    )

    # Test the position restraint for protein center
    with open(f"{process.workDir()}/posre_0001.itp", "r") as f:
        posres = f.read().split("\n")
        posres = [tuple(line.split()) for line in posres]

    assert ("9", "1", "4184.0", "4184.0", "4184.0") in posres

    # Test the position restraint for alchemical ion
    with open(f"{process.workDir()}/test.top", "r") as f:
        top = f.read()
    lines = top[top.index("Merged_Molecule") :].split("\n")
    assert lines[6] == '#include "posre_0002.itp"'

    with open(f"{process.workDir()}/posre_0002.itp", "r") as f:
        posres = f.read().split("\n")

    assert posres[2].split() == ["1", "1", "4184.0", "4184.0", "4184.0"]

    # Test if the original coordinate is correct
    with open(f"{process.workDir()}/test.gro", "r") as f:
        gro = f.read().splitlines()
    assert gro[2].split() == ["1ACE", "HH31", "1", "1.791", "1.610", "2.058"]

    # Test if the reference coordinate is passed
    with open(f"{process.workDir()}/test_ref.gro", "r") as f:
        gro = f.read().splitlines()
    assert gro[2].split() == ["1ACE", "HH31", "1", "0.000", "0.000", "0.000"]


@pytest.mark.parametrize(
    ("restraint", "protocol", "target"),
    [
        (
            "backbone",
            BSS.Protocol.FreeEnergyEquilibration,
            "@5-7,9,15-17 | @9,6442,6443",
        ),
        (
            "heavy",
            BSS.Protocol.FreeEnergyEquilibration,
            "@2,5-7,9,11,15-17,19 | @9,6442,6443",
        ),
        ("all", BSS.Protocol.FreeEnergyEquilibration, "@1-22 | @9,6442,6443"),
        ("none", BSS.Protocol.FreeEnergyEquilibration, "@9,6442,6443"),
    ],
)
def test_amber_alchemical_ion(
    alchemical_ion_system, restraint, protocol, target, alchemical_ion_system_psores
):
    # Create an equilibration protocol with backbone restraints.
    protocol = protocol(restraint=restraint)

    # Create the process object.
    process = BSS.Process.Amber(
        alchemical_ion_system,
        protocol,
        name="test",
        reference_system=alchemical_ion_system_psores,
    )

    # Check that the correct restraint mask is in the config.
    config = " ".join(process.getConfig())
    assert target in config
    # Check is the reference file is passed to the cmd
    assert "-ref test_ref.rst7" in process.getArgString()

    # Test if the original coordinate is correct
    original = BSS.IO.readMolecules(
        [f"{process.workDir()}/test.prm7", f"{process.workDir()}/test.rst7"]
    )
    original_crd = original.getMolecule(0).coordinates()[0]
    assert str(original_crd) == "(17.9138 A, 16.0981 A, 20.5786 A)"
    # Test if the reference coordinate is passed
    ref = BSS.IO.readMolecules(
        [f"{process.workDir()}/test.prm7", f"{process.workDir()}/test_ref.rst7"]
    )
    ref_crd = ref.getMolecule(0).coordinates()[0]
    assert str(ref_crd) == "(0.0000e+00 A, 0.0000e+00 A, 0.0000e+00 A)"
