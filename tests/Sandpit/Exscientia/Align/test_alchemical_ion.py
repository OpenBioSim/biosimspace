import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align._alch_ion import (
    _get_protein_com_idx,
    _mark_alchemical_ion,
)
from BioSimSpace.Sandpit.Exscientia._SireWrappers import Atom
from tests.conftest import has_gromacs, root_fp


@pytest.fixture
def system():
    mol = BSS.IO.readMolecules(
        [f"{root_fp}/input/ala.top", f"{root_fp}/input/ala.crd"]
    )[0]
    system = BSS.Solvent.tip3p(mol, ion_conc=0.15, shell=2 * BSS.Units.Length.nanometer)
    return system


@pytest.fixture
def alchemical_ion_system(system):
    ion = system[-1]
    ion = _mark_alchemical_ion(ion)
    system.updateMolecules(ion)
    return system


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
@pytest.mark.parametrize(
    "input_system,isalchem", [("system", False), ("alchemical_ion_system", True)]
)
def test_isAlchemicalIon(input_system, isalchem, request):
    system = request.getfixturevalue(input_system)
    assert system[-1].isAlchemicalIon() is isalchem


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
@pytest.mark.parametrize(
    "input_system,isalchem", [("system", None), ("alchemical_ion_system", True)]
)
def test_getAlchemicalIon(input_system, isalchem, request):
    system = request.getfixturevalue(input_system)
    ion = system.getAlchemicalIon()
    if isalchem is None:
        assert ion is None
    else:
        assert isinstance(ion, Atom)


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_getAlchemicalIonIdx(alchemical_ion_system):
    index = alchemical_ion_system.getAlchemicalIonIdx()
    assert index == 2053


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_get_protein_com_idx(alchemical_ion_system):
    index = _get_protein_com_idx(alchemical_ion_system)
    assert index == 8
