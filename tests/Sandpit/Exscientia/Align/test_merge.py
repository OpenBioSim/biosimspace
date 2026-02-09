import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align._merge import _removeDummies
from tests.Sandpit.Exscientia.conftest import (
    get_energy,
    has_alchemlyb,
    has_amber,
    has_gromacs,
    has_openff,
)


@pytest.mark.skipif(
    has_amber is False
    or has_gromacs is False
    or has_openff is False
    or has_alchemlyb is False,
    reason="Requires that AMBER, GROMACS, OpenFF, and alchemlyb are installed.",
)
def test__removeDummies(
    benzene, pyrrole, mapping_benzene_pyrrole, merged_benzene_pyrrole
):
    benzene_merged = _removeDummies(merged_benzene_pyrrole, is_lambda1=False)
    pyrrole_merged = _removeDummies(merged_benzene_pyrrole, is_lambda1=True)
    energy0_unmerged = get_energy(benzene)
    energy1_unmerged = get_energy(pyrrole)
    energy0_merged = get_energy(benzene_merged)
    energy1_merged = get_energy(
        pyrrole_merged, coord_mol=pyrrole, coord_mol_mapping=mapping_benzene_pyrrole
    )
    assert energy0_unmerged == pytest.approx(energy0_merged)
    assert energy1_unmerged == pytest.approx(energy1_merged)
