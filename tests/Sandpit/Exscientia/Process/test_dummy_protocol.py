import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Process._process import Process
from tests.conftest import root_fp
from tests.Sandpit.Exscientia.conftest import has_amber, has_gromacs


@pytest.fixture(scope="session")
def system():
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules(
        [f"{root_fp}/input/ala.top", f"{root_fp}/input/ala.crd"]
    )


@pytest.mark.skipif(
    has_amber is False or has_gromacs is False,
    reason="Requires AMBER and GROMACS to be installed.",
)
@pytest.mark.parametrize("process", ["Gromacs", "Amber"])
def test_process(system, process):
    protocol = BSS.Protocol.Dummy()
    process_cls = getattr(BSS.Process, process)
    process = process_cls(system, protocol)
    assert isinstance(process, Process)
