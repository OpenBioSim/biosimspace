import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Process._process import Process


@pytest.fixture(scope="session")
def system():
    """Re-use the same molecuar system for each test."""
    return BSS.IO.readMolecules(["test/input/ala.top", "test/input/ala.crd"])


@pytest.mark.parametrize("process", ["Gromacs", "Amber"])
def test_process(system, process):
    protocol = BSS.Protocol.Dummy()
    process_cls = getattr(BSS.Process, process)
    process = process_cls(system, protocol)
    assert isinstance(process, Process)
