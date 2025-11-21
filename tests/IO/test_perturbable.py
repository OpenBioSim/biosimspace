import pytest

import BioSimSpace as BSS

from tests.conftest import url


@pytest.fixture(scope="module")
def perturbable_molecule():
    """Re-use the same perturbable system for each test."""
    return BSS.IO.readPerturbableSystem(
        f"{url}/perturbable_system0.prm7",
        f"{url}/perturbable_system0.rst7",
        f"{url}/perturbable_system1.prm7",
        f"{url}/perturbable_system1.rst7",
    ).getMolecules()[0]


def test_connectivity(perturbable_molecule):
    """
    Make sure there is a single connectivity property when the end states
    have the same bonding.
    """

    assert perturbable_molecule._sire_object.has_property("connectivity")
    assert not perturbable_molecule._sire_object.has_property("connectivity0")
    assert not perturbable_molecule._sire_object.has_property("connectivity1")
