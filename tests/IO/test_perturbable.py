import pytest

import BioSimSpace as BSS

from tests.conftest import url


def test_connectivity(perturbable_molecule):
    """
    Make sure there is a single connectivity property when the end states
    have the same bonding.
    """

    assert perturbable_molecule._sire_object.has_property("connectivity")
    assert not perturbable_molecule._sire_object.has_property("connectivity0")
    assert not perturbable_molecule._sire_object.has_property("connectivity1")
