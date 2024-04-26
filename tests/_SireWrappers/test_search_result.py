import BioSimSpace as BSS

import pytest

# Store the tutorial URL.
url = BSS.tutorialUrl()


@pytest.fixture(scope="session")
def molecule(system):
    """Re-use the same molecule for each test."""
    return system[0]


@pytest.fixture(scope="session")
def residue(molecule):
    """Re-use the same residue for each test."""
    return molecule.getResidues()[0]


# Parameterise the function with a set of searchable objects.
@pytest.mark.parametrize("fixture", ["system", "molecule", "residue"])
# Ensure that valid searches return a result.
def test_valid_search(fixture, request):
    obj = request.getfixturevalue(fixture)
    obj.search("element H")


# Parameterise the function with a set of searchable objects.
@pytest.mark.parametrize("fixture", ["system", "molecule", "residue"])
# Ensure that invalid searches raise a ValueError exception.
def test_invalid_search(fixture, request):
    obj = request.getfixturevalue(fixture)

    # Assert that a search for an invalid element raises a ValueError.
    with pytest.raises(ValueError):
        obj.search("element Z")

    # Assert that a search for a valid element using an invalid property
    # map raises a ValueError.
    with pytest.raises(ValueError):
        obj.search("element H", property_map={"element": "cheese"})
