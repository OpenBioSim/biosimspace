import random

import pytest

import BioSimSpace.Types as Types

# Create a list of the Type classes.
types = []
for var in dir(Types):
    if var[0] == "_":
        continue
    elif var[0].upper() == var[0]:
        # Exclude non unit based types.
        if var != "Vector" and var != "Coordinate" and var != "GeneralUnit":
            types.append(getattr(Types, var))


def insert_space(string, index):
    return string[:index] + " " + string[index:]


def lower_case(string, index):
    x = len(string)
    if index == 0:
        return string[0].lower() + string[1:]
    elif index == x - 1:
        return string[: x - 1] + string[x - 1].lower()
    else:
        return string[:index] + string[index].lower() + string[index + 1 :]


@pytest.mark.parametrize("Type", types)
def test_supported_units(Type):
    """Test that it's possible to instatiate objects for all supported units."""

    # Store a list of the supported units.
    units = list(Type._supported_units.keys())

    # Loop over all units.
    for unit in units:
        # Attempt to instantiate an object.
        my_type = Type(1.0, unit)

        # Swap the case of random characters.
        for x in range(0, 5):
            unit = lower_case(unit, random.randint(0, len(unit) - 1))

        # Add random whitespace.
        for x in range(0, 5):
            unit = insert_space(unit, random.randint(0, len(unit) - 1))

        # Attempt to instantiate an object.
        my_type = Type(1.0, unit)


@pytest.mark.parametrize("Type", types)
def test_abbreviations(Type):
    """Test that all unit abbreviations work."""

    # Store a list of the unit abbreviations.
    abbreviations = list(Type._abbreviations.keys())

    # Loop over all abbreviations.
    for abbrev in abbreviations:
        # Attempt to instantiate an object.
        my_type = Type(1.0, abbrev)

        # Swap the case of random characters.
        for x in range(0, 5):
            abbrev = lower_case(abbrev, random.randint(0, len(abbrev) - 1))

        # Add random whitespace.
        for x in range(0, 5):
            unit = insert_space(abbrev, random.randint(0, len(abbrev) - 1))

        # Attempt to instantiate an object.
        my_type = Type(1.0, abbrev)


@pytest.mark.parametrize("Type", types)
def test_round_trip(Type):
    """Test that round-trip unit conversion works."""

    # Store a list of the supported units.
    units = list(Type._supported_units.keys())

    # Loop over all units.
    for index, unit in enumerate(units):
        # Create the initial object.
        my_type = Type(1.0, unit)

        # Convert to the next unit.
        for x in range(1, len(units) + 1):
            my_type = my_type._convert_to(units[index - x])

        # Make sure the unit and value are correct.
        assert my_type.value() == pytest.approx(1.0)
        assert unit == my_type.unit()


@pytest.mark.parametrize("Type", types)
def test_from_sire_unit(Type):
    """Test that types can be instantiated from their Sire unit equivalents."""

    # Store a list of the supported Sire units.
    sire_units = Type._supported_units.values()

    # Loop over all units.
    for sire_unit in Type._supported_units.values():
        # Try to instantiate the type from the Sire unit.
        my_type = Type(sire_unit)


@pytest.mark.parametrize("Type", types)
def test_container_mul(Type):
    """Test that types can applied to containers."""

    # Attempt to instantiate an object.
    my_type = Type(1.0, Type._default_unit)

    # Apply the type to a list.
    container = [1, 2, 3] * my_type

    # Make sure the container is a list.
    assert isinstance(container, list)

    # Assert we have a container of the original type.
    assert all(isinstance(x, Type) for x in container)

    # Apply the type to a tuple.
    container = (1, 2, 3) * my_type

    # Make sure the container is a tuple.
    assert isinstance(container, tuple)

    # Assert we have a container of the original type.
    assert all(isinstance(x, Type) for x in container)


@pytest.mark.parametrize("Type", types)
def test_operators(Type):
    """Test equality/inequality operators."""

    t0 = Type(1.0, Type._default_unit)
    t1 = Type(2.0, Type._default_unit)
    t2 = Type(3.0, Type._default_unit)

    assert t0 == t0
    assert t0 <= t0
    assert t0 >= t0
    assert t0 != t1
    assert t0 < t1
    assert t1 > t0
    assert t1 < t2
    assert t2 > t1
