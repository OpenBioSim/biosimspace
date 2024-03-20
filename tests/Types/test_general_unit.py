import pytest

import BioSimSpace.Types as Types
import BioSimSpace.Units as Units

import sire as sr


@pytest.mark.parametrize(
    "string, dimensions",
    [
        (
            "kilo Cal oriEs per Mole / angstrom **2",
            tuple(sr.u("kcal_per_mol / angstrom**2").dimensions()),
        ),
        ("k Cal_per  _mOl / nm^2", tuple(sr.u("kcal_per_mol / nm**2").dimensions())),
        (
            "kj p  eR  moles / pico METERs2",
            tuple(sr.u("kJ_per_mol / pm**2").dimensions()),
        ),
        (
            "coul oMbs / secs * ATm os phereS",
            tuple(sr.u("coulombs / second / atm").dimensions()),
        ),
        ("pm**3 * rads * de grEE", tuple(sr.u("pm**3 * rad * degree").dimensions())),
    ],
)
def test_supported_units(string, dimensions):
    """Test that we can create GeneralUnit objects with the correct dimensions
    by evaluating strings as unit based algebraic expressions.
    """

    # Try to create the GeneralUnit from the string.
    general_unit = Types._GeneralUnit(string)

    # Assert that the dimensions match.
    assert general_unit.dimensions() == dimensions


@pytest.mark.parametrize(
    "string, matching_type",
    [
        ("radian * degree**2 / radian^2", Types.Angle),
        ("angstrom**3 / nanometer", Types.Area),
        ("coulombs * angstrom**-2 * nanometer**2", Types.Charge),
        ("(kcal_per_mol / angstrom**2) * nanometer**2", Types.Energy),
        ("angstrom**3 * nanometer^-1 / picometer", Types.Length),
        ("bar * kJ_per_mol**2 / (kcal_per_mol * kJ_per_mol)", Types.Pressure),
        ("coulomb * kelvin^-3 * celsius**2 * kelvin^2 / e_charge", Types.Temperature),
        ("nanoseconds^3 * kelvin^-3 * celsius**3 / milliseconds**2", Types.Time),
        ("angstroms cubed * atm^-3 * bar**3", Types.Volume),
    ],
)
def test_type_conversion(string, matching_type):
    """Test that GeneralUnit objects can be converted to a type with matching
    dimensions.
    """

    # Try to create the GeneralUnit from the string.
    general_unit = Types._GeneralUnit(string)

    # Assert that the types match.
    assert type(general_unit) is matching_type


@pytest.mark.parametrize(
    "string, default_unit",
    [
        ("degree", Units.Angle.radian),
        ("meters2", Units.Area.angstrom2),
        ("coulombs", Units.Charge.electron_charge),
        ("kJ_per_mol", Units.Energy.kcal_per_mol),
        ("nanometer", Units.Length.angstrom),
        ("bar", Units.Pressure.atm),
        ("fahrenheit", Units.Temperature.kelvin),
        ("days", Units.Time.nanosecond),
        ("picometers**3", Units.Volume.angstrom3),
    ],
)
def test_default_conversion(string, default_unit):
    """Test that GeneralUnit objects are always converted to the default
    unit for that type.
    """

    # Try to create the GeneralUnit from the string.
    general_unit = Types._GeneralUnit(string)

    # Assert that units match.
    assert general_unit.unit() == default_unit.unit()


@pytest.mark.parametrize(
    "unit_type",
    [
        Units.Angle.radian,
        Units.Area.angstrom2,
        Units.Charge.electron_charge,
        Units.Energy.kcal_per_mol,
        Units.Length.angstrom,
        Units.Pressure.atm,
        Units.Temperature.kelvin,
        Units.Time.nanosecond,
        Units.Volume.angstrom3,
    ],
)
def test_pos_pow(unit_type):
    """Test that unit-based types can be raised to positive powers."""

    # Store the dimensions associated with the original type.
    old_dimensions = unit_type.dimensions()

    # Square the unit-based type.
    unit_type = unit_type**2

    # Store the new dimensions.
    new_dimensions = unit_type.dimensions()

    # Each dimension entry should be twice the old value.
    for d0, d1 in zip(old_dimensions, new_dimensions):
        assert d1 == 2 * d0


@pytest.mark.parametrize(
    "unit_type",
    [
        Units.Angle.radian,
        Units.Area.angstrom2,
        Units.Charge.electron_charge,
        Units.Energy.kcal_per_mol,
        Units.Length.angstrom,
        Units.Pressure.atm,
        Units.Temperature.kelvin,
        Units.Time.nanosecond,
        Units.Volume.angstrom3,
    ],
)
def test_neg_pow(unit_type):
    """Test that unit-based types can be raised to negative powers."""

    # Store the dimensions associated with the original type.
    old_dimensions = unit_type.dimensions()

    # Invert the unit-based type.
    unit_type = unit_type**-1

    # Store the new dimensions.
    new_dimensions = unit_type.dimensions()

    # Each dimension entry should be the inverse of the old value.
    for d0, d1 in zip(old_dimensions, new_dimensions):
        assert d1 == -d0


def test_frac_pow():
    """Test that unit-based types can be raised to fractional powers."""

    # Create a base unit type.
    unit_type = 2 * Units.Length.angstrom

    # Store the original value and dimensions.
    value = unit_type.value()
    dimensions = unit_type.dimensions()

    # Square the type.
    unit_type = unit_type**2

    # Assert that we can't take the cube root.
    with pytest.raises(ValueError):
        unit_type = unit_type ** (1 / 3)

    # Now take the square root.
    unit_type = unit_type ** (1 / 2)

    # The value should be the same.
    assert unit_type.value() == value

    # The dimensions should be the same.
    assert unit_type.dimensions() == dimensions

    # Cube the type.
    unit_type = unit_type**3

    # Assert that we can't take the square root.
    with pytest.raises(ValueError):
        unit_type = unit_type ** (1 / 2)

    # Now take the cube root.
    unit_type = unit_type ** (1 / 3)

    # The value should be the same.
    assert unit_type.value() == value

    # The dimensions should be the same.
    assert unit_type.dimensions() == dimensions

    # Square the type again.
    unit_type = unit_type**2

    # Now take the negative square root.
    unit_type = unit_type ** (-1 / 2)

    # The value should be inverted.
    assert unit_type.value() == 1 / value

    # The dimensions should be negated.
    assert unit_type.dimensions() == tuple(-d for d in dimensions)


@pytest.mark.parametrize(
    "string",
    [
        "degree",
        "meters2",
        "coulombs",
        "kJ_per_mol",
        "nanometer",
        "bar",
        "fahrenheit",
        "days",
        "picometers**3",
    ],
)
def test_dimensionless(string):
    """Test that GeneralUnit objects convert to dimensionless float values
    when divided by themself.
    """

    # Try to create the GeneralUnit from the string.
    general_unit = Types._GeneralUnit(string)

    # Check that we get back a float when divided by itself.
    assert isinstance(general_unit / general_unit, float)


def test_dimensionless_value():
    """Check that conversion to a dimensionless unit preserves the value
    of the unit conversion.
    """

    value = (Units.Energy.kcal_per_mol / Units.Length.angstrom**2) / (
        Units.Energy.kj_per_mol / Units.Length.nanometer**2
    )

    assert value == pytest.approx(418.4)


def test_value_and_unit():
    """
    Regression test to make sure that a general unit with a value and unit can
    be parsed correctly.
    """
    general_unit = Types._GeneralUnit(2, "kcal per mol / angstrom**2")
