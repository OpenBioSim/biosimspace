######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2025
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BioSimSpace is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""A mass type."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Mass"]

from sire.legacy import Units as _SireUnits

from ._type import Type as _Type


class Mass(_Type):
    """A mass type."""

    # A list of the supported Sire unit names.
    _sire_units = [
        "kilogram",
        "gram",
        "milligram",
        "microgram",
        "nanogram",
        "picogram",
        "femtogram",
    ]

    # Dictionary of allowed units.
    _supported_units = {
        "KILOGRAM": _SireUnits.kilogram,
        "GRAM": _SireUnits.gram,
        "MILLIGRAM": _SireUnits.milligram,
        "MICROGRAM": _SireUnits.microgram,
        "NANOGRAM": _SireUnits.nanogram,
        "PICOGRAM": _SireUnits.picogram,
        "FEMTOGRAM": _SireUnits.femtogram,
    }

    # Map unit abbreviations to the full name.
    _abbreviations = {
        "KG": "KILOGRAM",
        "G": "GRAM",
        "MG": "MILLIGRAM",
        "UG": "MICROGRAM",
        "NG": "NANOGRAM",
        "PG": "PICOGRAM",
        "FG": "FEMTOGRAM",
    }

    # Print format.
    _print_format = {
        "GRAM": "g",
        "MILLIGRAM": "mg",
        "MICROGRAM": "ug",
        "NANOGRAM": "ng",
        "PICOGRAM": "pg",
        "FEMTOGRAM": "fg",
    }

    # Documentation strings.
    _doc_strings = {
        "KILOGRAM": "A mass in kilogram.",
        "GRAM": "A mass in gram.",
        "MILLIGRAM": "A mass in milligram.",
        "MICROGRAM": "A mass in microgram.",
        "NANOGRAM": "A mass in nanogram.",
        "PICOGRAM": "A mass in picogram.",
        "FEMTOGRAM": "A mass in femtogram.",
    }

    # Null type unit for avoiding issue printing configargparse help.
    _default_unit = "GRAM"

    # The dimension mask.
    _dimensions = tuple(list(_supported_units.values())[0].dimensions())

    def __init__(self, *args):
        """
        Constructor.

        ``*args`` can be a value and unit, or a string representation
        of the length, e.g. "12 grams".

        Parameters
        ----------

        value : float
            The value.

        unit : str
            The unit.

        string : str
            A string representation of the mass.

        Examples
        --------

        Create an object representing a mass of 12 grams then
        print the mass in kilograms.

        >>> import BioSimSpace as BSS
        >>> length = BSS.Types.Length(12, "G")
        >>> print(length.kilograms())

        The same as above, except passing a string representation of the
        length to the constructor.

        >>> import BioSimSpace as BSS
        >>> length = BSS.Types.Length("12 G")
        >>> print(length.kilograms())

        The string matching is extremeley flexible, so all of the following
        would be valid arguments: "12 G", "12 grams", "1.2e1 grams".
        """

        # Call the base class constructor.
        super().__init__(*args)

    def __mul__(self, other):
        """Multiplication operator."""

        # Handle containers by converting each item in the container to
        # this type.
        if isinstance(other, list):
            return [self.__mul__(item) for item in other]
        if isinstance(other, tuple):
            return tuple([self.__mul__(item) for item in other])

        # Convert int to float.
        if type(other) is int:
            other = float(other)

        # Multiplication by float.
        if isinstance(other, float):
            mag = self._value * other
            return Mass(mag, self._unit)

        # Multiplication by another type.
        elif isinstance(other, _Type):
            from ._general_unit import GeneralUnit as _GeneralUnit

            return _GeneralUnit(self._to_sire_unit() * other._to_sire_unit())

        # Multiplication by a string.
        elif isinstance(other, str):
            try:
                mass = Mass(other)
                return self * length
            except:
                raise ValueError(
                    "Could not convert the string to a 'BioSimSpace.Mass' type."
                )
        else:
            raise TypeError(
                "unsupported operand type(s) for *: '%s' and '%s'"
                % (type(self), type(other))
            )

    def __rmul__(self, other):
        """Multiplication operator."""

        # Multiplication is commutative: a*b = b*a
        return self.__mul__(other)

    def kilogram(self):
        """
        Return the mass in kilograms.

        Returns
        -------

        mass : :class:`Mass <BioSimSpace.Types.Mass>`
            The mass in kilograms.
        """
        return Mass(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.kilogram),
            "KILOGRAM",
        )

    def gram(self):
        """
        Return the mass in grams.

        Returns
        -------

        mass : :class:`Mass <BioSimSpace.Types.Mass>`
            The mass in grams.
        """
        return Mass(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.gram),
            "GRAM",
        )

    def milligram(self):
        """
        Return the mass in milligrams.

        Returns
        -------

        mass : :class:`Mass <BioSimSpace.Types.Mass>`
            The mass in milligrams.
        """
        return Mass(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.milligram),
            "MILLIGRAM",
        )

    def microgram(self):
        """
        Return the mass in micrograms.

        Returns
        -------

        mass : :class:`Mass <BioSimSpace.Types.Mass>`
            The mass in micrograms.
        """
        return Mass(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.microgram),
            "MICROGRAM",
        )

    def nanogram(self):
        """
        Return the mass in nanograms.

        Returns
        -------

        mass : :class:`Mass <BioSimSpace.Types.Mass>`
            The mass in nanograms.
        """
        return Mass(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.nanogram),
            "NANOGRAM",
        )

    def picogram(self):
        """
        Return the mass in picograms.

        Returns
        -------

        mass : :class:`Mass <BioSimSpace.Types.Mass>`
            The mass in picograms.
        """
        return Mass(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.picogram),
            "PICOGRAM",
        )

    def femtogram(self):
        """
        Return the mass in femtograms.

        Returns
        -------

        mass : :class:`Mass <BioSimSpace.Types.Mass>`
            The mass in femtograms.
        """
        return Mass(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.femtogram),
            "FEMTOGRAM",
        )

    def _to_default_unit(self, mag=None):
        """
        Internal method to return an object of the same type in the default unit.

        Parameters
        ----------

        mag : float
           The value (optional).

        Returns
        -------

        length : :class:`Mass <BioSimSpace.Types.Mass>`
            The mass in the default unit of grams.
        """
        if mag is None:
            return self.gram()
        else:
            return Mass(mag, "GRAM")

    def _convert_to(self, unit):
        """
        Return the mass in a different unit.

        Parameters
        ----------

        unit : str
            The unit to convert to.

        Returns
        -------

        length : :class:`Mass <BioSimSpace.Types.Mass>`
            The mass in the specified unit.
        """
        if unit == "KILOGRAM":
            return self.kilogram()
        elif unit == "GRAM":
            return self.gram()
        elif unit == "MILLIGRAM":
            return self.milligram()
        elif unit == "MICROGRAM":
            return self.microgram()
        elif unit == "NANOGRAM":
            return self.nanogram()
        elif unit == "PICOGRAM":
            return self.picogram()
        elif unit == "FEMTOGRAM":
            return self.femtogram()
        else:
            raise ValueError(
                "Supported units are: '%s'" % list(self._supported_units.keys())
            )

    @classmethod
    def _validate_unit(cls, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Strip any "S" characters.
        unit = unit.replace("S", "")

        # Check that the unit is supported.
        if unit in cls._supported_units:
            return unit
        elif unit in cls._abbreviations:
            return cls._abbreviations[unit]
        else:
            raise ValueError(
                "Supported units are: '%s'" % list(cls._supported_units.keys())
            )

    @staticmethod
    def _to_sire_format(unit):
        """
        Reformat the unit string so it adheres to the Sire unit formatting.

        Parameters
        ----------

        unit : str
            A string representation of the unit.

        Returns
        -------

        sire_unit : str
            The unit string in Sire compatible format.
        """

        unit = unit.replace("kg", "kilogram")
        unit = unit.replace("ug", "microgram")
        unit = unit.replace("mg", "milligram")
        unit = unit.replace("ng", "nanogram")
        unit = unit.replace("pg", "picogram")
        unit = unit.replace("fg", "femtogram")
        unit = unit.replace("g", "gram")

        # Convert powers.
        unit = unit.replace("kilogram-1", "(1/kilogram)")
        unit = unit.replace("gram-1", "(1/gram)")
        unit = unit.replace("milligram-1", "(1/milligram)")
        unit = unit.replace("microgram-1", "(1/microgram)")
        unit = unit.replace("nanogram-1", "(1/nanogram)")
        unit = unit.replace("picogram-1", "(1/picogram)")
        unit = unit.replace("femtogram-1", "(1/femtogram)")

        return unit
