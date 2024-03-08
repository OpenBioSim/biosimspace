######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2024
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

"""An energy type."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Energy"]

from sire.legacy import Units as _SireUnits

from ._type import Type as _Type


class Energy(_Type):
    """An energy type."""

    # A list of the supported Sire unit names.
    _sire_units = ["kcal_per_mol", "kJ_per_mol"]

    # Dictionary of allowed units.
    _supported_units = {
        "KILO CALORIES PER MOL": _SireUnits.kcal_per_mol,
        "KILO JOULES PER MOL": _SireUnits.kJ_per_mol,
        "KT": 2.479 * _SireUnits.kJ_per_mol,
    }

    # Map unit abbreviations to the full name.
    _abbreviations = {
        "KCAL/MOL": "KILO CALORIES PER MOL",
        "KJ/MOL": "KILO JOULES PER MOL",
        "KT": "KT",
    }

    # Print formatting.
    _print_format = {
        "KILO CALORIES PER MOL": "kcal/mol",
        "KILO JOULES PER MOL": "kJ/mol",
        "KT": "KT",
    }

    # Documentation strings.
    _doc_strings = {
        "KILO CALORIES PER MOL": "An energy in kcal per mol.",
        "KILO JOULES PER MOL": "An energy in kJ per mol.",
        "KT": "An energy in KT.",
    }

    # Null type unit for avoiding issue printing configargparse help.
    _default_unit = "KILO CALORIES PER MOL"

    # The dimension mask.
    _dimensions = tuple(list(_supported_units.values())[0].dimensions())

    def __init__(self, *args):
        """
        Constructor.

        ``*args`` can be a value and unit, or a string representation
        of the energy, e.g. "78.4 kcal/mol".

        Parameters
        ----------

        value : float
            The value.

        unit : str
            The unit.

        string : str
            A string representation of the energy.

        Examples
        --------

        Create an object representing an energy of -1038 kilo calories per
        mol and print the energy in kilo joules per mol.

        >>> import BioSimSpace as BSS
        >>> energy = BSS.Types.Energy(-1038, "kcal/mol")
        >>> print(energy.kj_per_mol())

        The same as above, except passing a string representation of the
        energy to the constructor.

        >>> import BioSimSpace as BSS
        >>> energy = BSS.Types.Energy("-1038 kcal/mol")
        >>> print(energy.kj_per_mol())

        The string matching is extremeley flexible, so all of the following
        would be valid arguments: "-1038 kcal/mol", "-1.038e3 kcal/mol",
        "-1038 kilo cal per mol".
        """

        # Call the base class constructor.
        super().__init__(*args)

    def kcal_per_mol(self):
        """
        Return the energy in kcal per mol.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The energy in kcal per mol.
        """
        return Energy(
            (self._value * self._supported_units[self._unit]).to(
                _SireUnits.kcal_per_mol
            ),
            "KILO CALORIES PER MOL",
        )

    def kj_per_mol(self):
        """
        Return the energy in kJ per mol.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The energy in kJ per mol.
        """
        return Energy(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.kJ_per_mol),
            "KILO JOULES PER MOL",
        )

    def kt(self):
        """
        Return the energy in KT.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The energy in KT.
        """
        return Energy(
            (self._value * self._supported_units[self._unit]).to(
                2.479 * _SireUnits.kJ_per_mol
            ),
            "KT",
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

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The energy in the default unit of kcal per mol.
        """
        if mag is None:
            return self.kcal_per_mol()
        else:
            return Energy(mag, "KILO CALORIES PER MOL")

    def _convert_to(self, unit):
        """
        Return the energy in a different unit.

        Parameters
        ----------

        unit : str
            The unit to convert to.

        Returns
        -------

        energy : :class:`Energy <BioSimSpace.Types.Energy>`
            The energy in the specified unit.
        """
        if unit == "KILO CALORIES PER MOL":
            return self.kcal_per_mol()
        elif unit == "KILO JOULES PER MOL":
            return self.kj_per_mol()
        elif unit == "KT":
            return self.kt()
        else:
            raise ValueError(
                "Supported units are: '%s'" % list(self._supported_units.keys())
            )

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Replace all instances of KILO with K.
        unit = unit.replace("KILO", "K")

        # Replace all instances of PER with the / character.
        unit = unit.replace("PER", "/")

        # Replace all instance of MOLE with MOL.
        unit = unit.replace("MOLE", "MOL")

        # Replace all instances of CALORIES with CAL.
        unit = unit.replace("CALORIES", "CAL")

        # Replace all instances of JOULES with J.
        unit = unit.replace("JOULES", "J")

        # Check that the unit is supported.
        if unit in self._abbreviations:
            return self._abbreviations[unit]
        else:
            raise ValueError(
                "Supported units are: '%s'" % list(self._supported_units.keys())
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

        unit = unit.replace("mole", "mol")
        unit = unit.replace("mols", "mol")
        unit = unit.replace("calories", "cal")
        unit = unit.replace("joules", "J")
        unit = unit.replace("kilocal", "kcal")
        unit = unit.replace("kiloJ", "kJ")
        unit = unit.replace("kj", "kJ")
        unit = unit.replace("kcals", "kcal")
        unit = unit.replace("kJs", "kJ")
        unit = unit.replace("kcalpermol", "kcal_per_mol")
        unit = unit.replace("kJpermol", "kJ_per_mol")

        # Convert powers.
        unit = unit.replace("kcal_per_mol2", "(kcal_per_mol*kcal_per_mol)")
        unit = unit.replace("kcal_per_mol3", "(kcal_per_mol*kcal_per_mol*kcal_per_mol)")
        unit = unit.replace("kcal_per_mol-1", "(1/kcal_per_mol)")
        unit = unit.replace("kcal_per_mol-2", "(1/(kcal_per_mol*kcal_per_mol))")
        unit = unit.replace(
            "kcal_per_mol-3", "(1/(kcal_per_mol*kcal_per_mol*kcal_per_mol))"
        )
        unit = unit.replace("kJ_per_mol2", "(kJ_per_mol*kJ_per_mol)")
        unit = unit.replace("kJ_per_mol3", "(kJ_per_mol*kJ_per_mol*kJ_per_mol)")
        unit = unit.replace("kJ_per_mol-1", "(1/kJ_per_mol)")
        unit = unit.replace("kJ_per_mol-2", "(1/(kJ_per_mol*kJ_per_mol))")
        unit = unit.replace("kJ_per_mol-3", "(1/(kJ_per_mol*kJ_per_mol*kJ_per_mol))")

        return unit
