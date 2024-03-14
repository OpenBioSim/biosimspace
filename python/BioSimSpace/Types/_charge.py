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

"""A charge type."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Charge"]

from sire.legacy import Units as _SireUnits

from ._type import Type as _Type


class Charge(_Type):
    """A charge type."""

    # A list of the supported Sire unit names.
    _sire_units = ["e_charge", "coulomb"]

    # Dictionary of allowed units.
    _supported_units = {
        "ELECTRON CHARGE": _SireUnits.e_charge,
        "COULOMB": _SireUnits.coulomb,
    }

    # Map unit abbreviations to the full name.
    _abbreviations = {"E": "ELECTRON CHARGE", "C": "COULOMB"}

    # Print formatting.
    _print_format = {"ELECTRON CHARGE": "|e|", "COULOMB": "C"}

    # Documentation strings.
    _doc_strings = {
        "ELECTRON CHARGE": "A charge in electron charge.",
        "COULOMB": "A charge in Coulomb.",
    }

    # Null type unit for avoiding issue printing configargparse help.
    _default_unit = "ELECTRON CHARGE"

    # The dimension mask:
    #     Angle, Charge, Length, Mass, Quantity, Temperature, Time
    _dimensions = (0, 1, 0, 0, 0, 0, 0)

    def __init__(self, *args):
        """
        Constructor.

        ``*args`` can be a value and unit, or a string representation
        of the charge, e.g. "0.2 e".

        Parameters
        ----------

        value : float
            The value.

        unit : str
            The unit.

        string : str
            A string representation of the charge.

        Examples
        --------

        Create an object representing a charge of 27.8 electron charge
        print the charge in Coulomb.

        >>> import BioSimSpace as BSS
        >>> charge = BSS.Types.Charge(27.8, "e")
        >>> print(charge.coulomb())

        The same as above, except passing a string representation of the
        charge to the constructor.

        >>> import BioSimSpace as BSS
        >>> charge = BSS.Types.Charge("3.1 atm")
        >>> print(charge.coulomb())

        The string matching is extremeley flexible, so all of the following
        would be valid arguments: "27.8 e", "27.8 electron charge",
        "2.78e1 e".
        """

        # Call the base class constructor.
        super().__init__(*args)

    def electron_charge(self):
        """
        Return the charge in electron charge.

        Returns
        -------

        charge : :class:`Charge <BioSimSpace.Types.Charge>`
            The charge in electron charge.
        """
        return Charge(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.e_charge),
            "ELECTRON CHARGE",
        )

    def coulomb(self):
        """
        Return the charge in Coulomb.

        Returns
        -------

        charge : :class:`Charge <BioSimSpace.Types.Charge>`
            The charge in Coulomb.
        """
        return Charge(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.coulomb),
            "COULOMB",
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

        charge : :class:`Charge <BioSimSpace.Types.Charge>`
            The charge in the default unit of electron charge.
        """
        if mag is None:
            return self.electron_charge()
        else:
            return Charge(mag, "ELECTRON CHARGE")

    def _convert_to(self, unit):
        """
        Return the charge in a different unit.

        Parameters
        ----------

        unit : str
            The unit to convert to.

        Returns
        -------

        charge : :class:`Charge <BioSimSpace.Types.Charge>`
            The charge in the specified unit.
        """
        if unit == "ELECTRON CHARGE":
            return self.electron_charge()
        elif unit == "COULOMB":
            return self.coulomb()
        else:
            raise ValueError(
                "Supported units are: '%s'" % list(self._supported_units.keys())
            )

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Replace all instances of ELECTRON with E.
        unit = unit.replace("ELECTRON", "E")

        # Replace all instances of ELECT with E.
        unit = unit.replace("ELECT", "E")

        # Replace all instances of ELEC with E.
        unit = unit.replace("ELEC", "E")

        # Strip all instances of CHARGE.
        unit = unit.replace("CHARGE", "")

        # Strip all instances of -.
        unit = unit.replace("-", "")

        # Strip all instances of S.
        unit = unit.replace("S", "")

        # Replace all instance of COULOMB with C.
        unit = unit.replace("COULOMB", "C")

        # Replace all instance of COUL with C.
        unit = unit.replace("COUL", "C")

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

        unit = unit.replace("coulombs", "coulomb")
        unit = unit.replace("electroncharge", "e_charge")
        unit = unit.replace("echarge", "e_charge")

        # Convert powers. (Just 2nd and third for now.)
        unit = unit.replace("coloumb2", "(coulomb*coulomb)")
        unit = unit.replace("coloumb3", "(coulomb*coulomb*coulomb)")
        unit = unit.replace("coloumb-1", "(1/coulomb)")
        unit = unit.replace("coloumb-2", "(1/(coulomb*coulomb))")
        unit = unit.replace("coloumb3", "(1/(coulomb*coulomb*coulomb))")
        unit = unit.replace("e_charge2", "(e_charge*e_charge)")
        unit = unit.replace("e_charge3", "(e_charge*e_charge*e_charge)")
        unit = unit.replace("e_charge-1", "(1/e_charge)")
        unit = unit.replace("e_charge-2", "(1/(e_charge*e_charge))")
        unit = unit.replace("e_charge-3", "(1/(e_charge*e_charge*e_charge))")

        return unit
