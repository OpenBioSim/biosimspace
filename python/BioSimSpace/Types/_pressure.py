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

"""A pressure type."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["Pressure"]

from sire.legacy import Units as _SireUnits

from ._type import Type as _Type


class Pressure(_Type):
    """A pressure type."""

    # A list of the supported Sire unit names.
    _sire_units = ["atm", "bar"]

    # Dictionary of allowed units.
    _supported_units = {"ATMOSPHERE": _SireUnits.atm, "BAR": _SireUnits.bar}

    # Map unit abbreviations to the full name.
    _abbreviations = {"ATM": "ATMOSPHERE", "BAR": "BAR"}

    # Print formatting.
    _print_format = {"ATMOSPHERE": "atm", "BAR": "bar"}

    # Documentation strings.
    _doc_strings = {
        "ATMOSPHERE": "A pressure in atmosphere.",
        "BAR": "A pressure in bar.",
    }

    # Null type unit for avoiding issue printing configargparse help.
    _default_unit = "ATMOSPHERE"

    # The dimension mask:
    #     Angle, Charge, Length, Mass, Quantity, Temperature, Time
    _dimensions = (0, 0, -1, 1, 0, 0, -2)

    def __init__(self, *args):
        """
        Constructor.

        ``*args`` can be a value and unit, or a string representation
        of the pressure, e.g. "1.38 bar".

        Parameters
        ----------

        value : float
            The value.

        unit : str
            The unit.

        string : str
            A string representation of the pressure.

        Examples
        --------

        Create an object representing a pressure of 3.1 atmosphere then
        print the pressure in bar.

        >>> import BioSimSpace as BSS
        >>> pressure = BSS.Types.Pressure(3.1, "atm")
        >>> print(pressure.bar())

        The same as above, except passing a string representation of the
        pressure to the constructor.

        >>> import BioSimSpace as BSS
        >>> pressure = BSS.Types.Pressure("3.1 atm")
        >>> print(pressure.nanoseconds())

        The string matching is extremeley flexible, so all of the following
        would be valid arguments: "3.1 atm", "3.1 atmosphere", "0.31e1 atm".
        """

        # Call the base class constructor.
        super().__init__(*args)

    def atm(self):
        """
        Return the atmospheric pressure.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure in atomospheres.
        """
        return Pressure(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.atm),
            "ATMOSPHERE",
        )

    def bar(self):
        """
        Return the pressure in bar.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure in bar.
        """
        return Pressure(
            (self._value * self._supported_units[self._unit]).to(_SireUnits.bar), "BAR"
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

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure in the default unit of atomospheres.
        """
        if mag is None:
            return self.atm()
        else:
            return Pressure(mag, "ATMOSPHERE")

    def _convert_to(self, unit):
        """
        Return the pressure in a different unit.

        Parameters
        ----------

        unit : str
            The unit to convert to.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure in the specified unit.
        """
        if unit == "ATMOSPHERE":
            return self.atm()
        elif unit == "BAR":
            return self.bar()
        else:
            raise ValueError(
                "Supported units are: '%s'" % list(self._supported_units.keys())
            )

    def _validate_unit(self, unit):
        """Validate that the unit are supported."""

        # Strip whitespace and convert to upper case.
        unit = unit.replace(" ", "").upper()

        # Strip all instances of PRESSURE.
        unit = unit.replace("PRESSURE", "")
        unit = unit.replace("PRESS", "")
        unit = unit.replace("PRES", "")

        # Replace all instance of ATMOSPHERIC/ATMOSPHERE with ATM.
        unit = unit.replace("ATMOSPHERIC", "ATM")
        unit = unit.replace("ATMOSPHERE", "ATM")

        # Strip all "S" characters.
        unit = unit.replace("S", "")

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

        unit = unit.replace("pressure", "")
        unit = unit.replace("press", "")
        unit = unit.replace("pres", "")
        unit = unit.replace("atmospheric", "atm")
        unit = unit.replace("atmospheres", "atm")
        unit = unit.replace("atmosphere", "atm")

        # Convert powers. (Just 2nd and third for now.)
        unit = unit.replace("bar2", "(bar*bar)")
        unit = unit.replace("bar3", "(bar*bar*bar)")
        unit = unit.replace("bar-1", "(1/bar)")
        unit = unit.replace("bar-2", "(1/(bar*bar))")
        unit = unit.replace("bar3", "(1/(bar*bar*bar))")
        unit = unit.replace("atm2", "(atm*atm)")
        unit = unit.replace("atm3", "(atm*atm*atm)")
        unit = unit.replace("atm-1", "(1/atm)")
        unit = unit.replace("atm-2", "(1/(atm*atm))")
        unit = unit.replace("atm-3", "(1/(atm*atm*atm))")

        return unit
