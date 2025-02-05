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

"""A mixin for position restraints."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["_PositionRestraintMixin"]

from .. import Types as _Types
from .. import Units as _Units


class _PositionRestraintMixin:
    """A mixin for applying position restraints to other protocols."""

    # Supported restraint keywords.
    _restraints = ["backbone", "heavy", "all", "none"]

    def __init__(
        self,
        restraint=None,
        force_constant=10 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
    ):
        """Constructor.

        Parameters
        ----------

        restraint : str, [int]
            The type of restraint to perform. This should be one of the
            following options:
                "backbone"
                     Protein backbone atoms. The matching is done by a name
                     template, so is unreliable on conversion between
                     molecular file formats.
                "heavy"
                     All non-hydrogen atoms that aren't part of water
                     molecules or free ions.
                "all"
                     All atoms that aren't part of water molecules or free
                     ions.
            Alternatively, the user can pass a list of atom indices for
            more fine-grained control. If None, then no restraints are used.

        force_constant : :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`, float
            The force constant for the restraint potential. If a 'float' is
            passed, then default units of 'kcal_per_mol / angstrom**2' will
            be used.
        """
        # Set the restraint.
        if restraint is not None:
            self.setRestraint(restraint)
        else:
            self._restraint = None

        # Set the force constant.
        self.setForceConstant(force_constant)

    def _get_parm(self):
        """Return a string representation of the parameters."""
        if self._restraint is not None:
            return (
                f"restraint='{self._restraint}', "
                f"force_constant={self._force_constant.value()} kcal_per_mol/angstrom**2"
            )
        else:
            return (
                f"restraint={self._restraint}, "
                f"force_constant={self._force_constant.value()} kcal_per_mol/angstrom**2"
            )

    def __str__(self):
        print(self._get_parm())
        """Return a human readable string representation of the object."""
        return f"<BioSimSpace.Protocol._PositionRestraintMixin: {self._get_parm()}>"

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return f"BioSimSpace.Protocol._PositionRestraintMixin({self._get_parm()})"

    def __eq__(self, other):
        """Equality operator."""

        if not isinstance(other, _PositionRestraintMixin):
            return False

        return (
            self._restraint == other._restraint
            and self._force_constant == other._force_constant
        )

    def getRestraint(self):
        """
        Return the type of restraint.

        Returns
        -------

        restraint : str, [int]
            The type of restraint, either a keyword or a list of atom indices.
        """
        return self._restraint

    def setRestraint(self, restraint):
        """Set the type of restraint.

        Parameters
        ----------

        restraint : str, [int]
            The type of restraint to perform. This should be one of the
            following options:
                "backbone"
                     Protein backbone atoms. The matching is done by a name
                     template, so is unreliable on conversion between
                     molecular file formats.
                "heavy"
                     All non-hydrogen atoms that aren't part of water
                     molecules or free ions.
                "all"
                     All atoms that aren't part of water molecules or free
                     ions.
            Alternatively, the user can pass a list of atom indices for
            more fine-grained control.
        """

        if isinstance(restraint, str):
            # Convert to lower case and strip whitespace.
            restraint = restraint.lower().replace(" ", "")
            if restraint not in self._restraints:
                raise ValueError(f"'restraint' must be one of: {self._restraints}")
            # Set to NoneType if equal to "none", since this makes checking
            # whether a restraint is set elsewhere much easier.
            if restraint == "none":
                restraint = None

        elif isinstance(restraint, (list, tuple)):
            if not all(type(x) is int for x in restraint):
                raise ValueError("'restraint' must be a list of 'int' types!")
            # Create a set to sort and ensure no duplicates, then convert back to a list.
            restraint = list(set(restraint))
            restraint.sort()

        else:
            raise TypeError(
                "'restraint' must be of type 'str', or a list of 'int' types."
            )

        self._restraint = restraint

    def getForceConstant(self):
        """Return the force constant for the restraint.

        Returns
        -------

        force_constant : float, str, class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`
            The force constant for the restraint, in units of
            kcal_per_mol/angstrom**2.
        """
        return self._force_constant

    def setForceConstant(self, force_constant):
        """
        Set the force constant for the restraint.

        Parameters
        ----------

        force_constant : int, float, str, :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`, float
            The force constant for the restraint, in units of
            kcal_per_mol/angstrom**2.
        """

        # Convert int to float.
        if type(force_constant) is int:
            force_constant = float(force_constant)

        if isinstance(force_constant, float):
            # Use default units.
            force_constant *= _Units.Energy.kcal_per_mol / _Units.Area.angstrom2

        else:
            if isinstance(force_constant, str):
                try:
                    force_constant = _Types._GeneralUnit(force_constant)
                except Exception:
                    raise ValueError(
                        "Unable to parse 'force_constant' string."
                    ) from None

            elif not isinstance(force_constant, _Types._GeneralUnit):
                raise TypeError(
                    "'force_constant' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                )

            # Validate the dimensions.
            if force_constant.dimensions() != (1, 0, -2, 0, 0, -1, 0):
                raise ValueError(
                    "'force_constant' has invalid dimensions! "
                    f"Expected dimensions are 'M Q-1 T-2', found '{force_constant.unit()}'"
                )

        self._force_constant = force_constant

    @classmethod
    def restraints(cls):
        """Return a list of the supported restraint keywords.

        Returns
        -------

        restraints : [str]
            A list of the supported restraint keywords.
        """
        return cls._restraints.copy()
