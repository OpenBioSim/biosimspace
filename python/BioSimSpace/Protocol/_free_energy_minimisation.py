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

"""Functionality for a minimisation alchecmical free-energy protocol."""

__all__ = ["FreeEnergyMinimisation"]

from ._free_energy_mixin import _FreeEnergyMixin
from ._minimisation import Minimisation as _Minimisation
from .. import Units as _Units


class FreeEnergyMinimisation(_Minimisation, _FreeEnergyMixin):
    """A class for storing free energy minimisation protocols."""

    def __init__(
        self,
        lam=0.0,
        lam_vals=None,
        min_lam=0.0,
        max_lam=1.0,
        num_lam=11,
        steps=10000,
        perturbation_type="full",
        restraint=None,
        force_constant=10 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
    ):
        """Constructor.

        Parameters
        ----------

        lam : float
            The perturbation parameter: [0.0, 1.0]

        lam_vals : [float]
            The list of lambda parameters.

        min_lam : float
            The minimum lambda value.

        max_lam : float
            The maximum lambda value.

        num_lam : int
            The number of lambda values.

        steps : int
            The maximum number of steps to perform.

        perturbation_type : str
            The type of perturbation to perform. Options are:
             "full" : A full perturbation of all terms (default option).
             "discharge_soft" : Perturb all discharging soft atom charge terms (i.e. value->0.0).
             "vanish_soft" : Perturb all vanishing soft atom LJ terms (i.e. value->0.0).
             "flip" : Perturb all hard atom terms as well as bonds/angles.
             "grow_soft" : Perturb all growing soft atom LJ terms (i.e. 0.0->value).
             "charge_soft" : Perturb all charging soft atom LJ terms (i.e. 0.0->value).

             Currently perturubation_type != "full" is only supported by
             BioSimSpace.Process.Somd.

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

        # Call the base class constructors.

        _Minimisation.__init__(
            self,
            steps=steps,
            restraint=restraint,
            force_constant=force_constant,
        )

        _FreeEnergyMixin.__init__(
            self,
            lam=lam,
            lam_vals=lam_vals,
            min_lam=min_lam,
            max_lam=max_lam,
            num_lam=num_lam,
            perturbation_type=perturbation_type,
        )

    def _get_parm(self):
        """Return a string representation of the parameters."""

        return ", ".join(
            [_Minimisation._get_parm(self), _FreeEnergyMixin._get_parm(self)]
        )

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return f"<BioSimSpace.Protocol.FreeEnergyMinimisation: {self._get_parm()}>"

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._is_customised:
            return "BioSimSpace.Protocol.Custom"
        else:
            return f"BioSimSpace.Protocol.FreeEnergyMinimisation({self._get_parm()})"

    def __eq__(self, other):
        """Equality operator."""

        if not isinstance(other, FreeEnergyMinimisation):
            return False

        if self._is_customised or other._is_customised:
            return False

        return _Minimisation.__eq__(self, other) and _FreeEnergyMixin.__eq__(
            self, other
        )

    def _to_regular_protocol(self):
        """
        Convert to a regular minimisation protocol. This can be used to run
        'normal' minimistion at a specific lambda end state using SOMD.
        """
        return _Minimisation(
            steps=self.getSteps(),
            restraint=self.getRestraint(),
            force_constant=self.getForceConstant(),
        )
