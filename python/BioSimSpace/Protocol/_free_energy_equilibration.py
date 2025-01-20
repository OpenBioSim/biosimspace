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

"""Functionality for an equilibration alchecmical free-energy protocol."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["FreeEnergyEquilibration"]

from .. import Types as _Types
from .. import Units as _Units

from ._free_energy_mixin import _FreeEnergyMixin
from ._equilibration import Equilibration as _Equilibration


class FreeEnergyEquilibration(_Equilibration, _FreeEnergyMixin):
    """A class for storing free energy equilibration protocols."""

    def __init__(
        self,
        lam=0.0,
        lam_vals=None,
        min_lam=0.0,
        max_lam=1.0,
        num_lam=11,
        timestep=_Types.Time(2, "femtosecond"),
        runtime=_Types.Time(1, "nanoseconds"),
        temperature_start=_Types.Temperature(300, "kelvin"),
        temperature_end=_Types.Temperature(300, "kelvin"),
        temperature=None,
        pressure=None,
        thermostat_time_constant=_Types.Time(1, "picosecond"),
        report_interval=200,
        restart_interval=1000,
        restart=False,
        restraint=None,
        force_constant=10 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
        perturbation_type="full",
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

        timestep : :class:`Time <BioSimSpace.Types.Time>`
            The integration timestep.

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The running time.

        temperature_start : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The starting temperature.

        temperature_end : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The final temperature.

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The equilibration temperature. This takes precedence of over
            the other temperatures, i.e. to run at fixed temperature.

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure. If this argument is omitted then the simulation
            is run using the NVT ensemble.

        thermostat_time_constant : :class:`Time <BioSimSpace.Types.Time>`
            Time constant for thermostat coupling.

        report_interval : int
            The frequency at which statistics are recorded. (In integration steps.)

        restart_interval : int
            The frequency at which restart configurations and trajectory
            frames are saved. (In integration steps.)

        restart : bool
            Whether this is a continuation of a previous simulation.

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
        """

        # Call the base class constructors.

        _Equilibration.__init__(
            self,
            timestep=timestep,
            runtime=runtime,
            temperature_start=temperature_start,
            temperature_end=temperature_end,
            temperature=temperature,
            pressure=pressure,
            thermostat_time_constant=thermostat_time_constant,
            report_interval=report_interval,
            restart_interval=restart_interval,
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
            [_Equilibration._get_parm(self), _FreeEnergyMixin._get_parm(self)]
        )

    def __str__(self):
        """Return a human readable string representation of the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return f"<BioSimSpace.Protocol.FreeEnergyEquilibration: {self._get_parm()}>"

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        if self._is_customised:
            return "<BioSimSpace.Protocol.Custom>"
        else:
            return f"BioSimSpace.Protocol.FreeEnergyEquilibration({self._get_parm()})"

    def __eq__(self, other):
        """Equality operator."""

        if not isinstance(other, FreeEnergyEquilibration):
            return False

        if self._is_customised or other._is_customised:
            return False

        return _Equilibration.__eq__(self, other) and _FreeEnergyMixin.__eq__(
            self, other
        )

    def _to_regular_protocol(self):
        """
        Convert to a regular equilibration protocol. This can be used to run
        'normal' equilibration at a specific lambda end state using SOMD.
        """
        return _Equilibration(
            timestep=self.getTimeStep(),
            runtime=self.getRunTime(),
            temperature_start=self.getStartTemperature(),
            temperature_end=self.getEndTemperature(),
            pressure=self.getPressure(),
            thermostat_time_constant=self.getThermostatTimeConstant(),
            report_interval=self.getReportInterval(),
            restart_interval=self.getRestartInterval(),
            restraint=self.getRestraint(),
            force_constant=self.getForceConstant(),
            restart=self.isRestart(),
        )
