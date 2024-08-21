from BioSimSpace.Units.Area import angstrom2
from BioSimSpace.Units.Energy import kcal_per_mol
from BioSimSpace.Units.Length import angstrom
from .._SireWrappers import System as _System
from .. import Types as _Types
from ._protocol import Protocol as _Protocol
from ._position_restraint_mixin import _PositionRestraintMixin
from .. import Units as _Units
import math as _math
import numpy as _np
import warnings as _warnings
import json as _json
from ..Types import Vector as _Vector

__all__ = ["AToMMinimisation", "AToMEquilibration", "AToMAnnealing", "AToMProduction"]


# When placed in to BSS this needs to be AToM_protocol(protocol):
class _AToM(_Protocol, _PositionRestraintMixin):
    def __init__(
        self,
        system=None,
        data=None,
        core_alignment=True,
        com_distance_restraint=True,
        com_k=25.0 * kcal_per_mol / angstrom2,
        com_restraint_width=5.0 * angstrom,
        restraint=None,
        force_constant=10 * kcal_per_mol / angstrom2,
        positional_restraint_width=0.5 * angstrom,
        align_kf_sep=2.5 * kcal_per_mol / angstrom2,
        align_k_theta=10.0 * kcal_per_mol,
        align_k_psi=10.0 * kcal_per_mol,
        SC_umax=1000.0 * kcal_per_mol,
        SC_u0=500.0 * kcal_per_mol,
        SC_a=0.0625,
    ):
        # Call the base class constructor.
        super().__init__()

        # first check that EITHER system or data is passed
        if system is None and data is None:
            raise ValueError(
                "Either 'system' or 'data' must be passed to the AToM protocol."
            )

        if system is not None and not isinstance(system, _System):
            raise TypeError("'system' must be of type 'BioSimSpace.System'")

        if data is not None and not isinstance(data, dict):
            raise TypeError("'data' must be of type 'dict'")

        if isinstance(system, _System) and data is None:
            try:
                sdata = _json.loads(system._sire_object.property("atom_data").value())
            except Exception as e:
                raise ValueError(
                    f"Unable to extract AToM data from the system object. The following error was raised: {e}."
                )
            # convert the "displacement" key back to a vector
            d = sdata["displacement"]
            displacement = _Vector(*d)
            sdata["displacement"] = displacement
            self._system_data = sdata

        elif system is not None and data is not None:
            _warnings.warn(
                "Both 'system' and 'data' were passed. Using 'data' and ignoring data from 'system'."
            )

        # Store the AToM system.
        if isinstance(data, dict):
            self._system_data = data
        elif data is not None:
            raise TypeError("'data' must be of type 'dict'")

        # Whether or not to use alignment restraints.
        self.setCoreAlignment(core_alignment)

        # Whether or not to use the CMCM restraint.
        self.setCOMDistanceRestraint(com_distance_restraint)

        # Store com_k value.
        self.setCOMk(com_k)

        # Store com_restraint_width value.
        self.setCOMWidth(com_restraint_width)

        # Store the width of the coordinate restraint.
        self.setPosRestWidth(positional_restraint_width)

        # Store the align_kf_sep value.
        self.setAlignKfSep(align_kf_sep)

        # Store the align_k_theta value.
        self.setAlignKTheta(align_k_theta)

        # Store the align_k_psi value.
        self.setAlignKPsi(align_k_psi)

        # Store the SC_umax value.
        self.setSCUmax(SC_umax)

        # Store the SC_u0 value.
        self.setSCU0(SC_u0)

        # Store the SC_a value.
        self.setSCa(SC_a)

        # Set the postition restraint.
        _PositionRestraintMixin.__init__(self, restraint, force_constant)

    def __str__(self):
        d = self.getData()
        """Return a string representation of the protocol."""
        string = "<BioSimSpace.Protocol.AToM>: "
        string += "timestep=%s " % self.getTimeStep()
        string += ", runtime=%s " % self.getRunTime()
        string += ", temperature=%s " % self.getTemperature()
        if self._pressure is not None:
            string += ", pressure=%s, " % self.getPressure()
        string += ", lambda1=%s " % self.getLambda1()
        string += ", lambda2=%s " % self.getLambda2()
        string += ", ligand1 core atoms=%s" % d["ligand1_rigid_core"]
        string += ", ligand2 core atoms=%s" % d["ligand2_rigid_core"]
        string += ", report_interval=%s " % self.getReportInterval()
        string += ", restart_interval=%s " % self.getRestartInterval()
        string += ">"

        return string

    def __repr__(self):
        """Return a string showing how to instantiate the object."""
        return self.__str__()

    def getData(self):
        """
        Return the AToM data dictionary.

        Returns
        -------

        data : dict
            The AToM data dictionary.
        """
        return self._system_data

    def getCoreAlignment(self):
        """
        Return core alignment boolean.

        Returns
        -------

        core_alignment : bool
            Whether to use core alignment.
        """
        return self._core_alignment

    def setCoreAlignment(self, core_alignment):
        """
        Set the core alignment flag.

        Parameters
        ----------

        core_alignment : bool
            Whether to use core alignment.
        """
        if isinstance(core_alignment, bool):
            self._core_alignment = core_alignment
        else:
            _warnings.warn("Non-boolean core alignment flag. Defaulting to True!")
            self._core_alignment = True

    def getCOMDistanceRestraint(self):
        """
        Return CMCM restraint boolean.

        Returns
        -------

        com_distance_restraint : bool
            Whether to use the CMCM restraint.
        """
        return self._com_distance_restraint

    def setCOMDistanceRestraint(self, com_distance_restraint):
        """
        Set the CMCM restraint flag.

        Parameters
        ----------

        com_distance_restraint : bool
            Whether to use the CMCM restraint.
        """
        if isinstance(com_distance_restraint, bool):
            self._com_distance_restraint = com_distance_restraint
        else:
            _warnings.warn(
                "Non-boolean com distance restraint flag. Defaulting to True!"
            )
            self._com_distance_restraint = True

    def getPosRestWidth(self):
        """
        Return the width of the position restraint.

        Returns
        -------

        positional_restraint_width : :class:`Length <BioSimSpace.Types.Length>`
            The width of the position restraint.
        """
        return self._positional_restraint_width

    def setPosRestWidth(self, positional_restraint_width):
        """
        Set the width of the position restraint.

        Parameters
        ----------

        positional_restraint_width : int, float, str, :class:`Length <BioSimSpace.Types.Length>`
            The width of the position restraint.
        """
        # Convert int to float.
        if type(positional_restraint_width) is int:
            positional_restraint_width = float(positional_restraint_width)

        if isinstance(positional_restraint_width, float):
            # Use default units.
            positional_restraint_width *= _Units.Length.angstrom

        else:
            if isinstance(positional_restraint_width, str):
                try:
                    positional_restraint_width = _Types.Length(
                        positional_restraint_width
                    )
                except Exception:
                    raise ValueError(
                        "Unable to parse 'positional_restraint_width' string."
                    ) from None

            elif not isinstance(positional_restraint_width, _Types.Length):
                raise TypeError(
                    "'positional_restraint_width' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                )

            # Validate the dimensions.
            if positional_restraint_width.dimensions() != (0, 1, 0, 0, 0, 0, 0):
                raise ValueError(
                    "'positional_restraint_width' has invalid dimensions! "
                    f"Expected dimensions of Length, found '{positional_restraint_width.unit()}'"
                )
        self._positional_restraint_width = positional_restraint_width

    def getAlignKfSep(self):
        """
        Return the align_kf_sep value.

        Returns
        -------

        align_kf_sep : :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`
            The align_kf_sep value in kcal/mol angstrom**2.
        """
        return self._align_kf_sep

    def setAlignKfSep(self, align_kf_sep):
        """
        Set the align_kf_sep value.

        Parameters
        ----------

        align_kf_sep : int, float, str, :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`, float
            Length value for the alignment restraint kcal/mol angstrom**2.
        """
        # Convert int to float.
        if type(align_kf_sep) is int:
            align_kf_sep = float(align_kf_sep)

        if isinstance(align_kf_sep, float):
            # Use default units.
            align_kf_sep *= _Units.Energy.kcal_per_mol / _Units.Area.angstrom2

        else:
            if isinstance(align_kf_sep, str):
                try:
                    align_kf_sep = _Types._GeneralUnit(align_kf_sep)
                except Exception:
                    raise ValueError("Unable to parse 'align_kf_sep' string.") from None

            elif not isinstance(align_kf_sep, _Types._GeneralUnit):
                raise TypeError(
                    "'align_kf_sep' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                )

            # Validate the dimensions.
            if align_kf_sep.dimensions() != (1, 0, -2, 0, 0, -1, 0):
                raise ValueError(
                    "'align_kf_sep' has invalid dimensions! "
                    f"Expected dimensions of energy density/area (e.g. kcal/molA^2), found '{align_kf_sep.unit()}'"
                )
        self._align_kf_sep = align_kf_sep

    def getAlignKTheta(self):
        """
        Return the align_k_theta value.

        Returns
        -------

        align_k_theta : :class:`Energy <BioSimSpace.Types.Energy>`
            The align_k_theta value in kcal/mol.
        """
        return self._align_k_theta

    def setAlignKTheta(self, align_k_theta):
        """
        Set the align_k_theta value.

        Parameters
        ----------

        align_k_theta : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            Force constant for the alignment angular constraint in kcal/mol.

        """
        # Convert int to float.
        if type(align_k_theta) is int:
            align_k_theta = float(align_k_theta)

        if isinstance(align_k_theta, float):
            # Use default units.
            align_k_theta *= _Units.Energy.kcal_per_mol

        else:
            if isinstance(align_k_theta, str):
                try:
                    align_k_theta = _Types._GeneralUnit(align_k_theta)
                except Exception:
                    raise ValueError(
                        "Unable to parse 'align_k_theta' string."
                    ) from None

            elif not isinstance(align_k_theta, _Types.Energy):
                raise TypeError(
                    "'align_k_theta' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                )

            # Validate the dimensions.
            if align_k_theta.dimensions() != (1, 2, -2, 0, 0, -1, 0):
                raise ValueError(
                    "'align_k_theta' has invalid dimensions! "
                    f"Expected dimensions of energy density (e.g. kcal/mol), found '{align_k_theta.unit()}'"
                )
        self._align_k_theta = align_k_theta

    def getAlignKPsi(self):
        """
        Return the align_k_psi value.

        Returns
        -------

        align_k_psi: :class:`Energy <BioSimSpace.Types.Energy>`
            The align_k_psi value in kcal/mol.
        """
        return self._align_k_psi

    def setAlignKPsi(self, align_k_psi):
        """
        Set the align_k_psi value.

        Parameters
        ----------

        align_k_psi : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            Force constant for the alignment dihedral constraint in kcal/mol.
        """
        # Convert int to float.
        if type(align_k_psi) is int:
            align_k_psi = float(align_k_psi)

        if isinstance(align_k_psi, float):
            # Use default units.
            align_k_psi *= _Units.Energy.kcal_per_mol

        else:
            if isinstance(align_k_psi, str):
                try:
                    align_k_psi = _Types._GeneralUnit(align_k_psi)
                except Exception:
                    raise ValueError("Unable to parse 'align_k_psi' string.") from None

            elif not isinstance(align_k_psi, _Types.Energy):
                raise TypeError(
                    "'align_k_psi' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                )

            # Validate the dimensions.
            if align_k_psi.dimensions() != (1, 2, -2, 0, 0, -1, 0):
                raise ValueError(
                    "'align_k_psi' has invalid dimensions! "
                    f"Expected dimensions of energy density (e.g. kcal/mol), found '{align_k_psi.unit()}'"
                )
        self._align_k_psi = align_k_psi

    def getSCUmax(self):
        """
        Return the SC_umax value.

        Returns
        -------

        SC_umax : :class:`Energy <BioSimSpace.Types.Energy>`
            The SC_umax value in kcal/mol.
        """
        return self._SC_umax

    def setSCUmax(self, SC_umax):
        """
        Set the SC_umax value.

        Parameters
        ----------

        SC_umax : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The softcore Umax value in kcal/mol.
        """
        # Convert int to float.
        if type(SC_umax) is int:
            SC_umax = float(SC_umax)

        if isinstance(SC_umax, float):
            # Use default units.
            SC_umax *= _Units.Energy.kcal_per_mol

        else:
            if isinstance(SC_umax, str):
                try:
                    SC_umax = _Types._GeneralUnit(SC_umax)
                except Exception:
                    raise ValueError("Unable to parse 'SC_umax' string.") from None

            elif not isinstance(SC_umax, _Types.Energy):
                raise TypeError(
                    "'SC_umax' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                )

            # Validate the dimensions.
            if SC_umax.dimensions() != (1, 2, -2, 0, 0, -1, 0):
                raise ValueError(
                    "'align_k_theta' has invalid dimensions! "
                    f"Expected dimensions of energy density (e.g. kcal/mol), found '{SC_umax.unit()}'"
                )
        self._SC_umax = SC_umax

    def getSCU0(self):
        """
        Return the SC_u0 value.

        Returns
        -------

        SC_u0 : :class:`Energy <BioSimSpace.Types.Energy>`
            The SC_u0 value in kcal/mol.
        """
        return self._SC_u0

    def setSCU0(self, SC_u0):
        """
        Set the SC_u0 value.

        Parameters
        ----------

        SC_u0 : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The softcore u0 value in kcal/mol.
        """
        # Convert int to float.
        if type(SC_u0) is int:
            SC_u0 = float(SC_u0)

        if isinstance(SC_u0, float):
            # Use default units.
            SC_u0 *= _Units.Energy.kcal_per_mol

        else:
            if isinstance(SC_u0, str):
                try:
                    SC_u0 = _Types._GeneralUnit(SC_u0)
                except Exception:
                    raise ValueError("Unable to parse 'SC_u0' string.") from None

            elif not isinstance(SC_u0, _Types.Energy):
                raise TypeError(
                    "'SC_u0' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                )

            # Validate the dimensions.
            if SC_u0.dimensions() != (1, 2, -2, 0, 0, -1, 0):
                raise ValueError(
                    "'align_k_theta' has invalid dimensions! "
                    f"Expected dimensions of energy density (e.g. kcal/mol), found '{SC_u0.unit()}'"
                )
        self._SC_u0 = SC_u0

    def getSCa(self):
        """
        Return the SC_a value.

        Returns
        -------

        SC_a : float
            The SC_a value.
        """
        return self._SC_a

    def setSCa(self, SC_a):
        """
        Set the SC_a value.

        Parameters
        ----------

        SC_a : float
            The softcore a value.
        """
        if isinstance(SC_a, (int, float)):
            self._SC_a = float(SC_a)
        else:
            raise TypeError("'SC_a' must be of type 'float'")

    def getCOMk(self):
        """
        Return the com_k value.

        Returns
        -------

        com_k : :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`
            The com_k value in kcal/mol A**2.
        """
        return self._com_k

    def setCOMk(self, com_k):
        """
        Set the com_k value.

        Parameters
        ----------

        com_k : int, float, str, :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>
            The force constant for the CM-CM force in kcal/mol A**2.
        """
        # Convert int to float.
        if type(com_k) is int:
            com_k = float(com_k)

        if isinstance(com_k, float):
            # Use default units.
            com_k *= _Units.Energy.kcal_per_mol / _Units.Area.angstrom2

        else:
            if isinstance(com_k, str):
                try:
                    com_k = _Types._GeneralUnit(com_k)
                except Exception:
                    raise ValueError("Unable to parse 'com_k' string.") from None

            elif not isinstance(com_k, _Types._GeneralUnit):
                raise TypeError(
                    "'com_k' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                )

            # Validate the dimensions.
            if com_k.dimensions() != (1, 0, -2, 0, 0, -1, 0):
                raise ValueError(
                    "'align_k_theta' has invalid dimensions! "
                    f"Expected dimensions of energy density/area (e.g. kcal/molA^2), found '{com_k.unit()}'"
                )
        self._com_k = com_k

    def getCOMWidth(self):
        """
        Return the com_restraint_width value.

        Returns
        -------

        com_restraint_width : :class:`Length <BioSimSpace.Types.Length>`
            The com_restraint_width value in angstroms.
        """
        return self._com_restraint_width

    def setCOMWidth(self, com_restraint_width):
        """
        Set the com_restraint_width value.

        Parameters
        ----------

        com_restraint_width : int, float, str, :class:`Length <BioSimSpace.Types.Length>
            The com_restraint_width value in angstroms.
        """
        # Convert int to float.
        if type(com_restraint_width) is int:
            com_restraint_width = float(com_restraint_width)

        if isinstance(com_restraint_width, float):
            # Use default units.
            com_restraint_width *= _Units.Length.angstrom

        else:
            if isinstance(com_restraint_width, str):
                try:
                    com_restraint_width = _Types.Length(com_restraint_width)
                except Exception:
                    raise ValueError(
                        "Unable to parse 'com_restraint_width' string."
                    ) from None

            elif not isinstance(com_restraint_width, _Types.Length):
                raise TypeError(
                    "'com_restraint_width' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                )

            # Validate the dimensions.
            if com_restraint_width.dimensions() != (0, 1, 0, 0, 0, 0, 0):
                raise ValueError(
                    "'align_k_theta' has invalid dimensions! "
                    f"Expected dimensions of Length, found '{com_restraint_width.unit()}'"
                )
        self._com_restraint_width = com_restraint_width


class AToMMinimisation(_AToM):
    """
    Minimisation protocol for AToM simulations.
    """

    def __init__(
        self,
        system=None,
        data=None,
        steps=10000,
        core_alignment=True,
        com_distance_restraint=True,
        restraint=None,
        force_constant=10 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
        positional_restraint_width=0.5 * angstrom,
        align_kf_sep=2.5 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
        align_k_theta=10 * _Units.Energy.kcal_per_mol,
        align_k_psi=10 * _Units.Energy.kcal_per_mol,
        SC_umax=1000 * _Units.Energy.kcal_per_mol,
        SC_u0=500 * _Units.Energy.kcal_per_mol,
        SC_a=0.0625,
        com_k=25 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
        com_restraint_width=5 * _Units.Length.angstrom,
    ):
        """
        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
                A prepared AToM system.
        data : dict
            The AToM data dictionary.

        core_alignment : bool
            Whether to use rigid core restraints to align the two ligands.

        com_distance_restraint : bool
            Whether to use a center of mass distance restraint.
            This restraint applies to the protein/host and both ligands, and
            is used to maintain the relative positions of all of them.

        com_k : int, float, str, :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`
            The force constant for the center of mass distance restraint (kcal/mol/A^2).

        com_restraint_width : int, float, str, :class:`Length <BioSimSpace.Types.Length>
            The width (tolerance) of the center of mass distance restraint (A).

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

        positional_restraint_width : :class:`Length <BioSimSpace.Types.Length>`, float
                The width of the flat-bottom potential used for coordinate restraint in Angstroms.

        pos_restrained_atoms : [int]
            The atoms to be restrained.

        align_kf_sep : int, float, str, :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`
            The force constant for the distance portion of the alignment restraint (kcal/(mol A^2)).

        align_k_theta : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The force constant for the angular portion of the alignment restaint (kcal/mol).

        align_k_psi : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The force constant for the dihedral portion of the alignment restraint (kcal/mol).

        SC_umax : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The Umax value for the ATM softcore potential (kcal/mol).

        SC_u0 : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The uh value for the ATM softcore potential (kcal/mol).

        SC_a : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The a value for the ATM softcore potential."""

        super().__init__(
            system,
            data,
            core_alignment,
            com_distance_restraint,
            com_k,
            com_restraint_width,
            restraint,
            force_constant,
            positional_restraint_width,
            align_kf_sep,
            align_k_theta,
            align_k_psi,
            SC_umax,
            SC_u0,
            SC_a,
        )
        # Store the number of minimisation steps.
        self.setSteps(steps)

    def getSteps(self):
        """
        Return the number of minimisation steps.

        Returns
        -------

        steps : int
            The number of minimisation steps.
        """
        return self._steps

    def setSteps(self, steps):
        """
        Set the number of minimisation steps.

        Parameters
        ----------

        steps : int
            The number of minimisation steps.
        """
        if isinstance(steps, int):
            self._steps = steps
        else:
            raise TypeError("'steps' must be of type 'int'")


class AToMEquilibration(_AToM):
    """Equilibration protocol for AToM simulations."""

    def __init__(
        self,
        system=None,
        data=None,
        timestep=2 * _Units.Time.femtosecond,
        runtime=0.2 * _Units.Time.nanosecond,
        temperature_start=300 * _Units.Temperature.kelvin,
        temperature_end=300 * _Units.Temperature.kelvin,
        temperature=None,
        pressure=1 * _Units.Pressure.atm,
        thermostat_time_constant=1 * _Units.Time.picosecond,
        report_interval=100,
        restart_interval=100,
        core_alignment=True,
        com_distance_restraint=True,
        com_k=25 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
        com_restraint_width=5 * _Units.Length.angstrom,
        restraint=None,
        force_constant=10 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
        positional_restraint_width=0.5 * angstrom,
        align_kf_sep=2.5 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
        align_k_theta=10 * _Units.Energy.kcal_per_mol,
        align_k_psi=10 * _Units.Energy.kcal_per_mol,
        SC_umax=1000 * _Units.Energy.kcal_per_mol,
        SC_u0=500 * _Units.Energy.kcal_per_mol,
        SC_a=0.0625,
        use_atm_force=False,
        direction=1,
        lambda1=0.0,
        lambda2=0.0,
        alpha=0.0 * _Units.Energy.kcal_per_mol,
        uh=0.0 * _Units.Energy.kcal_per_mol,
        W0=0.0 * _Units.Energy.kcal_per_mol,
    ):
        """
        Create a new equilibration protocol.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>``
            A prepared AToM system.

        data : dict
            The AToM data dictionary.

        timestep : str, :class:`Time <BioSimSpace.Types.Time>`
            The integration timestep.

        runtime : str, :class:`Time <BioSimSpace.Types.Time>`
            The running time.

        temperature_start : str, :class:`Temperature <BioSimSpace.Types.Temperature>`
            The starting temperature.

        temperature_end : str, :class:`Temperature <BioSimSpace.Types.Temperature>`
            The ending temperature.

        temperature : str, :class:`Temperature <BioSimSpace.Types.Temperature>`
             The equilibration temperature. This takes precedence of over the other temperatures, i.e. to run at fixed temperature.

        pressure : str, :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure. Pass pressure=None to use the NVT ensemble.

        thermostat_time_constant : str, :class:`Time <BioSimSpace.Types.Time>`
            Time constant for thermostat coupling.

        report_interval : int
            The frequency at which statistics are recorded. (In integration steps.)

        restart_interval : int
            The frequency at which restart configurations and trajectory

        core_alignment : bool
            Whether to use rigid core restraints to align the two ligands.

        com_distance_restraint : bool
            Whether to use a center of mass distance restraint.
            This restraint applies to the protein/host and both ligands, and
            is used to maintain the relative positions of all of them.

        com_k : int, float, str, :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`
            The force constant for the center of mass distance restraint (kcal/mol/A^2).

        com_restraint_width : int, float, str, :class:`Length <BioSimSpace.Types.Length>
            The width (tolerance) of the center of mass distance restraint (A).

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

        force_constant : float, :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`
            The force constant for the restraint potential (kcal/(mol A^2).

        positional_restraint_width :  float, :class:`Length <BioSimSpace.Types.Length>`
            The width of the flat-bottom potential used for coordinate restraint in Angstroms.

        pos_restrained_atoms : [int]
            The atoms to be restrained.

        align_kf_sep : int, float, str, :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`
            The force constant for the distance portion of the alignment restraint (kcal/(mol A^2).

        align_k_theta : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The force constant for the angular portion of the alignment restaint (kcal/mol).

        align_k_psi : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The force constant for the dihedral portion of the alignment restraint (kcal/mol).

        SC_umax : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The Umax value for the ATM softcore potential (kcal/mol).

        SC_u0 : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The uh value for the ATM softcore potential (kcal/mol).

        SC_a : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The a value for the ATM softcore potential.


        use_atm_force : bool
            Whether to apply the ATM force within the equilibration protocol.

        direction : str
            The direction of the equilibration. Ignored if use_atm_force is False.

        lambda1 : float
            The lambda1 value for the ATM force. Ignored if use_atm_force is False.

        lambda2 : float
            The lambda2 value for the ATM force. Ignored if use_atm_force is False.

        alpha : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The alpha value for the ATM force. Ignored if use_atm_force is False.
            Value in kcal/mol.

        uh : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The uh value for the ATM force. Ignored if use_atm_force is False.
            Value in kcal/mol.

        W0 : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The W0 value for the ATM force. Ignored if use_atm_force is False.
            Value in kcal/mol.
        """
        super().__init__(
            system=system,
            data=data,
            core_alignment=core_alignment,
            com_distance_restraint=com_distance_restraint,
            com_k=com_k,
            com_restraint_width=com_restraint_width,
            restraint=restraint,
            force_constant=force_constant,
            positional_restraint_width=positional_restraint_width,
            align_kf_sep=align_kf_sep,
            align_k_theta=align_k_theta,
            align_k_psi=align_k_psi,
            SC_umax=SC_umax,
            SC_u0=SC_u0,
            SC_a=SC_a,
        )
        # Store
        self.setTimestep(timestep)

        self.setRuntime(runtime)
        # Constant temperature equilibration.
        if temperature is not None:
            self.setStartTemperature(temperature)
            self.setEndTemperature(temperature)
            self._is_const_temp = True

        # Heating / cooling simulation.
        else:
            self._is_const_temp = False

            # Set the start temperature.
            self.setStartTemperature(temperature_start)

            # Set the final temperature.
            self.setEndTemperature(temperature_end)

            # Constant temperature simulation.
            if self._temperature_start == self._temperature_end:
                self._is_const_temp = True

        # Set the system pressure.
        if pressure is not None:
            self.setPressure(pressure)
        else:
            self._pressure = None

        self.setThermostatTimeConstant(thermostat_time_constant)

        self.setReportInterval(report_interval)

        self.setRestartInterval(restart_interval)

        self.setUseATMForce(use_atm_force)

        self.setDirection(direction)

        self.setLambda1(lambda1)

        self.setLambda2(lambda2)

        self.setAlpha(alpha)

        self.setUh(uh)

        self.setW0(W0)

    def getTimeStep(self):
        """
        Return the time step.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The integration time step.
        """
        return self._timestep

    def setTimestep(self, timestep):
        """
        Set the time step.

        Parameters
        ----------

        time : str, :class:`Time <BioSimSpace.Types.Time>`
            The integration time step.
        """
        if isinstance(timestep, str):
            try:
                self._timestep = _Types.Time(timestep)
            except:
                raise ValueError("Unable to parse 'timestep' string.") from None
        elif isinstance(timestep, _Types.Time):
            self._timestep = timestep
        else:
            raise TypeError(
                "'timestep' must be of type 'str' or 'BioSimSpace.Types.Time'"
            )

    def getRunTime(self):
        """
        Return the running time.

        Returns
        -------

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The simulation run time.
        """
        return self._runtime

    def setRuntime(self, runtime):
        """
        Set the running time.

        Parameters
        ----------

        runtime : str, :class:`Time <BioSimSpace.Types.Time>`
            The simulation run time.
        """
        if isinstance(runtime, str):
            try:
                self._runtime = _Types.Time(runtime)
            except:
                raise ValueError("Unable to parse 'runtime' string.") from None
        elif isinstance(runtime, _Types.Time):
            self._runtime = runtime
        else:
            raise TypeError(
                "'runtime' must be of type 'str' or 'BioSimSpace.Types.Time'"
            )

    def getStartTemperature(self):
        """
        Return the starting temperature.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The starting temperature.
        """
        return self._temperature_start

    def setStartTemperature(self, temperature):
        """
        Set the starting temperature.

        Parameters
        ----------

        temperature : str, :class:`Temperature <BioSimSpace.Types.Temperature>`
            The starting temperature.
        """

        if isinstance(temperature, str):
            try:
                temperature = _Types.Temperature(temperature)
            except:
                raise ValueError("Unable to parse 'temperature' string.") from None
        elif not isinstance(temperature, _Types.Temperature):
            raise TypeError(
                "'temperature' must be of type 'str' or 'BioSimSpace.Types.Temperature'"
            )

        if _math.isclose(temperature.kelvin().value(), 0, rel_tol=1e-6):
            temperature._value = 0.01
        self._temperature_start = temperature

    def getEndTemperature(self):
        """
        Return the final temperature.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The final temperature.
        """
        return self._temperature_end

    def setEndTemperature(self, temperature):
        """
        Set the final temperature.

        Parameters
        ----------

        temperature : str, :class:`Temperature <BioSimSpace.Types.Temperature>`
            The final temperature.
        """
        if isinstance(temperature, str):
            try:
                temperature = _Types.Temperature(temperature)
            except:
                raise ValueError("Unable to parse 'temperature' string.") from None
        elif not isinstance(temperature, _Types.Temperature):
            raise TypeError(
                "'temperature' must be of type 'str' or 'BioSimSpace.Types.Temperature'"
            )

        if _math.isclose(temperature.kelvin().value(), 0, rel_tol=1e-6):
            temperature._value = 0.01
        self._temperature_end = temperature

    def getPressure(self):
        """
        Return the pressure.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure.
        """
        return self._pressure

    def setPressure(self, pressure):
        """
        Set the pressure.

        Parameters
        ----------

        pressure : str, :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure.
        """
        if isinstance(pressure, str):
            try:
                self._pressure = _Types.Pressure(pressure)
            except:
                raise ValueError("Unable to parse 'pressure' string.") from None
        elif isinstance(pressure, _Types.Pressure):
            self._pressure = pressure
        else:
            raise TypeError(
                "'pressure' must be of type 'str' or 'BioSimSpace.Types.Pressure'"
            )

    def getThermostatTimeConstant(self):
        """
        Return the time constant for the thermostat.

        Returns
        -------

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The time constant for the thermostat.
        """
        return self._thermostat_time_constant

    def setThermostatTimeConstant(self, thermostat_time_constant):
        """
        Set the time constant for the thermostat.

        Parameters
        ----------

        thermostat_time_constant : str, :class:`Time <BioSimSpace.Types.Time>`
            The time constant for the thermostat.
        """
        if isinstance(thermostat_time_constant, str):
            try:
                self._thermostat_time_constant = _Types.Time(thermostat_time_constant)
            except:
                raise ValueError(
                    "Unable to parse 'thermostat_time_constant' string."
                ) from None
        elif isinstance(thermostat_time_constant, _Types.Time):
            self._thermostat_time_constant = thermostat_time_constant
        else:
            raise TypeError(
                "'thermostat_time_constant' must be of type 'BioSimSpace.Types.Time'"
            )

    def getReportInterval(self):
        """
        Return the interval between reporting statistics. (In integration steps.).
        Returns
        -------
        report_interval : int
            The number of integration steps between reporting statistics.
        """
        return self._report_interval

    def setReportInterval(self, report_interval):
        """
        Set the interval at which statistics are reported. (In integration steps.).

        Parameters
        ----------

        report_interval : int
            The number of integration steps between reporting statistics.
        """
        if not type(report_interval) is int:
            raise TypeError("'report_interval' must be of type 'int'")

        if report_interval <= 0:
            _warnings.warn("'report_interval' must be positive. Using default (100).")
            report_interval = 100

        self._report_interval = report_interval

    def getRestartInterval(self):
        """
        Return the interval between saving restart confiugrations, and/or
        trajectory frames. (In integration steps.).

        Returns
        -------

        restart_interval : int
            The number of integration steps between saving restart
            configurations and/or trajectory frames.
        """
        return self._restart_interval

    def setRestartInterval(self, restart_interval):
        """
        Set the interval between saving restart confiugrations, and/or
        trajectory frames. (In integration steps.).

        Parameters
        ----------

        restart_interval : int
            The number of integration steps between saving restart
            configurations and/or trajectory frames.
        """
        if not type(restart_interval) is int:
            raise TypeError("'restart_interval' must be of type 'int'")

        if restart_interval <= 0:
            _warnings.warn("'restart_interval' must be positive. Using default (500).")
            restart_interval = 500

        self._restart_interval = restart_interval

    def getUseATMForce(self):
        """
        Return the use_atm_force flag.

        Returns
        -------

        use_atm_force : bool
            Whether to apply the ATM force within the equilibration protocol.
        """
        return self._use_atm_force

    def setUseATMForce(self, use_atm_force):
        """
        Set the use_atm_force flag.

        Parameters
        ----------

        use_atm_force : bool
            Whether to apply the ATM force within the equilibration protocol.
        """
        if not isinstance(use_atm_force, bool):
            raise TypeError("'use_atm_force' must be of type 'bool'")
        self._use_atm_force = use_atm_force

    def getDirection(self):
        """
        Return the direction of the equilibration.

        Returns
        -------

        direction : str
            The direction of the equilibration. Ignored if use_atm_force is False.
        """
        return self._direction

    def setDirection(self, direction):
        """
        Set the direction of the equilibration.

        Parameters
        ----------

        direction : str
            The direction of the equilibration. Ignored if use_atm_force is False.
        """
        if int(direction) != 1 and int(direction) != -1:
            raise TypeError("'direction' must have a value of 1 or -1")
        self._direction = int(direction)

    def getLambda1(self):
        """
        Return the lambda1 value for the ATM force.

        Returns
        -------

        lambda1 : float
            The lambda1 value for the ATM force. Ignored if use_atm_force is False.
        """
        return self._lambda1

    def setLambda1(self, lambda1):
        """
        Set the lambda1 value for the ATM force.

        Parameters
        ----------

        lambda1 : float
            The lambda1 value for the ATM force. Ignored if use_atm_force is False.
        """
        if not isinstance(lambda1, (float, int)):
            raise TypeError("'lambda1' must be of type 'float'")
        if not 0 <= float(lambda1) <= 0.5:
            raise ValueError("lambda1 must be between 0 and 0.5")
        self._lambda1 = float(lambda1)

    def getLambda2(self):
        """
        Return the lambda2 value for the ATM force.

        Returns
        -------

        lambda2 : float
            The lambda2 value for the ATM force. Ignored if use_atm_force is False.
        """
        return self._lambda2

    def setLambda2(self, lambda2):
        """
        Set the lambda2 value for the ATM force.

        Parameters
        ----------

        lambda2 : float
            The lambda2 value for the ATM force. Ignored if use_atm_force is False.
        """
        if not isinstance(lambda2, (float, int)):
            raise TypeError("'lambda2' must be of type 'float'")
        if not 0 <= float(lambda2) <= 0.5:
            raise ValueError("lambda2 must be between 0 and 0.5")
        self._lambda2 = float(lambda2)

    def getAlpha(self):
        """
        Return the alpha value for the ATM force.

        Returns
        -------

        alpha : :class:`Energy <BioSimSpace.Types.Energy>`
            The alpha value for the ATM force in kcal/mol. Ignored if use_atm_force is False.
        """
        return self._alpha

    def setAlpha(self, alpha):
        """
        Set the alpha value for the ATM force.

        Parameters
        ----------

        alpha : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The alpha value for the ATM force in kcal/mol. Ignored if use_atm_force is False.
        """
        # Convert int to float.
        if type(alpha) is int:
            alpha = float(alpha)

        if isinstance(alpha, float):
            # Use default units.
            alpha *= _Units.Energy.kcal_per_mol

        else:
            if isinstance(alpha, str):
                try:
                    alpha = _Types._GeneralUnit(alpha)
                except Exception:
                    raise ValueError("Unable to parse 'alpha' string.") from None

            elif not isinstance(alpha, _Types.Energy):
                raise TypeError(
                    "'alpha' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                )

            # Validate the dimensions.
            if alpha.dimensions() != (1, 2, -2, 0, 0, -1, 0):
                raise ValueError(
                    "'align_k_theta' has invalid dimensions! "
                    f"Expected dimensions of energy density (e.g. kcal/mol), found '{alpha.unit()}'"
                )
        self._alpha = alpha

    def getUh(self):
        """
        Return the uh value for the ATM force.

        Returns
        -------

        uh : :class:`Energy <BioSimSpace.Types.Energy>`
            The uh value for the ATM force in kcal/mol. Ignored if use_atm_force is False.
        """
        return self._uh

    def setUh(self, uh):
        """
        Set the uh value for the ATM force.

        Parameters
        ----------

        uh : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The uh value for the ATM force in kcal/mol. Ignored if use_atm_force is False.
        """
        # Convert int to float.
        if type(uh) is int:
            uh = float(uh)

        if isinstance(uh, float):
            # Use default units.
            uh *= _Units.Energy.kcal_per_mol

        else:
            if isinstance(uh, str):
                try:
                    uh = _Types._GeneralUnit(uh)
                except Exception:
                    raise ValueError("Unable to parse 'uh' string.") from None

            elif not isinstance(uh, _Types.Energy):
                raise TypeError(
                    "'uh' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                )

            # Validate the dimensions.
            if uh.dimensions() != (1, 2, -2, 0, 0, -1, 0):
                raise ValueError(
                    "'align_k_theta' has invalid dimensions! "
                    f"Expected dimensions of energy density (e.g. kcal/mol), found '{uh.unit()}'"
                )
        self._uh = uh

    def getW0(self):
        """
        Return the W0 value for the ATM force.

        Returns
        -------

        W0 : :class:`Energy <BioSimSpace.Types.Energy>`
            The W0 value for the ATM force in kcal/mol. Ignored if use_atm_force is False.
        """
        return self._W0

    def setW0(self, W0):
        """
        Set the W0 value for the ATM force.

        Parameters
        ----------

        W0 :int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The W0 value for the ATM force in kcal/mol. Ignored if use_atm_force is False.
        """
        # Convert int to float.
        if type(W0) is int:
            W0 = float(W0)

        if isinstance(W0, float):
            # Use default units.
            W0 *= _Units.Energy.kcal_per_mol

        else:
            if isinstance(W0, str):
                try:
                    W0 = _Types._GeneralUnit(W0)
                except Exception:
                    raise ValueError("Unable to parse 'W0' string.") from None

            elif not isinstance(W0, _Types.Energy):
                raise TypeError(
                    "'W0' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                )

            # Validate the dimensions.
            if W0.dimensions() != (1, 2, -2, 0, 0, -1, 0):
                raise ValueError(
                    "'align_k_theta' has invalid dimensions! "
                    f"Expected dimensions of energy density (e.g. kcal/mol), found '{W0.unit()}'"
                )
        self._W0 = W0

    def isConstantTemp(self):
        """
        Return whether the protocol has a constant temperature.

        Returns
        -------

        is_const_temp : bool
            Whether the temperature is fixed.
        """
        return self._temperature_start == self._temperature_end

    @classmethod
    def restraints(cls):
        """
        Return a list of the supported restraint keywords.

        Returns
        -------

        restraints : [str]
            A list of the supported restraint keywords.
        """
        return cls._restraints.copy()


class AToMAnnealing(_AToM):
    """Annealing protocol for AToM simulations."""

    def __init__(
        self,
        system=None,
        data=None,
        timestep=2 * _Units.Time.femtosecond,
        runtime=0.2 * _Units.Time.nanosecond,
        temperature=300 * _Units.Temperature.kelvin,
        pressure=1 * _Units.Pressure.atm,
        thermostat_time_constant=1 * _Units.Time.picosecond,
        report_interval=100,
        restart_interval=100,
        core_alignment=True,
        com_distance_restraint=True,
        com_k=25 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
        com_restraint_width=5 * _Units.Length.angstrom,
        restraint=None,
        force_constant=10 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
        positional_restraint_width=0.5 * angstrom,
        align_kf_sep=2.5 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
        align_k_theta=10 * _Units.Energy.kcal_per_mol,
        align_k_psi=10 * _Units.Energy.kcal_per_mol,
        SC_umax=1000 * _Units.Energy.kcal_per_mol,
        SC_u0=500 * _Units.Energy.kcal_per_mol,
        SC_a=0.0625,
        direction=1,
        lambda1=0.0,
        lambda2=0.0,
        alpha=0.0,
        uh=0.0,
        W0=0.0,
        anneal_values="default",
        anneal_numcycles=100,
    ):
        """
        Create a new annealing protocol.

        Parameters
        ----------
        system : :class:`System <BioSimSpace._SireWrappers.System>`
            A prepared AToM system.

        data : dict
            The AToM data dictionary.

        timestep : str, :class:`Time <BioSimSpace.Types.Time>`
            The integration timestep.

        runtime : str, :class:`Time <BioSimSpace.Types.Time>`
            The running time.

        temperature : str, :class:`Temperature <BioSimSpace.Types.Temperature>`
            The temperature.

        pressure : str, :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure. Pass pressure=None to use the NVT ensemble.

        thermostat_time_constant : str, :class:`Time <BioSimSpace.Types.Time>`
            Time constant for thermostat coupling.

        report_interval : int
            The frequency at which statistics are recorded. (In integration steps.)

        restart_interval : int
            The frequency at which restart configurations and trajectory

        core_alignment : bool
            Whether to use rigid core restraints to align the two ligands.

        com_distance_restraint : bool
            Whether to use a center of mass distance restraint.
            This restraint applies to the protein/host and both ligands, and
            is used to maintain the relative positions of all of them.

        com_k : int, float, str, :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`
            The force constant for the center of mass distance restraint (kcal/mol/A^2).

        com_restraint_width : int, float, str, :class:`Length <BioSimSpace.Types.Length>
            The width (tolerance) of the center of mass distance restraint (A).

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

        force_constant : float, :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`
            The force constant for the restraint potential. If a 'float' is
            passed, then default units of 'kcal_per_mol / angstrom**2' will
            be used.

        positional_restraint_width : float, :class:`Length <BioSimSpace.Types.Length>`
            The width of the flat-bottom potential used for coordinate restraint in Angstroms.

        pos_restrained_atoms : [int]
            The atoms to be restrained.

        align_kf_sep : int, float, str, :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`
            The force constant for the distance portion of the alignment restraint (kcal/(mol A^2)).

        align_k_theta : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The force constant for the angular portion of the alignment restaint (kcal/mol).

        align_k_psi : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The force constant for the dihedral portion of the alignment restraint (kcal/mol).

        SC_umax : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The Umax value for the ATM softcore potential (kcal/mol).

        SC_u0 : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The uh value for the ATM softcore potential (kcal/mol).

        SC_a : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The a value for the ATM softcore potential.

        direction : str
            The direction of the Annealing.

        lambda1 : float
            The lambda1 value for the ATM force.
            Superceded by any values defined in anneal_values.

        lambda2 : float
            The lambda2 value for the ATM force.
            Superceded by any values defined in anneal_values.

        alpha : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The alpha value for the ATM force.
            Value in kcal/mol.
            Superceded by any values defined in anneal_values.

        uh : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The uh value for the ATM force.
            Value in kcal/mol.
            Superceded by any values defined in anneal_values.

        W0 : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The W0 value for the ATM force.
            Value in kcal/mol.
            Superceded by any values defined in anneal_values.

        anneal_values : dict, None, "default"
            If None, then no annealing will be performed.
            If "default", then lambda values will be annealed from 0 to 0.5.
            If more complex annealing is required, then
            a dictionary with some or all of the following keys should be given:
                "lambda1_start" : float
                    The starting value for lambda1.
                "lambda1_end" : float
                    The ending value for lambda1.
                "lambda2_start" : float
                    The starting value for lambda2.
                "lambda2_end" : float
                    The ending value for lambda2.
                "alpha_start" : float
                    The starting value for alpha.
                "alpha_end" : float
                    The ending value for alpha.
                "uh_start" : float
                    The starting value for uh.
                "uh_end" : float
                    The ending value for uh.
                "W0_start" : float
                    The starting value for W0.
                "W0_end" : float
                    The ending value for W0
            Any unspecified values will use their default lambda=0 value.

        anneal_numcycles : int
            The number of annealing cycles to perform, defines the rate at which values are incremented. Default 100.
        """
        super().__init__(
            system=system,
            data=data,
            core_alignment=core_alignment,
            com_distance_restraint=com_distance_restraint,
            com_k=com_k,
            com_restraint_width=com_restraint_width,
            restraint=restraint,
            force_constant=force_constant,
            positional_restraint_width=positional_restraint_width,
            align_kf_sep=align_kf_sep,
            align_k_theta=align_k_theta,
            align_k_psi=align_k_psi,
            SC_umax=SC_umax,
            SC_u0=SC_u0,
            SC_a=SC_a,
        )

        self.setTimestep(timestep)

        self.setRuntime(runtime)

        self.setTemperature(temperature)

        # Set the system pressure.
        if pressure is not None:
            self.setPressure(pressure)
        else:
            self._pressure = None

        self.setThermostatTimeConstant(thermostat_time_constant)

        self.setReportInterval(report_interval)

        self.setRestartInterval(restart_interval)

        self.setDirection(direction)

        self.setLambda1(lambda1)

        self.setLambda2(lambda2)

        self.setAlpha(alpha)

        self.setUh(uh)

        self.setW0(W0)

        # Store the anneal values.
        self.setAnnealValues(anneal_values)

        # Set the number of annealing cycles.
        self.setAnnealNumCycles(anneal_numcycles)

    def getTimeStep(self):
        """
        Return the time step.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The integration time step.
        """
        return self._timestep

    def setTimestep(self, timestep):
        """
        Set the time step.

        Parameters
        ----------

        time : str, :class:`Time <BioSimSpace.Types.Time>`
            The integration time step.
        """
        if isinstance(timestep, str):
            try:
                self._timestep = _Types.Time(timestep)
            except:
                raise ValueError("Unable to parse 'timestep' string.") from None
        elif isinstance(timestep, _Types.Time):
            self._timestep = timestep
        else:
            raise TypeError(
                "'timestep' must be of type 'str' or 'BioSimSpace.Types.Time'"
            )

    def getRunTime(self):
        """
        Return the running time.

        Returns
        -------

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The simulation run time.
        """
        return self._runtime

    def setRuntime(self, runtime):
        """
        Set the running time.

        Parameters
        ----------

        runtime : str, :class:`Time <BioSimSpace.Types.Time>`
            The simulation run time.
        """
        if isinstance(runtime, str):
            try:
                self._runtime = _Types.Time(runtime)
            except:
                raise ValueError("Unable to parse 'runtime' string.") from None
        elif isinstance(runtime, _Types.Time):
            self._runtime = runtime
        else:
            raise TypeError(
                "'runtime' must be of type 'str' or 'BioSimSpace.Types.Time'"
            )

    def getTemperature(self):
        """
        Return temperature.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The simulation temperature.
        """
        return self._temperature

    def setTemperature(self, temperature):
        """
        Set the temperature.

        Parameters
        ----------

        temperature : str, :class:`Temperature <BioSimSpace.Types.Temperature>`
            The simulation temperature.
        """
        if isinstance(temperature, str):
            try:
                self._temperature = _Types.Temperature(temperature)
            except:
                raise ValueError("Unable to parse 'temperature' string.") from None
        elif isinstance(temperature, _Types.Temperature):
            self._temperature = temperature
        else:
            raise TypeError(
                "'temperature' must be of type 'str' or 'BioSimSpace.Types.Temperature'"
            )

    def getPressure(self):
        """
        Return the pressure.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure.
        """
        return self._pressure

    def setPressure(self, pressure):
        """
        Set the pressure.

        Parameters
        ----------

        pressure : str, :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure.
        """
        if isinstance(pressure, str):
            try:
                self._pressure = _Types.Pressure(pressure)
            except:
                raise ValueError("Unable to parse 'pressure' string.") from None
        elif isinstance(pressure, _Types.Pressure):
            self._pressure = pressure
        else:
            raise TypeError(
                "'pressure' must be of type 'str' or 'BioSimSpace.Types.Pressure'"
            )

    def getThermostatTimeConstant(self):
        """
        Return the time constant for the thermostat.

        Returns
        -------

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The time constant for the thermostat.
        """
        return self._thermostat_time_constant

    def setThermostatTimeConstant(self, thermostat_time_constant):
        """
        Set the time constant for the thermostat.

        Parameters
        ----------

        thermostat_time_constant : str, :class:`Time <BioSimSpace.Types.Time>`
            The time constant for the thermostat.
        """
        if isinstance(thermostat_time_constant, str):
            try:
                self._thermostat_time_constant = _Types.Time(thermostat_time_constant)
            except:
                raise ValueError(
                    "Unable to parse 'thermostat_time_constant' string."
                ) from None
        elif isinstance(thermostat_time_constant, _Types.Time):
            self._thermostat_time_constant = thermostat_time_constant
        else:
            raise TypeError(
                "'thermostat_time_constant' must be of type 'BioSimSpace.Types.Time'"
            )

    def getReportInterval(self):
        """
        Return the interval between reporting statistics. (In integration steps.).
        Returns
        -------
        report_interval : int
            The number of integration steps between reporting statistics.
        """
        return self._report_interval

    def setReportInterval(self, report_interval):
        """
        Set the interval at which statistics are reported. (In integration steps.).

        Parameters
        ----------

        report_interval : int
            The number of integration steps between reporting statistics.
        """
        if not type(report_interval) is int:
            raise TypeError("'report_interval' must be of type 'int'")

        if report_interval <= 0:
            _warnings.warn("'report_interval' must be positive. Using default (100).")
            report_interval = 100

        self._report_interval = report_interval

    def getRestartInterval(self):
        """
        Return the interval between saving restart confiugrations, and/or
        trajectory frames. (In integration steps.).

        Returns
        -------

        restart_interval : int
            The number of integration steps between saving restart
            configurations and/or trajectory frames.
        """
        return self._restart_interval

    def setRestartInterval(self, restart_interval):
        """
        Set the interval between saving restart confiugrations, and/or
        trajectory frames. (In integration steps.).

        Parameters
        ----------

        restart_interval : int
            The number of integration steps between saving restart
            configurations and/or trajectory frames.
        """
        if not type(restart_interval) is int:
            raise TypeError("'restart_interval' must be of type 'int'")

        if restart_interval <= 0:
            _warnings.warn("'restart_interval' must be positive. Using default (500).")
            restart_interval = 500

        self._restart_interval = restart_interval

    def getDirection(self):
        """
        Return the direction of the equilibration.

        Returns
        -------

        direction : str
            The direction of the equilibration. Ignored if use_atm_force is False.
        """
        return self._direction

    def setDirection(self, direction):
        """
        Set the direction of the equilibration.

        Parameters
        ----------

        direction : str
            The direction of the equilibration. Ignored if use_atm_force is False.
        """
        if int(direction) != 1 and int(direction) != -1:
            raise TypeError("'direction' must have a value of 1 or -1")
        self._direction = int(direction)

    def getLambda1(self):
        """
        Return the lambda1 value for the ATM force.

        Returns
        -------

        lambda1 : float
            The lambda1 value for the ATM force. Ignored if use_atm_force is False.
        """
        return self._lambda1

    def setLambda1(self, lambda1):
        """
        Set the lambda1 value for the ATM force.

        Parameters
        ----------

        lambda1 : float
            The lambda1 value for the ATM force. Ignored if use_atm_force is False.
        """
        if not isinstance(lambda1, (float, int)):
            raise TypeError("'lambda1' must be of type 'float'")
        if not 0 <= float(lambda1) <= 0.5:
            raise ValueError("lambda1 must be between 0 and 0.5")
        self._lambda1 = float(lambda1)

    def getLambda2(self):
        """
        Return the lambda2 value for the ATM force.

        Returns
        -------

        lambda2 : float
            The lambda2 value for the ATM force. Ignored if use_atm_force is False.
        """
        return self._lambda2

    def setLambda2(self, lambda2):
        """
        Set the lambda2 value for the ATM force.

        Parameters
        ----------

        lambda2 : float
            The lambda2 value for the ATM force. Ignored if use_atm_force is False.
        """
        if not isinstance(lambda2, (float, int)):
            raise TypeError("'lambda2' must be of type 'float'")
        if not 0 <= float(lambda2) <= 0.5:
            raise ValueError("lambda2 must be between 0 and 0.5")
        self._lambda2 = float(lambda2)

    def getAlpha(self):
        """
        Return the alpha value for the ATM force.

        Returns
        -------

        alpha : :class:`Energy <BioSimSpace.Types.Energy>`
            The alpha value for the ATM force in kcal/mol. Ignored if use_atm_force is False.
        """
        return self._alpha

    def setAlpha(self, alpha):
        """
        Set the alpha value for the ATM force.

        Parameters
        ----------

        alpha : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The alpha value for the ATM force in kcal/mol. Ignored if use_atm_force is False.
        """
        # Convert int to float.
        if type(alpha) is int:
            alpha = float(alpha)

        if isinstance(alpha, float):
            # Use default units.
            alpha *= _Units.Energy.kcal_per_mol

        else:
            if isinstance(alpha, str):
                try:
                    alpha = _Types._GeneralUnit(alpha)
                except Exception:
                    raise ValueError("Unable to parse 'alpha' string.") from None

            elif not isinstance(alpha, _Types.Energy):
                raise TypeError(
                    "'alpha' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                )

            # Validate the dimensions.
            if alpha.dimensions() != (1, 2, -2, 0, 0, -1, 0):
                raise ValueError(
                    "'align_k_theta' has invalid dimensions! "
                    f"Expected dimensions of energy density (e.g. kcal/mol), found '{alpha.unit()}'"
                )
        self._alpha = alpha

    def getUh(self):
        """
        Return the uh value for the ATM force.

        Returns
        -------

        uh : :class:`Energy <BioSimSpace.Types.Energy>`
            The uh value for the ATM force in kcal/mol. Ignored if use_atm_force is False.
        """
        return self._uh

    def setUh(self, uh):
        """
        Set the uh value for the ATM force.

        Parameters
        ----------

        uh : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The uh value for the ATM force in kcal/mol. Ignored if use_atm_force is False.
        """
        # Convert int to float.
        if type(uh) is int:
            uh = float(uh)

        if isinstance(uh, float):
            # Use default units.
            uh *= _Units.Energy.kcal_per_mol

        else:
            if isinstance(uh, str):
                try:
                    uh = _Types._GeneralUnit(uh)
                except Exception:
                    raise ValueError("Unable to parse 'uh' string.") from None

            elif not isinstance(uh, _Types.Energy):
                raise TypeError(
                    "'uh' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                )

            # Validate the dimensions.
            if uh.dimensions() != (1, 2, -2, 0, 0, -1, 0):
                raise ValueError(
                    "'align_k_theta' has invalid dimensions! "
                    f"Expected dimensions of energy density (e.g. kcal/mol), found '{uh.unit()}'"
                )
        self._uh = uh

    def getW0(self):
        """
        Return the W0 value for the ATM force.

        Returns
        -------

        W0 : :class:`Energy <BioSimSpace.Types.Energy>`
            The W0 value for the ATM force in kcal/mol. Ignored if use_atm_force is False.
        """
        return self._W0

    def setW0(self, W0):
        """
        Set the W0 value for the ATM force.

        Parameters
        ----------

        W0 :int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The W0 value for the ATM force in kcal/mol. Ignored if use_atm_force is False.
        """
        # Convert int to float.
        if type(W0) is int:
            W0 = float(W0)

        if isinstance(W0, float):
            # Use default units.
            W0 *= _Units.Energy.kcal_per_mol

        else:
            if isinstance(W0, str):
                try:
                    W0 = _Types._GeneralUnit(W0)
                except Exception:
                    raise ValueError("Unable to parse 'W0' string.") from None

            elif not isinstance(W0, _Types.Energy):
                raise TypeError(
                    "'W0' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                )

            # Validate the dimensions.
            if W0.dimensions() != (1, 2, -2, 0, 0, -1, 0):
                raise ValueError(
                    "'align_k_theta' has invalid dimensions! "
                    f"Expected dimensions of energy density (e.g. kcal/mol), found '{W0.unit()}'"
                )
        self._W0 = W0

    def getAnnealValues(self):
        """
        Return the anneal protocol.

        Returns
        -------

        anneal_protocol : dict
            The anneal protocol.
        """
        return self._anneal_values

    def setAnnealValues(self, anneal_values):
        """
        Set the anneal protocol.

        Parameters
        ----------

        anneal_values : dict
            The anneal values.
        """

        def capitalise_keys(input_dict):
            # The first letter of each key needs to be captilised
            # so that it can be properly passed to openMM later
            capitalized_dict = {}
            for key, value in input_dict.items():
                capitalized_key = key.capitalize()
                capitalized_dict[capitalized_key] = value
            return capitalized_dict

        if anneal_values == "default":
            self._anneal_values = capitalise_keys(
                {
                    "lambda1_start": 0,
                    "lambda1_end": 0.5,
                    "lambda2_start": 0,
                    "lambda2_end": 0.5,
                }
            )
        elif isinstance(anneal_values, dict):
            # check that the given keys are valid
            keys = [
                "lambda1_start",
                "lambda1_end",
                "lambda2_start",
                "lambda2_end",
                "alpha_start",
                "alpha_end",
                "uh_start",
                "uh_end",
                "W0_start",
                "W0_end",
            ]
            for key in anneal_values:
                if key not in keys:
                    raise ValueError(
                        f"The anneal values can only contain the following keys: 'lambda1_start', 'lambda1_end', 'lambda2_start', 'lambda2_end', 'alpha_start', 'alpha_end', 'uh_start', 'uh_end', 'W0_start', 'W0_end', 'runtime'. The following keys are invalid: {key}"
                    )
                if key == "lambda1_start" or key == "lambda1_end":
                    if not 0 <= float(anneal_values[key]) <= 0.5:
                        raise ValueError("lambda1 must be between 0 and 0.5")
                if key == "lambda2_start" or key == "lambda2_end":
                    if not 0 <= float(anneal_values[key]) <= 0.5:
                        raise ValueError("lambda2 must be between 0 and 0.5")
                # check that none of the other keys are negative
                if (
                    key != "lambda1_start"
                    and key != "lambda1_end"
                    and key != "lambda2_start"
                    and key != "lambda2_end"
                ):
                    if float(anneal_values[key]) < 0:
                        raise ValueError(f"{key} must be greater than or equal to 0")
                # also check that they are floats
                if not isinstance(anneal_values[key], (float, int)):
                    raise TypeError(f"{key} must be of type 'float'")
            self._anneal_values = capitalise_keys(anneal_values)
        elif anneal_values is None:
            self._anneal_values = None

        else:
            raise TypeError(
                "'anneal_values' must be of type 'dict', 'None', or 'default'"
            )

    def getAnnealNumCycles(self):
        """
        Return the number of annealing cycles.

        Returns
        -------

        anneal_numcycles : int
            The number of annealing cycles.
        """
        return self._anneal_numcycles

    def setAnnealNumCycles(self, anneal_numcycles):
        """
        Set the number of annealing cycles.

        Parameters
        ----------

        anneal_numcycles : int
            The number of annealing cycles.
        """
        if isinstance(anneal_numcycles, int):
            self._anneal_numcycles = anneal_numcycles
        else:
            raise TypeError("'anneal_numcycles' must be of type 'int'")

    def _set_current_index(self, index):
        """
        The current index of the window.
        In annealing protocols this should not be touched by the user.

        Parameters
        ----------
        index : int
            The index of the current lambda window.
        """
        if index < 0:
            raise ValueError("index must be positive")
        if not isinstance(index, int):
            raise TypeError("index must be an integer")
        self._current_index = index

    def _get_window_index(self):
        """
        A function to get the index of the current lambda window.

        Returns
        -------
        index : int
            The index of the current lambda window.
        """
        try:
            return self._current_index
        except:
            return None


class AToMProduction(_AToM):
    """Production protocol for AToM simulations."""

    def __init__(
        self,
        system=None,
        data=None,
        timestep=2 * _Units.Time.femtosecond,
        runtime=1.0 * _Units.Time.nanosecond,
        temperature=300 * _Units.Temperature.kelvin,
        pressure=1 * _Units.Pressure.atm,
        thermostat_time_constant=1 * _Units.Time.picosecond,
        report_interval=100,
        restart_interval=100,
        restart=False,
        core_alignment=True,
        com_distance_restraint=True,
        com_k=25 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
        com_restraint_width=5 * _Units.Length.angstrom,
        restraint=None,
        force_constant=10 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
        positional_restraint_width=0.5 * angstrom,
        num_lambda=22,
        direction=None,
        lambda1=None,
        lambda2=None,
        alpha=None,
        uh=None,
        W0=None,
        align_kf_sep=2.5 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
        align_k_theta=10 * _Units.Energy.kcal_per_mol,
        align_k_psi=10 * _Units.Energy.kcal_per_mol,
        SC_umax=100 * _Units.Energy.kcal_per_mol,
        SC_u0=50 * _Units.Energy.kcal_per_mol,
        SC_a=0.0625,
        analysis_method="UWHAM",
    ):
        """
        Create a new production protocol.

        Parameters
        ----------
        system : :class:`System <BioSimSpace._SireWrappers.System>`
            A prepared AToM system.

        data : dict
            The AToM data dictionary.

        timestep : str, :class:`Time <BioSimSpace.Types.Time>`
            The integration timestep.

        runtime : str, :class:`Time <BioSimSpace.Types.Time>`
            The running time.

        temperature : str, :class:`Temperature <BioSimSpace.Types.Temperature>`
            The temperature.

        pressure : str, :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure. Pass pressure=None to use the NVT ensemble.

        thermostat_time_constant : str, :class:`Time <BioSimSpace.Types.Time>`
            Time constant for thermostat coupling.

        report_interval : int
            The frequency at which statistics are recorded. (In integration steps.)

        restart_interval : int
            The frequency at which restart configurations and trajectory

        core_alignment : bool
            Whether to use rigid core restraints to align the two ligands.

        com_distance_restraint : bool
            Whether to use a center of mass distance restraint.
            This restraint applies to the protein/host and both ligands, and
            is used to maintain the relative positions of all of them.

        com_k : int, float, str, :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`
            The force constant for the center of mass distance restraint (kcal/mol/A^2).

        com_restraint_width : int, float, str, :class:`Length <BioSimSpace.Types.Length>
            The width (tolerance) of the center of mass distance restraint (A).

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

        force_constant : float, :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`
            The force constant for the restraint potential. If a 'float' is
            passed, then default units of 'kcal_per_mol / angstrom**2' will
            be used.

        positional_restraint_width : :class:`Length <BioSimSpace.Types.Length>`, float
            The width of the flat-bottom potential used for coordinate restraint in Angstroms.

        pos_restrained_atoms : [int]
            The atoms to be restrained.

        align_kf_sep : int, float, str, :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>`
            The force constant for the distance portion of the alignment restraint (kcal/(mol A^2)).

        align_k_theta : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The force constant for the angular portion of the alignment restaint (kcal/mol).

        align_k_psi : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The force constant for the dihedral portion of the alignment restraint (kcal/mol).

        SC_umax : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The Umax value for the ATM softcore potential (kcal/mol).

        SC_u0 : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The uh value for the ATM softcore potential (kcal/mol).

        SC_a : int, float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The a value for the ATM softcore potential.

        restart : bool
            Whether this is a continuation of a previous simulation.

        num_lambda : int
            The number of lambda values. This will be used to set the window-dependent
            AToM parameters, unless they are explicitly set by the user.

        lambdas : [float]
            The lambda values.

        direction : [int]
            The direction values. Must be either 1 (forwards) or -1 (backwards).

        lambda1 : [float]
            The lambda1 values.

        lambda2 : [float]
            The lambda2 values.

        alpha : [int], float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The alpha values.

        uh : [int], float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The uh values.

        W0 : [int], float, str, :class:`Energy <BioSimSpace.Types.Energy>`
            The W0 values.

        analysis_method : str
            The method to use for analysis. Options are "UWHAM", "MBAR" or "both"
            This affects the output files and the analysis that is performed.
            USE of "UWHAM" is strongly recommended, "MBAR" analysis is still experimental.
        """
        super().__init__(
            system=system,
            data=data,
            core_alignment=core_alignment,
            com_distance_restraint=com_distance_restraint,
            com_k=com_k,
            com_restraint_width=com_restraint_width,
            restraint=restraint,
            force_constant=force_constant,
            positional_restraint_width=positional_restraint_width,
            align_kf_sep=align_kf_sep,
            align_k_theta=align_k_theta,
            align_k_psi=align_k_psi,
            SC_umax=SC_umax,
            SC_u0=SC_u0,
            SC_a=SC_a,
        )

        self.setTimestep(timestep)

        self.setRuntime(runtime)

        self.setTemperature(temperature)

        # Set the system pressure.
        if pressure is not None:
            self.setPressure(pressure)
        else:
            self._pressure = None

        self.setThermostatTimeConstant(thermostat_time_constant)

        self.setReportInterval(report_interval)

        self.setRestartInterval(restart_interval)

        # Set the restart flag.
        self.setRestart(restart)
        # Set the number of lambda values.
        # If other window-dependent parameters are not set, then set them to
        # sensible defaults.
        self.setNumLambda(num_lambda)

        # Store the direction values.
        self.setDirection(direction)

        # Store the lambda1 values.
        self.setLambda1(lambda1)

        # Store the lambda2 values.
        self.setLambda2(lambda2)

        # Store the alpha values.
        self.setAlpha(alpha)

        # Store the uh values.
        self.setUh(uh)

        # Store the W0 values.
        self.setW0(W0)

        self._set_lambda_values()

        self.setAnalysisMethod(analysis_method)

    def getTimeStep(self):
        """
        Return the time step.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The integration time step.
        """
        return self._timestep

    def setTimestep(self, timestep):
        """
        Set the time step.

        Parameters
        ----------

        time : str, :class:`Time <BioSimSpace.Types.Time>`
            The integration time step.
        """
        if isinstance(timestep, str):
            try:
                self._timestep = _Types.Time(timestep)
            except:
                raise ValueError("Unable to parse 'timestep' string.") from None
        elif isinstance(timestep, _Types.Time):
            self._timestep = timestep
        else:
            raise TypeError(
                "'timestep' must be of type 'str' or 'BioSimSpace.Types.Time'"
            )

    def getRunTime(self):
        """
        Return the running time. Set the same as other OpenMM protocols - really should be Runtime not RunTime.

        Returns
        -------

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The simulation run time.
        """
        return self._runtime

    def setRuntime(self, runtime):
        """
        Set the running time.

        Parameters
        ----------

        runtime : str, :class:`Time <BioSimSpace.Types.Time>`
            The simulation run time.
        """
        if isinstance(runtime, str):
            try:
                self._runtime = _Types.Time(runtime)
            except:
                raise ValueError("Unable to parse 'runtime' string.") from None
        elif isinstance(runtime, _Types.Time):
            self._runtime = runtime
        else:
            raise TypeError(
                "'runtime' must be of type 'str' or 'BioSimSpace.Types.Time'"
            )

    def getTemperature(self):
        """
        Return temperature.

        Returns
        -------

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The simulation temperature.
        """
        return self._temperature

    def setTemperature(self, temperature):
        """
        Set the temperature.

        Parameters
        ----------

        temperature : str, :class:`Temperature <BioSimSpace.Types.Temperature>`
            The simulation temperature.
        """
        if isinstance(temperature, str):
            try:
                self._temperature = _Types.Temperature(temperature)
            except:
                raise ValueError("Unable to parse 'temperature' string.") from None
        elif isinstance(temperature, _Types.Temperature):
            self._temperature = temperature
        else:
            raise TypeError(
                "'temperature' must be of type 'str' or 'BioSimSpace.Types.Temperature'"
            )

    def getPressure(self):
        """
        Return the pressure.

        Returns
        -------

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure.
        """
        return self._pressure

    def setPressure(self, pressure):
        """
        Set the pressure.

        Parameters
        ----------

        pressure : str, :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure.
        """
        if isinstance(pressure, str):
            try:
                self._pressure = _Types.Pressure(pressure)
            except:
                raise ValueError("Unable to parse 'pressure' string.") from None
        elif isinstance(pressure, _Types.Pressure):
            self._pressure = pressure
        else:
            raise TypeError(
                "'pressure' must be of type 'str' or 'BioSimSpace.Types.Pressure'"
            )

    def getThermostatTimeConstant(self):
        """
        Return the time constant for the thermostat.

        Returns
        -------

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The time constant for the thermostat.
        """
        return self._thermostat_time_constant

    def setThermostatTimeConstant(self, thermostat_time_constant):
        """
        Set the time constant for the thermostat.

        Parameters
        ----------

        thermostat_time_constant : str, :class:`Time <BioSimSpace.Types.Time>`
            The time constant for the thermostat.
        """
        if isinstance(thermostat_time_constant, str):
            try:
                self._thermostat_time_constant = _Types.Time(thermostat_time_constant)
            except:
                raise ValueError(
                    "Unable to parse 'thermostat_time_constant' string."
                ) from None
        elif isinstance(thermostat_time_constant, _Types.Time):
            self._thermostat_time_constant = thermostat_time_constant
        else:
            raise TypeError(
                "'thermostat_time_constant' must be of type 'BioSimSpace.Types.Time'"
            )

    def getReportInterval(self):
        """
        Return the interval between reporting statistics. (In integration steps.).
        Returns
        -------
        report_interval : int
            The number of integration steps between reporting statistics.
        """
        return self._report_interval

    def setReportInterval(self, report_interval):
        """
        Set the interval at which statistics are reported. (In integration steps.).

        Parameters
        ----------

        report_interval : int
            The number of integration steps between reporting statistics.
        """
        if not type(report_interval) is int:
            raise TypeError("'report_interval' must be of type 'int'")

        if report_interval <= 0:
            _warnings.warn("'report_interval' must be positive. Using default (100).")
            report_interval = 100

        self._report_interval = report_interval

    def getRestartInterval(self):
        """
        Return the interval between saving restart confiugrations, and/or
        trajectory frames. (In integration steps.).

        Returns
        -------

        restart_interval : int
            The number of integration steps between saving restart
            configurations and/or trajectory frames.
        """
        return self._restart_interval

    def setRestartInterval(self, restart_interval):
        """
        Set the interval between saving restart confiugrations, and/or
        trajectory frames. (In integration steps.).

        Parameters
        ----------

        restart_interval : int
            The number of integration steps between saving restart
            configurations and/or trajectory frames.
        """
        if not type(restart_interval) is int:
            raise TypeError("'restart_interval' must be of type 'int'")

        if restart_interval <= 0:
            _warnings.warn("'restart_interval' must be positive. Using default (500).")
            restart_interval = 500

        self._restart_interval = restart_interval

    def isRestart(self):
        """
        Return whether this restart simulation.

        Returns
        -------

        is_restart : bool
            Whether this is a restart simulation.
        """
        return self._restart

    def setRestart(self, restart):
        """
        Set the restart flag.

        Parameters
        ----------

        restart : bool
            Whether this is a restart simulation.
        """
        if isinstance(restart, bool):
            self._restart = restart
        else:
            _warnings.warn("Non-boolean restart flag. Defaulting to False!")
            self._restart = False

    def getNumLambda(self):
        """
        Return the number of lambda values.

        Returns
        -------

        num_lambda : int
            The number of lambda values.
        """
        return self._num_lambda

    def setNumLambda(self, num_lambda):
        """
        Set the number of lambda values.

        Parameters
        ----------

        num_lambda : int
            The number of lambda values.
        """
        if isinstance(num_lambda, int) and num_lambda > 0:
            if num_lambda % 2 != 0:
                _warnings.warn(
                    "Warning: The AToM protocol is optimised for an even number of lambda values. Unknown behaviour may occur if using an odd number of lambda values."
                )
            self._num_lambda = num_lambda
            self._set_lambda_values()
        else:
            raise TypeError("'num_lambda' must be of type 'int'")

    def getDirection(self):
        """
        Return the direction values.

        Returns
        -------

        lambdas : [float]
            The directions.
        """
        return self._directions

    def setDirection(self, directions):
        """
        Set the direction values.

        Parameters
        ----------

        directions : [int]
            The directions.
        """
        if isinstance(directions, list):
            if len(directions) != self._num_lambda:
                raise ValueError(
                    "'directions' must have the same length as 'num_lambda'"
                )
            if all(item == 1 or item == -1 for item in directions):
                self._directions = directions
            else:
                raise ValueError("all entries in 'directions' must be either 1 or -1")
        elif directions is None:
            self._directions = [1] * _math.floor(self._num_lambda / 2) + [
                -1
            ] * _math.ceil(self._num_lambda / 2)
        else:
            raise TypeError("'directions' must be of type 'list' or 'None'")

    def getLambda1(self):
        """
        Return the lambda1 values.

        Returns
        -------

        lambda1 : [float]
            The lambda1 values.
        """
        return self._lambda1

    def setLambda1(self, lambda1):
        """
        Set the lambda1 values.

        Parameters
        ----------

        lambda1 : [float]
            The lambda1 values.
        """
        if isinstance(lambda1, list):
            if len(lambda1) != self._num_lambda:
                raise ValueError("'lambda1' must have the same length as 'num_lambda'")
            if all(isinstance(item, float) for item in lambda1) and all(
                item <= 0.5 for item in lambda1
            ):
                self._lambda1 = lambda1
            else:
                raise ValueError(
                    "all entries in 'lambda1' must be floats with a value less than or equal to 0.5"
                )
        elif lambda1 is None:
            # use numpy to create a [float]s
            self._lambda1 = _np.concatenate(
                [
                    _np.linspace(0, 0.5, _math.floor(self._num_lambda / 2)),
                    _np.linspace(0.5, 0, _math.ceil(self._num_lambda / 2)),
                ]
            ).tolist()
            # Round the floats to 5 decimal places
            self._lambda1 = [round(num, 5) for num in self._lambda1]
        else:
            raise TypeError("'lambda1' must be of type 'list'")

    def getLambda2(self):
        """
        Return the lambda2 values.

        Returns
        -------

        lambda2 : [float]
            The lambda2 values.
        """
        return self._lambda2

    def setLambda2(self, lambda2):
        """
        Set the lambda2 values.

        Parameters
        ----------

        lambda2 : [float]
            The lambda2 values.
        """
        if isinstance(lambda2, list):
            if len(lambda2) != self._num_lambda:
                raise ValueError("'lambda2' must have the same length as 'num_lambda'")
            if all(isinstance(item, float) for item in lambda2) and all(
                item <= 0.5 for item in lambda2
            ):
                if len(lambda2) != len(self._lambda1):
                    raise ValueError(
                        "'lambda2' and 'lambda1' must have the same length"
                    )
                self._lambda2 = lambda2
            else:
                raise ValueError("all entries in 'lambda2' must be floats")
        elif lambda2 is None:
            # use numpy to create a [float]s
            self._lambda2 = _np.concatenate(
                [
                    _np.linspace(0, 0.5, _math.floor(self._num_lambda / 2)),
                    _np.linspace(0.5, 0, _math.ceil(self._num_lambda / 2)),
                ]
            ).tolist()
            # Round the floats to 5 decimal places
            self._lambda2 = [round(num, 5) for num in self._lambda2]
        else:
            raise TypeError("'lambda2' must be of type 'list'")

    def getAlpha(self):
        """
        Return the alpha values.

        Returns
        -------

        alpha : [:class:`Energy <BioSimSpace.Types.Energy>]
            The alpha values in kcal/mol.
        """
        return self._alpha

    def setAlpha(self, alpha):
        """
        Set the alpha values.

        Parameters
        ----------

        alpha : [`Energy <BioSimSpace.Types.Energy>] or [int], [float], [str]
            The alpha values in kcal/mol.
        """
        if isinstance(alpha, list):
            if len(alpha) != self._num_lambda:
                raise ValueError("'alpha' must have the same length as 'num_lambda'")
            alpha_fin = []
            for a in alpha:
                # Convert int to float.
                if type(a) is int:
                    a = float(a)
                    a *= _Units.Energy.kcal_per_mol

                elif isinstance(a, float):
                    # Use default units.
                    a *= _Units.Energy.kcal_per_mol

                else:
                    if isinstance(a, str):
                        try:
                            a = _Types._GeneralUnit(a)
                        except Exception:
                            raise ValueError(
                                "Unable to parse 'alpha' string."
                            ) from None

                    elif not isinstance(a, _Types.Energy):
                        raise TypeError(
                            "'alpha' must be of type 'BioSimSpace.Types.Energy', 'str', or 'float'."
                        )

                    # Validate the dimensions.
                    if a.dimensions() != (1, 2, -2, 0, 0, -1, 0):
                        raise ValueError(
                            "'alpha' has invalid dimensions! "
                            f"Expected dimensions of energy density (e.g. kcal/mol), found '{a.unit()}'"
                        )
                alpha_fin.append(a)
            self._alpha = alpha_fin
        elif alpha is None:
            self._alpha = [0.00 * _Units.Energy.kcal_per_mol] * self._num_lambda
        else:
            raise TypeError("'alpha' must be of type 'list' or None")

    def getUh(self):
        """
        Return the uh values.

        Returns
        -------

        uh : [:class:`Energy <BioSimSpace.Types.Energy>]
            The uh values in kcal/mol.
        """
        return self._uh

    def setUh(self, uh):
        """
        Set the uh values.

        Parameters
        ----------

        uh : [:class:`Energy <BioSimSpace.Types.Energy>]
            The uh values in kcal/mol.
        """
        if isinstance(uh, list):
            if len(uh) != self._num_lambda:
                raise ValueError("'uh' must have the same length as 'num_lambda'")
            uh_fin = []
            for u in uh:
                # Convert int to float.
                if type(u) is int:
                    u = float(u)
                    u *= _Units.Energy.kcal_per_mol

                if isinstance(u, float):
                    # Use default units.
                    u *= _Units.Energy.kcal_per_mol

                else:
                    if isinstance(u, str):
                        try:
                            u = _Types._GeneralUnit(u)
                        except Exception:
                            raise ValueError(
                                "Unable to parse 'alpha' string."
                            ) from None

                    elif not isinstance(u, _Types.Energy):
                        raise TypeError(
                            "'alpha' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                        )

                    # Validate the dimensions.
                    if u.dimensions() != (1, 2, -2, 0, 0, -1, 0):
                        raise ValueError(
                            "'alpha' has invalid dimensions! "
                            f"Expected dimensions of energy density (e.g. kcal/mol), found '{u.unit()}'"
                        )
                uh_fin.append(u)
            self._uh = uh_fin
        elif uh is None:
            self._uh = [0.00 * _Units.Energy.kcal_per_mol] * self._num_lambda
        else:
            raise TypeError("'uh' must be of type 'list'")

    def getW0(self):
        """
        Return the W0 values.

        Returns
        -------

        W0 : [:class:`Energy <BioSimSpace.Types.Energy>]
            The W0 values in kcal/mol.
        """
        return self._W0

    def setW0(self, W0):
        """
        Set the W0 values.

        Parameters
        ----------

        W0 : [:class:`Energy <BioSimSpace.Types.Energy>] or [int], [float], [str]
            The W0 values in kcal/mol.
        """
        if isinstance(W0, list):
            if len(W0) != self._num_lambda:
                raise ValueError("'W0' must have the same length as 'num_lambda'")
            W0_fin = []
            for w in W0:
                # Convert int to float.
                if type(w) is int:
                    w = float(w)
                    w *= _Units.Energy.kcal_per_mol

                if isinstance(w, float):
                    # Use default units.
                    w *= _Units.Energy.kcal_per_mol

                else:
                    if isinstance(w, str):
                        try:
                            w = _Types._GeneralUnit(w)
                        except Exception:
                            raise ValueError(
                                "Unable to parse 'alpha' string."
                            ) from None

                    elif not isinstance(w, _Types.Energy):
                        raise TypeError(
                            "'alpha' must be of type 'BioSimSpace.Types._GeneralUnit', 'str', or 'float'."
                        )

                    # Validate the dimensions.
                    if w.dimensions() != (1, 2, -2, 0, 0, -1, 0):
                        raise ValueError(
                            "'alpha' has invalid dimensions! "
                            f"Expected dimensions of energy density (e.g. kcal/mol), found '{w.unit()}'"
                        )
                W0_fin.append(w)
            self._W0 = W0_fin
        elif W0 is None:
            self._W0 = [0.00 * _Units.Energy.kcal_per_mol] * self._num_lambda
        else:
            raise TypeError("'W0' must be of type 'list'")

    def _set_lambda_values(self):
        # Internal function to set the 'master lambda'
        # This lambda value serves as the master for all other window-dependent parameters
        self._lambda_values = _np.linspace(0, 1, self._num_lambda).tolist()

    def _get_lambda_values(self):
        # Internal function to get the 'master lambda'
        # This lambda value serves as the master for all other window-dependent parameters
        try:
            return self._lambda_values
        except:
            return None

    def setAnalysisMethod(self, analysis_method):
        """Set the method that will be used for analysis of the simulation results.
        This will change the output files that are generated.

        Parameters
        ----------
        analysis_method : str
            The method to use for analysis. Options are "UWHAM", "MBAR" or "both"
            This affects the output files and the analysis that is performed.
            USE of "UWHAM" is strongly recommended, "MBAR" analysis is still experimental.
        """
        allowed_methods = ["UWHAM", "MBAR", "both"]
        if analysis_method in allowed_methods:
            self._analysis_method = analysis_method
        else:
            raise ValueError(f"analysis_method must be one of {allowed_methods}")

    def getAnalysisMethod(self):
        return self._analysis_method

    def set_current_index(self, index):
        """
        A function to set the index of the current lambda window.
        Used internally to set the values for all lambda-dependent parameters.
        Take care when using this function as it can lead to unexpected behaviour if not used correctly.

        Parameters
        ----------
        index : int
            The index of the current lambda window.
        """
        if index < 0:
            raise ValueError("index must be positive")
        if index >= len(self._lambda1):
            raise ValueError(
                "index must be less than the number of lambda1 values (len(lambda1))"
            )
        if not isinstance(index, int):
            raise TypeError("index must be an integer")
        self._current_index = index

    def get_window_index(self):
        """
        A function to get the index of the current lambda window.

        Returns
        -------
        index : int
            The index of the current lambda window.
        """
        try:
            return self._current_index
        except:
            return None
