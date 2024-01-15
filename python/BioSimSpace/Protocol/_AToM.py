from .. import Types as _Types
from ._protocol import Protocol as _Protocol
from ._position_restraint_mixin import _PositionRestraintMixin
from .. import Units as _Units
import math as _math
import numpy as _np

__all__ = ["AToM"]


# When placed in to BSS this needs to be AToM_protocol(protocol):
class AToM(_Protocol, _PositionRestraintMixin):
    def __init__(
        self,
        data,
        timestep=_Types.Time(2.0, "femtosecond"),
        runtime=_Types.Time(1, "picosecond"),
        temperature=_Types.Temperature(300, "kelvin"),
        pressure=_Types.Pressure(1, "atmosphere"),
        thermostat_time_constant=_Types.Time(1.0, "picosecond"),
        report_interval=100,
        restart_interval=100,
        restart=False,
        core_alignment=True,
        restraint=None,
        force_constant=10 * _Units.Energy.kcal_per_mol / _Units.Area.angstrom2,
        num_lambda=22,
        directions=None,
        lambda1=None,
        lambda2=None,
        alpha=None,
        U0=None,
        W0=None,
        align_kf_sep=2.5,
        align_k_theta=10.0,
        align_k_psi=10.0,
        SC_umax=100.0,
        SC_u0=50.0,
        sc_a=0.0625,
    ):
        """
        Create a protocol object.

        Parameters
        ----------

        data : dict
            The AToM data dictionary.

        timestep : :class:`Time <BioSimSpace.Types.Time>`
            The integration timestep.

        runtime : :class:`Time <BioSimSpace.Types.Time>`
            The running time.

        temperature : :class:`Temperature <BioSimSpace.Types.Temperature>`
            The temperature.

        pressure : :class:`Pressure <BioSimSpace.Types.Pressure>`
            The pressure. Pass pressure=None to use the NVT ensemble.

        thermostat_time_constant : :class:`Time <BioSimSpace.Types.Time>`
            Time constant for thermostat coupling.

        report_interval : int
            The frequency at which statistics are recorded. (In integration steps.)

        restart_interval : int
            The frequency at which restart configurations and trajectory

        restart : bool
            Whether this is a continuation of a previous simulation.

        use_core_alignment : bool
            Whether to use rigid core restraints to align the two ligands.

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
            be used.Restraint <BioSimSpace.Types.Restraint>`
            The restraint.

        num_lambda : int
            The number of lambda values. This will be used to set the window-dependent
            AToM parameters, unless they are explicitly set by the user.

        lambdas : list of float
            The lambda values.

        direction : list of int
            The direction values. Must be either 1 (forwards) or -1 (backwards).

        lambda1 : list of float
            The lambda1 values.

        lambda2 : list of float
            The lambda2 values.

        alpha : list of float
            The alpha values.

        U0 : list of float
            The U0 values.

        W0 : list of float
            The W0 values.

        pos_restrained_atoms : list of int
            The atoms to be restrained.

        align_kf_sep : float
            The force constant for the distance portion of the alignment restraint (kcal/mol/A^2).

        align_k_theta : float
            The force constant for the angular portion of the alignment restaint (kcal/mol).

        align_k_psi : float
            The force constant for the dihedral portion of the alignment restraint (kcal/mol).

        SC_umax : float
            The Umax value for the ATM softcore potential (kcal/mol).

        SC_u0 : float
            The U0 value for the ATM softcore potential (kcal/mol).

        sc_a : float
            The a value for the ATM softcore potential.

        """
        # TODO Make the master num_lambda functional
        # i.e use it to set other arrays of values to some sensible numbers
        # more complex than it may appear due to the spread of values
        # Call the base class constructor.
        super().__init__()

        # Store the AToM system.
        if isinstance(data, dict):
            self._system_data = data
        else:
            raise TypeError("'data' must be of type 'dict'")

        self.setTimeStep(timestep)

        self.setRunTime(runtime)

        self.setTemperature(temperature)

        # Set the system pressure.
        if pressure is not None:
            self.setPressure(pressure)
        else:
            self._pressure = None

        self.setThermostatTimeConstant(thermostat_time_constant)

        self.setReportInterval(report_interval)

        self.setRestartInterval(restart_interval)

        self.setRestart(restart)

        # Whether or not to use alignment restraints.
        self.setCoreAlignment(core_alignment)

        # Set the number of lambda values.
        # If other window-dependent parameters are not set, then set them to
        # sensible defaults.
        self.setNumLambda(num_lambda)

        # Store the direction values.
        self.setDirections(directions)

        # Store the lambda1 values.
        self.setLambda1(lambda1)

        # Store the lambda2 values.
        self.setLambda2(lambda2)

        # Store the alpha values.
        self.setAlpha(alpha)

        # Store the U0 values.
        self.setU0(U0)

        # Store the W0 values.
        self.setW0(W0)

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

        # Store the sc_a value.
        self.setSCa(sc_a)

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

    def getTimeStep(self):
        """
        Return the time step.

        Returns
        -------

        time : :class:`Time <BioSimSpace.Types.Time>`
            The integration time step.
        """
        return self._timestep

    def setTimeStep(self, timestep):
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

    def setRunTime(self, runtime):
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
            print("'report_interval' must be positive. Using default (100).")
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
            print("'restart_interval' must be positive. Using default (500).")
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
            print("Non-boolean restart flag. Defaulting to False!")
            self._restart = False

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
            print("Non-boolean core alignment flag. Defaulting to True!")
            self._core_alignment = True

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
                # Swap this print statement with a warning
                print(
                    "Warning: The AToM protocol is optimised for an even number of lambda values. Unknown behaviour may occur if using an odd number of lambda values."
                )
            self._num_lambda = num_lambda
            self._set_lambda_values()
        else:
            raise TypeError("'num_lambda' must be of type 'int'")

    def getDirections(self):
        """
        Return the direction values.

        Returns
        -------

        lambdas : list of float
            The directions.
        """
        return self._directions

    def setDirections(self, directions):
        """
        Set the direction values.

        Parameters
        ----------

        directions : list of int
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

        lambda1 : list of float
            The lambda1 values.
        """
        return self._lambda1

    def setLambda1(self, lambda1):
        """
        Set the lambda1 values.

        Parameters
        ----------

        lambda1 : list of float
            The lambda1 values.
        """
        if isinstance(lambda1, list):
            if len(lambda1) != self._num_lambda:
                raise ValueError("'lambda1' must have the same length as 'num_lambda'")
            if all(isinstance(item, float) for item in lambda1):
                self._lambda1 = lambda1
            else:
                raise ValueError("all entries in 'lambda1' must be floats")
        elif lambda1 is None:
            # use numpy to create a list of floats
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

        lambda2 : list of float
            The lambda2 values.
        """
        return self._lambda2

    def setLambda2(self, lambda2):
        """
        Set the lambda2 values.

        Parameters
        ----------

        lambda2 : list of float
            The lambda2 values.
        """
        if isinstance(lambda2, list):
            if len(lambda2) != self._num_lambda:
                raise ValueError("'lambda2' must have the same length as 'num_lambda'")
            if all(isinstance(item, float) for item in lambda2):
                if len(lambda2) != len(self._lambda1):
                    raise ValueError(
                        "'lambda2' and 'lambda1' must have the same length"
                    )
                self._lambda2 = lambda2
            else:
                raise ValueError("all entries in 'lambda2' must be floats")
        elif lambda2 is None:
            # use numpy to create a list of floats
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

        alpha : list of float
            The alpha values.
        """
        return self._alpha

    def setAlpha(self, alpha):
        """
        Set the alpha values.

        Parameters
        ----------

        alpha : list of float
            The alpha values.
        """
        if isinstance(alpha, list):
            if len(alpha) != self._num_lambda:
                raise ValueError("'alpha' must have the same length as 'num_lambda'")
            if all(isinstance(item, float) for item in alpha):
                self._alpha = alpha
            else:
                raise ValueError("all entries in 'alpha' must be floats")
        elif alpha is None:
            self._alpha = [0.00] * self._num_lambda
        else:
            raise TypeError("'alpha' must be of type 'list'")

    def getU0(self):
        """
        Return the U0 values.

        Returns
        -------

        U0 : list of float
            The U0 values.
        """
        return self._U0

    def setU0(self, U0):
        """
        Set the U0 values.

        Parameters
        ----------

        U0 : list of float
            The U0 values.
        """
        if isinstance(U0, list):
            if len(U0) != self._num_lambda:
                raise ValueError("'U0' must have the same length as 'num_lambda'")
            if all(isinstance(item, float) for item in U0):
                self._U0 = U0
            else:
                raise ValueError("all entries in 'U0' must be floats")
        elif U0 is None:
            self._U0 = [0.00] * self._num_lambda
        else:
            raise TypeError("'U0' must be of type 'list'")

    def getW0(self):
        """
        Return the W0 values.

        Returns
        -------

        W0 : list of float
            The W0 values.
        """
        return self._W0

    def setW0(self, W0):
        """
        Set the W0 values.

        Parameters
        ----------

        W0 : list of float
            The W0 values.
        """
        if isinstance(W0, list):
            if len(W0) != self._num_lambda:
                raise ValueError("'W0' must have the same length as 'num_lambda'")
            if all(isinstance(item, float) for item in W0):
                self._W0 = W0
            else:
                raise ValueError("all entries in 'W0' must be floats")
        elif W0 is None:
            self._W0 = [0.00] * self._num_lambda
        else:
            raise TypeError("'W0' must be of type 'list'")

    def getAlignKfSep(self):
        """
        Return the align_kf_sep value.

        Returns
        -------

        align_kf_sep : str, BSS.Types.Length
            The align_kf_sep value.
        """
        return self._align_kf_sep

    def setAlignKfSep(self, align_kf_sep):
        """
        Set the align_kf_sep value.

        Parameters
        ----------

        align_kf_sep : str, BSS.Types.Length
            Length value for the alignment restraint.
        """
        if isinstance(align_kf_sep, float):
            self._align_kf_sep = align_kf_sep
        else:
            raise TypeError("'align_kf_sep' must be of type 'float'")

    def getAlignKTheta(self):
        """
        Return the align_k_theta value.

        Returns
        -------

        align_k_theta : str, BSS.Types.Angle
            The align_k_theta value.
        """
        return self._align_k_theta

    def setAlignKTheta(self, align_k_theta):
        """
        Set the align_k_theta value.

        Parameters
        ----------

        align_k_theta : str, BSS.Types.Angle
            Angle value for the alignment angular constraint.
        """
        if isinstance(align_k_theta, float):
            self._align_k_theta = align_k_theta
        else:
            raise TypeError("'align_k_theta' must be of type 'float'")

    def getAlignKPsi(self):
        """
        Return the align_k_psi value.

        Returns
        -------

        align_k_psi : str, BSS.Types.Angle
            The align_k_psi value.
        """
        return self._align_k_psi

    def setAlignKPsi(self, align_k_psi):
        """
        Set the align_k_psi value.

        Parameters
        ----------

        align_k_psi : float
            Angle value for the alignment dihedral constraint.
        """
        if isinstance(align_k_psi, float):
            self._align_k_psi = align_k_psi
        else:
            raise TypeError("'align_k_psi' must be of type 'float'")

    def getSCUmax(self):
        """
        Return the SC_umax value.

        Returns
        -------

        SC_umax : float
            The SC_umax value.
        """
        return self._SC_umax

    def setSCUmax(self, SC_umax):
        """
        Set the SC_umax value.

        Parameters
        ----------

        SC_umax : float
            The SC_umax value.
        """
        if isinstance(SC_umax, float):
            self._SC_umax = SC_umax
        else:
            raise TypeError("'SC_umax' must be of type 'float'")

    def getSCU0(self):
        """
        Return the SC_u0 value.

        Returns
        -------

        SC_u0 : float
            The SC_u0 value.
        """
        return self._SC_u0

    def setSCU0(self, SC_u0):
        """
        Set the SC_u0 value.

        Parameters
        ----------

        SC_u0 : float
            The SC_u0 value.
        """
        if isinstance(SC_u0, float):
            self._SC_u0 = SC_u0
        else:
            raise TypeError("'SC_u0' must be of type 'float'")

    def getSCa(self):
        """
        Return the sc_a value.

        Returns
        -------

        sc_a : float
            The sc_a value.
        """
        return self._sc_a

    def setSCa(self, sc_a):
        """
        Set the sc_a value.

        Parameters
        ----------

        sc_a : float
            The sc_a value.
        """
        if isinstance(sc_a, float):
            self._sc_a = sc_a
        else:
            raise TypeError("'sc_a' must be of type 'float'")

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

    def _set_current_index(self, index):
        # Internal function to set index of the current simulation window
        # set using the master lambda list
        if index < 0:
            raise ValueError("index must be positive")
        if index >= len(self._lambda1):
            raise ValueError(
                "index must be less than the number of lambda1 values (len(lambda1))"
            )
        if not isinstance(index, int):
            raise TypeError("index must be an integer")
        self._current_index = index

    def _get_window_index(self):
        # Internal function to get the current window index
        try:
            return self._current_index
        except:
            return None
