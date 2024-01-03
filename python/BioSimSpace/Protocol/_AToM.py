from .. import Types as _Types
from ._protocol import Protocol as _Protocol
from ._position_restraint_mixin import _PositionRestraintMixin
from .. import Units as _Units

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
        directions=[
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
            -1,
        ],
        lambda1=[
            0.00,
            0.05,
            0.10,
            0.15,
            0.20,
            0.25,
            0.30,
            0.35,
            0.40,
            0.45,
            0.50,
            0.50,
            0.45,
            0.40,
            0.35,
            0.30,
            0.25,
            0.20,
            0.15,
            0.10,
            0.05,
            0.00,
        ],
        lambda2=[
            0.00,
            0.05,
            0.10,
            0.15,
            0.20,
            0.25,
            0.30,
            0.35,
            0.40,
            0.45,
            0.50,
            0.50,
            0.45,
            0.40,
            0.35,
            0.30,
            0.25,
            0.20,
            0.15,
            0.10,
            0.05,
            0.00,
        ],
        alpha=[
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
        ],
        U0=[
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
        ],
        w0coeff=[
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
            0.00,
        ],
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

        w0coeff : list of float
            The w0coeff values.

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
        # TODO ADD checks for list length consistency
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

        # Store the w0coeff values.
        self.setW0coeff(w0coeff)

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

    def getDirections(self):
        """
        Return the direction values.

        Returns
        -------

        lambdas : list of float
            The directions.
        """
        return self._direction

    def setDirections(self, directions):
        """
        Set the direction values.

        Parameters
        ----------

        directions : list of int
            The directions.
        """
        if isinstance(directions, list):
            if all(item == 1 or item == -1 for item in directions):
                self._directions = directions
            else:
                raise ValueError("all entries in 'directions' must be either 1 or -1")
        else:
            raise TypeError("'directions' must be of type 'list'")

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
            if all(isinstance(item, float) for item in lambda1):
                self._lambda1 = lambda1
            else:
                raise ValueError("all entries in 'lambda1' must be floats")
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
            if all(isinstance(item, float) for item in lambda2):
                if len(lambda2) != len(self._lambda1):
                    raise ValueError(
                        "'lambda2' and 'lambda1' must have the same length"
                    )
                self._lambda2 = lambda2
            else:
                raise ValueError("all entries in 'lambda2' must be floats")
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
            if all(isinstance(item, float) for item in alpha):
                self._alpha = alpha
            else:
                raise ValueError("all entries in 'alpha' must be floats")
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
            if all(isinstance(item, float) for item in U0):
                self._U0 = U0
            else:
                raise ValueError("all entries in 'U0' must be floats")
        else:
            raise TypeError("'U0' must be of type 'list'")

    def getW0coeff(self):
        """
        Return the w0coeff values.

        Returns
        -------

        w0coeff : list of float
            The w0coeff values.
        """
        return self._w0coeff

    def setW0coeff(self, w0coeff):
        """
        Set the w0coeff values.

        Parameters
        ----------

        w0coeff : list of float
            The w0coeff values.
        """
        if isinstance(w0coeff, list):
            if all(isinstance(item, float) for item in w0coeff):
                self._w0coeff = w0coeff
            else:
                raise ValueError("all entries in 'w0coeff' must be floats")
        else:
            raise TypeError("'w0coeff' must be of type 'list'")

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

    def _set_current_lambdas(self, lam1, lam2):
        # Internal function to set the current lambda values for a specific process
        if not (isinstance(lam1, float) and lam1 >= 0.0 and lam1 <= 1.0):
            raise TypeError(
                "'lam1' must be of type 'float' with value between 0.0 and 1.0"
            )
        if not (isinstance(lam2, float) and lam2 >= 0.0 and lam2 <= 1.0):
            raise TypeError(
                "'lam2' must be of type 'float' with value between 0.0 and 1.0"
            )
        self._current_lambda1 = lam1
        self._current_lambda2 = lam2

    def _get_current_lambdas(self):
        # Internal function to get the current lambda values for a specific process
        return self._current_lambda1, self._current_lambda2
