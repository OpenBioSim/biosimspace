import math as _math
import warnings as _warnings

from sire.legacy import Units as _SireUnits

from .. import Protocol as _Protocol, _gmx_version
from ..Align._alch_ion import _get_protein_com_idx
from ..Align._squash import _amber_mask_from_indices, _squashed_atom_mapping
from ..FreeEnergy._restraint import Restraint as _Restraint
from ..Units.Energy import kj_per_mol as _kj_per_mol
from ..Units.Length import nanometer as _nanometer
from .._Exceptions import IncompatibleError as _IncompatibleError


class ConfigFactory:
    # TODO: Integrate this class better into the other Protocols.
    """A class for generating a config based on a template protocol."""

    def __init__(self, system, protocol, explicit_dummies=False):
        """
        Constructor.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.

        protocol : :class:`Protocol <BioSimSpace.Protocol>`
            The input protocol.

        explicit_dummies : bool
            Whether to keep the dummy atoms explicit at the endstates or remove them.
            This option is only used for AMBER.
        """
        self.system = system
        self.protocol = protocol
        self.explicit_dummies = explicit_dummies

    @property
    def _has_box(self):
        """Return whether the current system has a box."""
        if "space" in self.system._sire_object.propertyKeys():
            try:
                # Make sure that we have a periodic box. The system will now have
                # a default cartesian space.
                box = self.system._sire_object.property("space")
                has_box = box.isPeriodic()
            except:
                has_box = False
        else:
            _warnings.warn("No simulation box found. Assuming gas phase simulation.")
            has_box = False
        return has_box

    @property
    def _has_water(self):
        """Return whether the current system has any water molecules."""
        return self.system.nWaterMolecules() > 0

    @property
    def _report_interval(self):
        """Return the report interval based on the protocol value."""
        if isinstance(self.protocol, _Protocol.Minimisation):
            report_interval = 100
        else:
            report_interval = self.protocol.getReportInterval()
            if self._steps == 0:
                # Deal with the case where we just want single point energy.
                report_interval = 1
            elif report_interval > self._steps:
                report_interval = self._steps
        return report_interval

    @property
    def _restart(self):
        """Return whether this is a restart simulation."""
        try:
            return self.protocol.isRestart()
        except:
            pass
        return False

    @property
    def _restart_interval(self):
        """Return the restart interval based on the protocol value."""
        if isinstance(self.protocol, _Protocol.Minimisation):
            restart_interval = None
        else:
            restart_interval = self.protocol.getRestartInterval()
            if restart_interval > self._steps:
                restart_interval = self._steps
        return restart_interval

    @property
    def _steps(self):
        # Return the number of steps based on the protocol value.
        if isinstance(self.protocol, _Protocol.Minimisation):
            steps = self.protocol.getSteps()
        else:
            steps = _math.ceil(self.protocol.getRunTime() / self.protocol.getTimeStep())
        return steps

    def _generate_amber_fep_masks(self, timestep):
        """Internal helper function which generates timasks and scmasks based on the system.

        Parameters
        ----------

        timestep : [float]
            The timestep in ps for the FEP perturbation. Generates a different mask based on this.

        Returns
        -------

        option_dict : dict
            A dictionary of AMBER-compatible options.
        """
        # Get the merged to squashed atom mapping of the whole system for both endpoints.
        kwargs = dict(environment=False, explicit_dummies=self.explicit_dummies)
        mcs_mapping0 = _squashed_atom_mapping(
            self.system, is_lambda1=False, common=True, dummies=False, **kwargs
        )
        mcs_mapping1 = _squashed_atom_mapping(
            self.system, is_lambda1=True, common=True, dummies=False, **kwargs
        )
        dummy_mapping0 = _squashed_atom_mapping(
            self.system, is_lambda1=False, common=False, dummies=True, **kwargs
        )
        dummy_mapping1 = _squashed_atom_mapping(
            self.system, is_lambda1=True, common=False, dummies=True, **kwargs
        )

        # Generate the TI and dummy masks.
        mcs0_indices, mcs1_indices, dummy0_indices, dummy1_indices = [], [], [], []
        for i in range(self.system.nAtoms()):
            if i in dummy_mapping0:
                dummy0_indices.append(dummy_mapping0[i])
            if i in dummy_mapping1:
                dummy1_indices.append(dummy_mapping1[i])
            if i in mcs_mapping0:
                mcs0_indices.append(mcs_mapping0[i])
            if i in mcs_mapping1:
                mcs1_indices.append(mcs_mapping1[i])
        ti0_indices = mcs0_indices + dummy0_indices
        ti1_indices = mcs1_indices + dummy1_indices

        # SHAKE should be used for timestep > 1 fs.
        if timestep >= 0.002:
            no_shake_mask = ""
        else:
            no_shake_mask = _amber_mask_from_indices(ti0_indices + ti1_indices)

        # Create an option dict with amber masks generated from the above indices.
        option_dict = {
            "timask1": f'"{_amber_mask_from_indices(ti0_indices)}"',
            "timask2": f'"{_amber_mask_from_indices(ti1_indices)}"',
            "scmask1": f'"{_amber_mask_from_indices(dummy0_indices)}"',
            "scmask2": f'"{_amber_mask_from_indices(dummy1_indices)}"',
            "noshakemask": f'"{no_shake_mask}"',
        }

        return option_dict

    def generateAmberConfig(self, extra_options=None, extra_lines=None):
        """
        Outputs the current protocol in a format compatible with AMBER.

        Parameters
        ----------

        extra_options : dict
            A dictionary containing extra options. Overrides the ones generated from the protocol.

        extra_lines : list
            A list of extra lines to be put at the end of the script.

        Returns
        -------

        config : list
            The generated config list in an AMBER format.
        """

        extra_options = extra_options if extra_options is not None else {}
        extra_lines = extra_lines if extra_lines is not None else []

        protocol_lines = []

        # Define some miscellaneous defaults.
        protocol_dict = {
            "ntpr": 200,  # Interval between reporting energies.
            "ntwr": self._restart_interval,  # Interval between saving restart files.
            "ntwx": self._restart_interval,  # Trajectory sampling frequency.
            "ntxo": 2,  # Output coordinates as NetCDF.
            "irest": int(self._restart),  # Whether to restart.
        }

        # Input.
        if self._restart:
            protocol_dict["ntx"] = 5  # Read coordinates and velocities.
        else:
            protocol_dict["ntx"] = 1  # Only read coordinates from file.

        # Minimisation.
        if isinstance(self.protocol, _Protocol.Minimisation):
            # Work out the number of steepest descent cycles.
            # This is 1000 or 10% of the number of steps, whichever is larger.
            if self._steps <= 1000:
                num_steep = self._steps
            else:
                num_steep = _math.ceil(self._steps / 10)
                if num_steep < 1000:
                    num_steep = 1000

            protocol_dict["imin"] = 1  # Minimisation simulation.
            protocol_dict["ntmin"] = 2  # Set the minimisation method to XMIN
            protocol_dict["maxcyc"] = self._steps  # Set the number of steps.
            protocol_dict["ncyc"] = (
                num_steep  # Set the number of steepest descent steps.
            )
            # FIX need to remove and fix this, only for initial testing
            timestep = 0.004
        else:
            # Define the timestep
            timestep = (
                self.protocol.getTimeStep().picoseconds().value()
            )  # Get the time step in ps
            protocol_dict["dt"] = f"{timestep:.3f}"  # Time step.
            protocol_dict["nstlim"] = self._steps  # Number of integration steps.

        # Constraints.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            protocol_dict["ntc"] = 2  # Enable SHAKE.
            protocol_dict["ntf"] = 2  # Don't calculate forces for constrained bonds.

        # PBC.
        if self._has_box:
            protocol_dict["cut"] = "8.0"  # Non-bonded cut-off.
            protocol_dict["iwrap"] = 1  # Wrap the coordinates.
        else:
            protocol_dict["ntb"] = 0  # No periodic box.
            protocol_dict["cut"] = "999."  # Non-bonded cut-off.

        # Restraints.
        if isinstance(self.protocol, _Protocol._PositionRestraintMixin):
            # Restrain the backbone.
            restraint = self.protocol.getRestraint()

            # Convert to a squashed representation, if needed
            if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
                atom_mapping0 = _squashed_atom_mapping(self.system, is_lambda1=False)
                atom_mapping1 = _squashed_atom_mapping(self.system, is_lambda1=True)

            if self.system.getAlchemicalIon():
                alchem_ion_idx = self.system.getAlchemicalIonIdx()
                protein_com_idx = _get_protein_com_idx(self.system)
                atom_idxs = [alchem_ion_idx, protein_com_idx]
                if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
                    atom_idxs = sorted(
                        {atom_mapping0[x] for x in atom_idxs if x in atom_mapping0}
                        | {atom_mapping1[x] for x in atom_idxs if x in atom_mapping1}
                    )
                # Convert to 1-based index
                alchemical_ion_mask = "@" + ",".join(
                    [str(atom_idx + 1) for atom_idx in atom_idxs]
                )
            else:
                alchemical_ion_mask = None

            if restraint is not None:
                # Get the indices of the atoms that are restrained.
                if type(restraint) is str:
                    atom_idxs = self.system.getRestraintAtoms(restraint)
                else:
                    atom_idxs = restraint

                # Convert to a squashed representation, if needed
                if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
                    atom_idxs = sorted(
                        {atom_mapping0[x] for x in atom_idxs if x in atom_mapping0}
                        | {atom_mapping1[x] for x in atom_idxs if x in atom_mapping1}
                    )

                # Don't add restraints if there are no atoms to restrain.
                if len(atom_idxs) > 0:
                    # Generate the restraint mask based on atom indices.
                    restraint_mask = _amber_mask_from_indices(atom_idxs)

                    # The restraintmask cannot be more than 256 characters.
                    if len(restraint_mask) > 256:
                        # AMBER has a limit on the length of the restraintmask
                        # so it's easy to overflow if we are matching by index
                        # on a large protein. As such, handle "backbone" and
                        # "heavy" restraints using a non-interoperable name mask.
                        if type(restraint) is str:
                            if restraint == "backbone":
                                # Determine wether the system contains protein, nucleic acid, or both.
                                restraint_atom_names = []
                                if self.system.nAminoAcids() > 0:
                                    restraint_atom_names += ["N", "CA", "C", "O"]
                                if self.system.nNucleotides() > 0:
                                    restraint_atom_names += [
                                        "P",
                                        "C5'",
                                        "C3'",
                                        "O3'",
                                        "O5'",
                                    ]
                                restraint_mask = "@" + ",".join(restraint_atom_names)
                            elif restraint == "heavy":
                                restraint_mask = "!:WAT & !@%NA,CL & !@H="
                            elif restraint == "all":
                                restraint_mask = "!:WAT & !@%NA,CL"

                        # We can't do anything about a custom restraint, since we don't
                        # know anything about the atoms.
                        else:
                            raise ValueError(
                                "AMBER atom 'restraintmask' exceeds 256 character limit!"
                            )

            else:
                restraint_mask = None

            if restraint_mask or alchemical_ion_mask:
                if restraint_mask and alchemical_ion_mask:
                    restraint_mask = f"{restraint_mask} | {alchemical_ion_mask}"
                elif alchemical_ion_mask:
                    restraint_mask = alchemical_ion_mask

                protocol_dict["ntr"] = 1
                force_constant = self.protocol.getForceConstant()._sire_unit
                force_constant = force_constant.to(
                    _SireUnits.kcal_per_mol / _SireUnits.angstrom2
                )
                protocol_dict["restraint_wt"] = force_constant
                protocol_dict["restraintmask"] = f'"{restraint_mask}"'

        # Pressure control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            if self.protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if self._has_box and self._has_water:
                    protocol_dict["ntp"] = 1  # Isotropic pressure scaling.
                    protocol_dict["pres0"] = (
                        f"{self.protocol.getPressure().bar().value():.5f}"  # Pressure in bar.
                    )
                    if isinstance(self.protocol, _Protocol.Equilibration):
                        protocol_dict["barostat"] = 1  # Berendsen barostat.
                    else:
                        protocol_dict["barostat"] = 2  # Monte Carlo barostat.
                else:
                    _warnings.warn(
                        "Cannot use a barostat for a vacuum or non-periodic simulation"
                    )
            else:
                protocol_dict["ntb"] = 1  # constant volume.

        # Temperature control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            protocol_dict["ntt"] = 3  # Langevin dynamics.
            protocol_dict["gamma_ln"] = "{:.5f}".format(
                1 / self.protocol.getTauT().picoseconds().value()
            )  # Collision frequency (ps^-1).
            if isinstance(self.protocol, _Protocol.Equilibration):
                temp0 = self.protocol.getStartTemperature().kelvin().value()
                temp1 = self.protocol.getEndTemperature().kelvin().value()
                if not self.protocol.isConstantTemp():
                    protocol_dict["tempi"] = f"{temp0:.2f}"  # Initial temperature.
                    protocol_dict["temp0"] = f"{temp1:.2f}"  # Final temperature.
                    protocol_dict["nmropt"] = 1
                    protocol_lines += [
                        f"&wt TYPE='TEMP0', istep1=0, istep2={self._steps}, value1={temp0:.2f}, value2={temp1:.2f} /"
                    ]
                else:
                    if not self._restart:
                        protocol_dict["tempi"] = f"{temp0:.2f}"  # Initial temperature.
                    protocol_dict["temp0"] = f"{temp0:.2f}"  # Constant temperature.
            else:
                temp = self.protocol.getTemperature().kelvin().value()
                if not self._restart:
                    protocol_dict["tempi"] = f"{temp:.2f}"  # Initial temperature.
                protocol_dict["temp0"] = f"{temp:.2f}"  # Final temperature.

        # Free energies.
        if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
            protocol_dict["icfe"] = 1  # Free energy mode.
            protocol_dict["ifsc"] = 1  # Use softcore potentials.
            protocol_dict["ntf"] = 1  # Remove SHAKE constraints.
            lambda_values = self.protocol.getLambdaValues(type="dataframe")
            protocol = [f"{lam:.5f}" for lam in lambda_values["fep"]]
            protocol_dict["mbar_states"] = len(protocol)  # Number of lambda values.
            protocol_dict["mbar_lambda"] = ", ".join(protocol)  # Lambda values.
            Lambda = self.protocol.getLambda(type="series")
            protocol_dict["clambda"] = "{:.5f}".format(
                Lambda["fep"]
            )  # Current lambda value.

            if isinstance(self.protocol, _Protocol.Production):
                protocol_dict["ifmbar"] = 1  # Calculate MBAR energies.
                protocol_dict["logdvdl"] = 1  # Output dVdl
            protocol_dict = {
                **protocol_dict,
                **self._generate_amber_fep_masks(timestep),
            }  # Atom masks.

        # Put everything together in a line-by-line format.
        total_dict = {**protocol_dict, **extra_options}
        dict_lines = [self.protocol.__class__.__name__, "&cntrl"]
        dict_lines += [
            f"   {k}={v}," for k, v in total_dict.items() if v is not None
        ] + ["/"]
        total_lines = protocol_lines + extra_lines
        if total_lines:
            total_lines += ["&wt TYPE='END' /"]
        total_lines = dict_lines + total_lines

        return total_lines

    def generateGromacsConfig(
        self,
        extra_options=None,
        extra_lines=None,
        restraint=None,
        perturbation_type=None,
    ):
        """
        Outputs the current protocol in a format compatible with GROMACS.

        Parameters
        ----------

        extra_options : dict
            A dictionary containing extra options. Overrides the ones generated from the protocol.

        extra_lines : list
            A list of extra lines to be put at the end of the script.

        restraint : :class:`Restraint <BioSimSpace.FreeEnergy.Restraint>`
            Restraint object that contains information for ABFE calculations.

        perturbation_type : str
            The type of perturbation to perform. Options are:
            "full" : A full perturbation of all terms (default option).
            "release_restraint" : Used with multiple distance restraints to release all
                                  restraints other than the "permanent" one when the ligand
                                  is fully decoupled. Note that lambda = 0.0 is the fully
                                  released state, and lambda = 1.0 is the fully restrained
                                  state (i.e. 0.0 -> value). The non-permanent restraints
                                  can be scaled with bonded-lambda.

        Returns
        -------

        config : list
            The generated config list in a GROMACS format.
        """
        if perturbation_type == "release_restraint" and restraint is None:
            raise ValueError(
                "Cannot use perturbation_type='release_restraint' without a Restraint object."
            )

        extra_options = extra_options if extra_options is not None else {}
        extra_lines = extra_lines if extra_lines is not None else []

        # Define some miscellaneous defaults.
        protocol_dict = {
            "nstlog": self._report_interval,  # Interval between writing to the log file.
            "nstenergy": self._report_interval,  # Interval between writing to the energy file.
            "nstxout-compressed": self._restart_interval,  # Interval between writing to the trajectory file.
        }

        # Minimisation.
        if isinstance(self.protocol, _Protocol.Minimisation):
            protocol_dict["integrator"] = "steep"  # Minimisation simulation.
        else:
            timestep = (
                self.protocol.getTimeStep().picoseconds().value()
            )  # Define the timestep in picoseconds
            protocol_dict["dt"] = f"{timestep:.4f}"  # Integration time step.
        protocol_dict["nsteps"] = self._steps  # Number of integration steps.

        # Constraints.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            protocol_dict["constraints"] = "h-bonds"  # Rigid bonded hydrogens.
            protocol_dict["constraint-algorithm"] = "LINCS"  # Linear constraint solver.

        # PBC.
        protocol_dict["pbc"] = "xyz"  # Simulate a fully periodic box.
        protocol_dict["cutoff-scheme"] = "Verlet"  # Use Verlet pair lists.
        if self._has_box and self._has_water:
            protocol_dict["ns-type"] = "grid"  # Use a grid to search for neighbours.
            protocol_dict["nstlist"] = (
                "20"  # Rebuild neighbour list every 20 steps. Recommended in the manual for parallel simulations and/or non-bonded force calculation on the GPU.
            )
            protocol_dict["rlist"] = "0.8"  # Set short-range cutoff.
            protocol_dict["rvdw"] = "0.8"  # Set van der Waals cutoff.
            protocol_dict["rcoulomb"] = "0.8"  # Set Coulomb cutoff.
            protocol_dict["coulombtype"] = "PME"  # Fast smooth Particle-Mesh Ewald.
            protocol_dict["DispCorr"] = (
                "EnerPres"  # Dispersion corrections for energy and pressure.
            )
        else:
            # Perform vacuum simulations by implementing pseudo-PBC conditions,
            # i.e. run calculation in a near-infinite box (333.3 nm).
            # c.f.: https://pubmed.ncbi.nlm.nih.gov/29678588
            protocol_dict["nstlist"] = (
                "1"  # Single neighbour list (all particles interact).
            )
            protocol_dict["rlist"] = "333.3"  # "Infinite" short-range cutoff.
            protocol_dict["rvdw"] = "333.3"  # "Infinite" van der Waals cutoff.
            protocol_dict["rcoulomb"] = "333.3"  # "Infinite" Coulomb cutoff.
            protocol_dict["coulombtype"] = "Cut-off"  # Plain cut-off.
        protocol_dict["vdwtype"] = "Cut-off"  # Twin-range van der Waals cut-off.

        # Pressure control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            if self.protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if self._has_box and self._has_water:
                    if _gmx_version >= 2021:
                        protocol_dict["pcoupl"] = "C-rescale"  # C-rescale barostat.
                    else:
                        protocol_dict["pcoupl"] = "Berendsen"  # Berendsen barostat.
                    # Do the MC move every 100 steps to be the same as AMBER.
                    protocol_dict["nstpcouple"] = 100
                    # 4ps time constant for pressure coupling.
                    # As the tau-p has to be 10 times larger than nstpcouple * dt (4 fs)
                    protocol_dict["tau-p"] = 4
                    protocol_dict["ref-p"] = (
                        f"{self.protocol.getPressure().bar().value():.5f}"  # Pressure in bar.
                    )
                    protocol_dict["compressibility"] = (
                        "4.5e-5"  # Compressibility of water.
                    )
                else:
                    _warnings.warn(
                        "Cannot use a barostat for a vacuum or non-periodic simulation"
                    )

        # Temperature control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
                protocol_dict["integrator"] = "sd"  # langevin dynamics.
            else:
                protocol_dict["integrator"] = "md"  # leap-frog dynamics.
                protocol_dict["tcoupl"] = "v-rescale"
            protocol_dict["tc-grps"] = (
                "system"  # A single temperature group for the entire system.
            )
            protocol_dict["tau-t"] = "{:.5f}".format(
                self.protocol.getTauT().picoseconds().value()
            )  # Collision frequency (ps).

            if isinstance(self.protocol, _Protocol.Equilibration):
                if self.protocol.isConstantTemp():
                    temp = "%.2f" % self.protocol.getStartTemperature().kelvin().value()
                else:
                    # still need a reference temperature for each group, even when heating/cooling
                    temp = "%.2f" % self.protocol.getEndTemperature().kelvin().value()
                    # Work out the final time of the simulation.
                    timestep = self.protocol.getTimeStep().picoseconds().value()
                    end_time = _math.floor(timestep * self._steps)

                    protocol_dict["annealing"] = (
                        "single"  # Single sequence of annealing points.
                    )
                    protocol_dict["annealing-npoints"] = (
                        2  # Two annealing points for "system" temperature group.
                    )

                    # Linearly change temperature between start and end times.
                    protocol_dict["annealing-time"] = "0 %d" % end_time
                    protocol_dict["annealing-temp"] = "%.2f %.2f" % (
                        self.protocol.getStartTemperature().kelvin().value(),
                        self.protocol.getEndTemperature().kelvin().value(),
                    )
            else:
                temp = "%.2f" % self.protocol.getTemperature().kelvin().value()
            protocol_dict["ref-t"] = temp

            # Regenerate the velocity if this is not a restart
            if not self.protocol.isRestart():
                protocol_dict["gen-vel"] = "yes"
                protocol_dict["gen-temp"] = temp

        # Free energies.
        if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
            protocol_dict["free-energy"] = "yes"  # Free energy mode.
            nDecoupledMolecules = self.system.nDecoupledMolecules()
            if nDecoupledMolecules == 1:
                [
                    mol,
                ] = self.system.getDecoupledMolecules()
                if not mol.isPerturbable():
                    decouple_dict = mol._sire_object.property("decouple")
                    protocol_dict["couple-moltype"] = mol._sire_object.name().value()

                    def tranform(charge, LJ):
                        if charge and LJ:
                            return "vdw-q"
                        elif charge and not LJ:
                            return "q"
                        elif not charge and LJ:
                            return "vdw"
                        else:
                            return "none"

                    protocol_dict["couple-lambda0"] = tranform(
                        decouple_dict["charge"][0], decouple_dict["LJ"][0]
                    )
                    protocol_dict["couple-lambda1"] = tranform(
                        decouple_dict["charge"][1], decouple_dict["LJ"][1]
                    )
                    if decouple_dict["intramol"].value():
                        # The intramol is being coupled to the lambda change and thus being annihilated.
                        protocol_dict["couple-intramol"] = "yes"
                    else:
                        protocol_dict["couple-intramol"] = "no"
                # Add the soft-core parameters for the ABFE
                protocol_dict["sc-alpha"] = 0.5
                protocol_dict["sc-power"] = 1
                protocol_dict["sc-sigma"] = 0.3
            elif nDecoupledMolecules > 1:
                raise ValueError(
                    "Gromacs cannot handle more than one decoupled molecule."
                )
            protocol_dict["calc-lambda-neighbors"] = -1  # Calculate MBAR energies.
            LambdaValues = self.protocol.getLambdaValues(type="dataframe")
            for name in [
                "fep",
                "bonded",
                "coul",
                "vdw",
                "restraint",
                "mass",
                "temperature",
            ]:
                if name in LambdaValues:
                    protocol_dict["{:<20}".format("{}-lambdas".format(name))] = (
                        " ".join(
                            list(map("{:.5f}".format, LambdaValues[name].to_list()))
                        )
                    )
            protocol_dict["init-lambda-state"] = (
                self.protocol.getLambdaIndex()
            )  # Current lambda value.
            protocol_dict["nstcalcenergy"] = (
                self._report_interval
            )  # Calculate energies every report_interval steps.
            protocol_dict["nstdhdl"] = (
                self._report_interval
            )  # Write gradients every report_interval steps.

            # Handle the combination of multiple distance restraints and perturbation type
            # of "release_restraint". In this case, the force constant of the "permanent"
            # restraint needs to be written to the mdp file.
            if perturbation_type == "release_restraint":
                if not isinstance(restraint, _Restraint):
                    raise ValueError(
                        "Cannot use perturbation_type='release_restraint' without a Restraint object."
                    )
                # Get force constant in correct units.
                force_constant = restraint._restraint_dict[
                    "permanent_distance_restraint"
                ]["kr"] / (_kj_per_mol / _nanometer**2)
                protocol_dict["disre-fc"] = force_constant

        # Put everything together in a line-by-line format.
        total_dict = {**protocol_dict, **extra_options}
        total_lines = [
            f"{k} = {v}" for k, v in total_dict.items() if v is not None
        ] + extra_lines

        return total_lines

    def generateSomdConfig(
        self,
        extra_options=None,
        extra_lines=None,
        restraint=None,
        perturbation_type=None,
    ):
        """Outputs the current protocol in a format compatible with SOMD.

        Parameters
        ----------

        extra_options : dict
            A dictionary containing extra options. Overrides the ones generated from the protocol.

        extra_lines : list
            A list of extra lines to be put at the end of the script.

        restraint : :class:`Restraint <BioSimSpace.FreeEnergy.Restraint>`
            The Restraint object that contains information for the ABFE
            calculations.

        perturbation_type : str
            The type of perturbation to perform. Options are:
            "full" : A full perturbation of all terms (default option).
            "discharge_soft" : Perturb all discharging soft atom charge terms (i.e. value->0.0).
            "vanish_soft" : Perturb all vanishing soft atom LJ terms (i.e. value->0.0).
            "flip" : Perturb all hard atom terms as well as bonds/angles.
            "grow_soft" : Perturb all growing soft atom LJ terms (i.e. 0.0->value).
            "charge_soft" : Perturb all charging soft atom LJ terms (i.e. 0.0->value).
            "restraint" : Perturb the receptor-ligand restraint strength by linearly
                          scaling the force constants (0.0->value).
            "release_restraint" : Used with multiple distance restraints to release all
                                  restraints other than the "permanent" one when the ligand
                                  is fully decoupled. Note that lambda = 0.0 is the fully
                                  released state, and lambda = 1.0 is the fully restrained
                                  state (i.e. 0.0 -> value).

        Returns
        -------

        config : list
            The generated config list in a SOMD format.
        """

        extra_options = extra_options if extra_options is not None else {}
        extra_lines = extra_lines if extra_lines is not None else []

        # Define some miscellaneous defaults.
        # Save molecular coordinates.
        protocol_dict = {"save coordinates": True}

        # Minimisation.
        if isinstance(self.protocol, _Protocol.Minimisation):
            # Minimisation simulation.
            protocol_dict["minimise"] = True
            # Maximum number of steps.
            protocol_dict["minimise maximum iterations"] = self._steps
            # Convergence tolerance.
            protocol_dict["minimise tolerance"] = 1
            # Perform a single SOMD cycle.
            protocol_dict["ncycles"] = 1
            # Perform a single MD move.
            protocol_dict["nmoves"] = 1
        else:
            # Get the report and restart intervals.
            report_interval = self._report_interval
            restart_interval = self._restart_interval
            runtime = self.protocol.getRunTime()

            # Get the timestep.
            timestep = self.protocol._timestep

            # Work out the number of cycles.
            ncycles = (runtime / timestep) / report_interval

            # If the number of cycles isn't integer valued, adjust the report
            # interval so that we match specified the run time.
            if ncycles - _math.floor(ncycles) != 0:
                ncycles = _math.floor(ncycles)
                if ncycles == 0:
                    ncycles = 1
                report_interval = _math.ceil((runtime / timestep) / ncycles)

            # For free energy simulations, the report interval must be a multiple
            # of the energy frequency which is 250 steps.
            if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
                if report_interval % 250 != 0:
                    report_interval = 250 * _math.ceil(report_interval / 250)

            # Work out the number of cycles per frame.
            cycles_per_frame = restart_interval / report_interval

            # Work out whether we need to adjust the buffer frequency.
            buffer_freq = 0
            if cycles_per_frame < 1:
                buffer_freq = cycles_per_frame * restart_interval
                cycles_per_frame = 1
                self._buffer_freq = buffer_freq
            else:
                cycles_per_frame = _math.floor(cycles_per_frame)

            # For free energy simulations, the buffer frequency must be an integer
            # multiple of the frequency at which free energies are written, which
            # is 250 steps. Round down to the closest multiple.
            if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
                if buffer_freq > 0:
                    buffer_freq = 250 * _math.floor(buffer_freq / 250)

            # The number of SOMD cycles.
            protocol_dict["ncycles"] = int(ncycles)
            # The number of moves per cycle.
            protocol_dict["nmoves"] = report_interval
            # Cycles per trajectory write.
            protocol_dict["ncycles_per_snap"] = cycles_per_frame
            # Buffering frequency.
            protocol_dict["buffered coordinates frequency"] = buffer_freq
            timestep = self.protocol.getTimeStep().femtoseconds().value()
            # Integration time step.
            protocol_dict["timestep"] = "%.2f femtosecond" % timestep

            # Use the Langevin Middle integrator if it is a 4 fs timestep
            if timestep >= 4.00:
                # Langevin middle integrator
                protocol_dict["integrator_type"] = "langevinmiddle"
            else:
                pass

        # PBC.
        if self._has_water:
            protocol_dict["reaction field dielectric"] = "78.3"  # Solvated box.
        if not self._has_box or not self._has_water:
            protocol_dict["cutoff type"] = "cutoffnonperiodic"  # No periodic box.
        else:
            protocol_dict["cutoff type"] = "cutoffperiodic"  # Periodic box.
        protocol_dict["cutoff distance"] = "10 angstrom"  # Non-bonded cut-off.

        # Restraints.
        if (
            isinstance(self.protocol, _Protocol._PositionRestraintMixin)
            and self.protocol.getRestraint() is not None
        ):
            raise _IncompatibleError("We currently don't support restraints with SOMD.")

        # Pressure control.
        protocol_dict["barostat"] = False  # Disable barostat (constant volume).
        if not isinstance(self.protocol, _Protocol.Minimisation):
            if self.protocol.getPressure() is not None:
                # Don't use barostat for vacuum simulations.
                if self._has_box and self._has_water:
                    protocol_dict["barostat"] = True  # Enable barostat.
                    pressure = self.protocol.getPressure().atm().value()
                    protocol_dict["pressure"] = (
                        "%.5f atm" % pressure
                    )  # Presure in atmosphere.
                else:
                    _warnings.warn(
                        "Cannot use a barostat for a vacuum or non-periodic simulation"
                    )

        # Temperature control.
        if not isinstance(self.protocol, _Protocol.Minimisation):
            if (
                isinstance(self.protocol, _Protocol.Equilibration)
                and not self.protocol.isConstantTemp()
            ):
                raise _IncompatibleError(
                    "SOMD only supports constant temperature equilibration."
                )

            protocol_dict["thermostat"] = "True"  # Turn on the thermostat.
            if not isinstance(self.protocol, _Protocol.Equilibration):
                protocol_dict["temperature"] = (
                    "%.2f kelvin" % self.protocol.getTemperature().kelvin().value()
                )
            else:
                protocol_dict["temperature"] = (
                    "%.2f kelvin" % self.protocol.getStartTemperature().kelvin().value()
                )

        # Free energies.
        if isinstance(self.protocol, _Protocol._FreeEnergyMixin):
            if not isinstance(self.protocol, _Protocol.Minimisation):
                protocol_dict["constraint"] = (
                    "hbonds-notperturbed"  # Handle hydrogen perturbations.
                )
                protocol_dict["energy frequency"] = (
                    250  # Write gradients every 250 steps.
                )

            protocol = [str(x) for x in self.protocol.getLambdaValues()]
            protocol_dict["lambda array"] = ", ".join(protocol)
            protocol_dict["lambda_val"] = (
                self.protocol.getLambda()
            )  # Current lambda value.

            try:  # RBFE
                res_num = (
                    self.system.search("perturbable")
                    .residues()[0]
                    ._sire_object.number()
                    .value()
                )
            except ValueError:  # No perturbable molecule - try ABFE
                res_num = (
                    self.system.getDecoupledMolecules()[0]
                    .getResidues()[0]
                    ._sire_object.number()
                    .value()
                )

            protocol_dict["perturbed residue number"] = (
                res_num  # Perturbed residue number.
            )

        # Put everything together in a line-by-line format.
        total_dict = {**protocol_dict, **extra_options}
        total_lines = [
            f"{k} = {v}" for k, v in total_dict.items() if v is not None
        ] + extra_lines

        # Restraint
        if restraint:
            # Handle Boresch restraints.
            if restraint._restraint_type == "boresch":
                if perturbation_type == "release_restraint":
                    raise _IncompatibleError(
                        "SOMD does not support releasing Boresch restraints, only "
                        "multiple distance restraints."
                    )
                total_lines.append("use boresch restraints = True")
                total_lines.append(restraint.toString(engine="SOMD"))

            # Handle multiple distance restraints.
            elif restraint._restraint_type == "multiple_distance":
                total_lines.append("use distance restraints = True")
                if perturbation_type != "restraint":
                    # In this case, we want to ensure that the "permanent" distance restraint is active
                    total_lines.append("use permanent distance restraints = True")
                total_lines.append(
                    restraint.toString(
                        engine="SOMD", perturbation_type=perturbation_type
                    )
                )

            # If we are turning on the restraint, need to specify this in the config file.
            if (
                perturbation_type == "restraint"
                or perturbation_type == "release_restraint"
            ):
                total_lines.append("turn on receptor-ligand restraints mode = True")

        return total_lines
