######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2023
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

"""
A class for holding restraints.
"""
from math import e
from scipy import integrate as _integrate
import numpy as _np
import warnings as _warnings

from sire.legacy.Units import angstrom3 as _Sire_angstrom3
from sire.legacy.Units import k_boltz as _k_boltz
from sire.legacy.Units import meter3 as _Sire_meter3
from sire.legacy.Units import mole as _Sire_mole

from .._SireWrappers import Atom as _Atom
from .._SireWrappers import System as _System
from ..Types import Angle as _Angle
from ..Types import Length as _Length
from ..Types import Temperature as _Temperature
from ..Units.Angle import degree as _degree
from ..Units.Angle import radian as _radian
from ..Units.Area import angstrom2 as _angstrom2
from ..Units.Energy import kcal_per_mol as _kcal_per_mol
from ..Units.Energy import kj_per_mol as _kj_per_mol
from ..Units.Length import angstrom as _angstrom
from ..Units.Length import nanometer as _nanometer
from ..Units.Temperature import kelvin as _kelvin


class Restraint:
    """
    The Restraint class which holds the restraint information for the ABFE
    calculations. Currently Boresch restaraints and multiple distance restraints
    (between pairs of atoms) are supported. For the multiple distance restraints,
    it is assumed that all restraints but one will be released after decoupling,
    and an analytical correction will be applied to account for releasing the last
    restraint.

    Boresch restraint is a set of harmonic restraints containing one bond, two
    angle and three dihedrals, which comes from three atoms in the ligand
    (l1, l2, l3) and three atoms in the protein (r1, r2, r3). The restraints
    are represented in the format of:
    atom1-atom2-... (equilibrium value, force constant)

    The nomenclature:
    Bonds: r1-l1 (r0, kr)
    Angles: r2-r1-l1 (thetaA0, kthetaA), r1-l1-l2 (thetaB0, kthetaB)
    Dihedrals: r3-r2-r1-l1 (phiA0, kphiA), r2-r1-l1-l2 (phiB0, kphiB), r1-l1-l2-l3 (phiC0, kphiC)

    The Boresch restraint_dict has the following compact format.

    restraint_dict = {
        "anchor_points":{"r1": BioSimSpace._SireWrappers.Atom,
                         "r2": BioSimSpace._SireWrappers.Atom,
                         "r3": BioSimSpace._SireWrappers.Atom,
                         "l1": BioSimSpace._SireWrappers.Atom,
                         "l2": BioSimSpace._SireWrappers.Atom,
                         "l3": BioSimSpace._SireWrappers.Atom},
        "equilibrium_values":{"r0": BioSimSpace.Types.Length,
                              "thetaA0": BioSimSpace.Types.Angle,
                              "thetaB0": BioSimSpace.Types.Angle,
                              "phiA0": BioSimSpace.Types.Angle,
                              "phiB0": BioSimSpace.Types.Angle,
                              "phiC0": BioSimSpace.Types.Angle},
        "force_constants":{"kr": BioSimSpace.Types.Energy / BioSimSpace.Types.Area,
                           "kthetaA": BioSimSpace.Types.Energy / (BioSimSpace.Types.Area * BioSimSpace.Types.Area),
                           "kthetaB": BioSimSpace.Types.Energy / (BioSimSpace.Types.Area * BioSimSpace.Types.Area),
                           "kphiA": BioSimSpace.Types.Energy / (BioSimSpace.Types.Area * BioSimSpace.Types.Area),
                           "kphiB": BioSimSpace.Types.Energy / (BioSimSpace.Types.Area * BioSimSpace.Types.Area),
                           "kphiC": BioSimSpace.Types.Energy / (BioSimSpace.Types.Area * BioSimSpace.Types.Area)}}

    The multiple distance restraints are flat-bottom restraints, which are represented as a dictionaries
    of the form:

    distance_restraint_dict = {r1: BioSimSpace._SireWrappers.Atom,
                               l1: BioSimSpace._SireWrappers.Atom,
                               r0: BioSimSpace.Types.Length,
                               kr: BioSimSpace.Types.Energy / BioSimSpace.Types.Area,
                               r_fb: BioSimSpace.Types.Length}

    where r1 and l1 are the atoms involved in the restraint in the receptor and the ligand respectively,
    r0 and kr are the equilibrium distance and the force constant, and r_fb is the flat bottom radius.

    The overall restraint_dict is a dictionary of the form:

    restraint_dict = {"distance_restraints": [distance_restraint_dict, ...],
                      "permanent_distance_restraint": distance_restraint_dict}

    where distance_restraints is a list of distance restraints which are released after decoupling,
    and permanent_distance_restraint is the distance restraint which is not released after decoupling.
    """

    # Create a dict of supported restraints and compatible engines.
    supported_restraints = {
        "boresch": ["gromacs", "somd"],
        "multiple_distance": ["gromacs", "somd"],
    }

    def __init__(self, system, restraint_dict, temperature, restraint_type="Boresch"):
        """
        Constructor.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.

        restraint_dict : dict
            The dict for holding the restraint.

        temperature : :class:`System <BioSimSpace.Types.Temperature>`
            The temperature of the system

        restraint_type : str
            The type of the restraint. (`Boresch`, `multiple_distance`)
        """
        if not isinstance(temperature, _Temperature):
            raise ValueError(
                "temperature must be of type 'BioSimSpace.Types.Temperature'"
            )
        else:
            self.T = temperature

        if restraint_type.lower() == "boresch":
            self._restraint_type = "boresch"
            # Test if the atoms are of BioSimSpace._SireWrappers.Atom
            for key in ["r3", "r2", "r1", "l1", "l2", "l3"]:
                if not isinstance(restraint_dict["anchor_points"][key], _Atom):
                    raise ValueError(
                        f"restraint_dict['anchor_points']['{key}'] "
                        "must be of type "
                        "'BioSimSpace._SireWrappers.Atom'"
                    )

            # Test if the equilibrium length of the bond r1-l1 is a length unit
            # Such as angstrom or nanometer
            if not isinstance(restraint_dict["equilibrium_values"]["r0"], _Length):
                raise ValueError(
                    "restraint_dict['equilibrium_values']['r0'] must be of type 'BioSimSpace.Types.Length'"
                )

            # Test if the equilibrium length of the angle and dihedral is a
            # angle unit such as radian or degree
            for key in ["thetaA0", "thetaB0", "phiA0", "phiB0", "phiC0"]:
                if not isinstance(restraint_dict["equilibrium_values"][key], _Angle):
                    raise ValueError(
                        f"restraint_dict['equilibrium_values']['{key}'] must be "
                        f"of type 'BioSimSpace.Types.Angle'"
                    )

            # Test if the force constant of the angle and dihedral is the correct unit
            # Such as kcal/mol/rad^2
            for key in ["kthetaA", "kthetaB", "kphiA", "kphiB", "kphiC"]:
                if restraint_dict["force_constants"][key] != 0:
                    dim = restraint_dict["force_constants"][key].dimensions()
                    if dim != (-2, 0, 2, 1, -1, 0, -2):
                        raise ValueError(
                            f"restraint_dict['force_constants']['{key}'] must be of type "
                            f"'BioSimSpace.Types.Energy'/'BioSimSpace.Types.Angle^2'"
                        )

            # Test for unstable combinations of force constants
            non_zero_force_const = [
                i[0] for i in restraint_dict["force_constants"].items() if i[1] != 0
            ]
            if "kr" not in non_zero_force_const:
                raise ValueError('"kr" cannot be zero')
            if "kthetaA" not in non_zero_force_const:
                if "kphiA" in non_zero_force_const or "kphiB" in non_zero_force_const:
                    raise ValueError(
                        "Restraining phiA or phiB without restraining thetaA "
                        "will produce unstable Boresch restraints."
                    )
            if "kthetaB" not in non_zero_force_const:
                if "kphiB" in non_zero_force_const or "kphiC" in non_zero_force_const:
                    raise ValueError(
                        "Restraining phiB or phiC without restraining thetaB "
                        "will produce unstable Boresch restraints."
                    )

            # Test if the force constant of the bond r1-l1 is the correct unit
            # Such as kcal/mol/angstrom^2
            dim = restraint_dict["force_constants"]["kr"].dimensions()
            if dim != (0, 0, 0, 1, -1, 0, -2):
                raise ValueError(
                    "restraint_dict['force_constants']['kr'] must be of type "
                    "'BioSimSpace.Types.Energy'/'BioSimSpace.Types.Length^2'"
                )

            # Ensure restrained angles are >= 10 kT from collinear
            R = (
                _k_boltz.value() * _kcal_per_mol / _kelvin
            ).value()  # molar gas constant in kcal mol-1 K-1
            T = self.T / _kelvin  # Temperature in Kelvin

            for angle in ["thetaA", "thetaB"]:
                force_const = restraint_dict["force_constants"][f"k{angle}"] / (
                    _kcal_per_mol / (_radian * _radian)
                )
                if force_const != 0:
                    equil_val = (
                        restraint_dict["equilibrium_values"][f"{angle}0"] / _radian
                    )

                    # Convert 10 kT to angle
                    R = (
                        _k_boltz.value() * _kcal_per_mol / _kelvin
                    ).value()  # molar gas constant in kcal mol-1 K-1
                    T = self.T / _kelvin  # Temperature in Kelvin
                    min_stable_dist = _np.sqrt((20 * R * T) / force_const)
                    min_dist = min([abs(equil_val - 0), abs(equil_val - _np.pi)])

                    if min_dist < min_stable_dist:
                        _warnings.warn(
                            f"The equilibrium value of {angle} is within 10 kT of"
                            "collinearity, which may result in unstable Boresch restraints."
                            " Consider increasing the force constants or selecting equilibrium"
                            " values further from 0 or pi radians."
                        )

        elif restraint_type.lower() == "multiple_distance":
            self._restraint_type = "multiple_distance"

            if not set(restraint_dict.keys()) == {
                "distance_restraints",
                "permanent_distance_restraint",
            }:
                raise ValueError(
                    "restraint_dict must have keys 'distance_restraints' and 'permanent_distance_restraint'"
                )

            # Warn the user if there are no distance restraints (although they may be deliberately supplying
            # only the permanent distance restraint)
            if len(restraint_dict["distance_restraints"]) == 0:
                _warnings.warn(
                    "No distance restraints have been specified other than the permanent distance restraint."
                )

            # Check each distance restraint is of the correct format
            all_restraints = restraint_dict["distance_restraints"] + [
                restraint_dict["permanent_distance_restraint"]
            ]
            for single_restraint_dict in all_restraints:
                if not set(single_restraint_dict.keys()) == {
                    "r1",
                    "l1",
                    "r0",
                    "r_fb",
                    "kr",
                }:
                    raise ValueError(
                        "distance_restraint_dict must have keys 'r1', 'l1', 'r0', 'r_fb', 'kr' "
                        f"but has keys {list(single_restraint_dict.keys())}"
                    )

                # Check that the atoms are of type BioSimSpace._SireWrappers.Atom
                for key in ["r1", "l1"]:
                    if not isinstance(single_restraint_dict[key], _Atom):
                        raise ValueError(
                            f"distance_restraint_dict['{key}'] "
                            "must be of type 'BioSimSpace._SireWrappers.Atom'"
                        )

                # Test that all quantities have the correct units
                for key in ["r0", "r_fb"]:
                    if not isinstance(single_restraint_dict[key], _Length):
                        raise ValueError(
                            f"distance_restraint_dict['{key}'] must be of type "
                            "'BioSimSpace.Types.Length'"
                        )
                if not single_restraint_dict["kr"].dimensions() == (
                    0,
                    0,
                    0,
                    1,
                    -1,
                    0,
                    -2,
                ):
                    raise ValueError(
                        "distance_restraint_dict['kr'] must be of type "
                        "'BioSimSpace.Types.Energy'/'BioSimSpace.Types.Length^2'"
                    )

        else:
            raise NotImplementedError(
                f"Restraint type {type} not implemented "
                f"yet. Only {Restraint.supported_restraints.keys()} "
                "are currently supported."
            )

        self._restraint_dict = restraint_dict
        self.system = system

    @property
    def system(self):
        return self._system

    @system.setter
    def system(self, system):
        """
        Update the system object.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.
        """
        if not isinstance(system, _System):
            raise TypeError(
                "'system' must be of type 'BioSimSpace._SireWrappers.System'"
            )

        if self._restraint_type == "boresch":
            # Check if the ligand atoms are decoupled.
            # Find the decoupled molecule, assume that only one can be
            # decoupled.
            (decoupled_mol,) = system.getDecoupledMolecules()
            for key in ["l1", "l2", "l3"]:
                atom = self._restraint_dict["anchor_points"][key]
                # Discussed in https://github.com/michellab/BioSimSpace/pull/337
                if (
                    atom._sire_object.molecule().number()
                    != decoupled_mol._sire_object.number()
                ):
                    raise ValueError(
                        f"The ligand atom {key} is not from decoupled molecule."
                    )
            for key in ["r1", "r2", "r3"]:
                atom = self._restraint_dict["anchor_points"][key]
                if not atom in system:
                    raise ValueError(f"The receptor atom {key} is not in the system.")

        if self._restraint_type == "multiple_distance":
            # Check if the ligand atoms are decoupled.
            # Find the decoupled molecule, assume that only one can be
            # decoupled.
            (decoupled_mol,) = system.getDecoupledMolecules()

            all_restraints = self._restraint_dict["distance_restraints"] + [
                self._restraint_dict["permanent_distance_restraint"]
            ]
            for single_restraint_dict in all_restraints:
                ligand_atom = single_restraint_dict["l1"]
                if (
                    ligand_atom._sire_object.molecule().number()
                    != decoupled_mol._sire_object.number()
                ):
                    raise ValueError(
                        f"The ligand atom {ligand_atom} is not from decoupled molecule."
                    )
                receptor_atom = single_restraint_dict["r1"]
                if not receptor_atom in system:
                    raise ValueError(
                        f"The protein atom {receptor_atom} is not in the system."
                    )

        # Store a copy of solvated system.
        self._system = system.copy()

    def _gromacs_boresch(self, perturbation_type=None, restraint_lambda=False):
        """Format the Gromacs string for boresch restraint."""

        # Format the atoms into index list
        def format_index(key_list):
            formated_index = []
            for key in key_list:
                formated_index.append(
                    "{:<10}".format(
                        self._system.getIndex(
                            self._restraint_dict["anchor_points"][key]
                        )
                        + 1
                    )
                )
            return " ".join(formated_index)

        parameters_string = "{eq0:<10} {fc0:<10} {eq1:<10} {fc1:<10}"
        # The Gromacs dihedral restraints has the format of
        # phi dphi fc and we don't want dphi for the restraint, it is hence zero
        dihedral_restraints_parameters_string = (
            "{eq0:<10} 0.00 {fc0:<10} {eq1:<10} 0.00 {fc1:<10}"
        )

        # Format the parameters for the bonds
        def format_bond(equilibrium_values, force_constants):
            """
            Format the bonds equilibrium values and force constant
            in into the Gromacs topology format.
            """
            converted_equ_val = (
                self._restraint_dict["equilibrium_values"][equilibrium_values]
                / _nanometer
            )
            converted_fc = self._restraint_dict["force_constants"][force_constants] / (
                _kj_per_mol / _nanometer**2
            )
            return parameters_string.format(
                eq0="{:.3f}".format(converted_equ_val),
                fc0="{:.2f}".format(0),
                eq1="{:.3f}".format(converted_equ_val),
                fc1="{:.2f}".format(converted_fc),
            )

        # Format the parameters for the angles and dihedrals
        def format_angle(equilibrium_values, force_constants, restraint_lambda):
            """
            Format the angle equilibrium values and force constant
            in into the Gromacs topology format.

            For Boresch restraint, we might want the dihedral to be stored
            under the [ dihedral_restraints ] and controlled by restraint-lambdas.
            Instead of under the [ dihedrals ] directive and controlled by bonded-lambdas.

            However, for dihedrals, the [ dihedral_restraints ] has a different function type
            compared with [ dihedrals ] and more values for the force constant, so we need
            to format them differently.

            When restraint_lambda is True, the dihedrals will be stored in the dihedral_restraints.
            """
            converted_equ_val = (
                self._restraint_dict["equilibrium_values"][equilibrium_values] / _degree
            )
            converted_fc = self._restraint_dict["force_constants"][force_constants] / (
                _kj_per_mol / (_radian * _radian)
            )
            par_string = (
                dihedral_restraints_parameters_string
                if restraint_lambda
                else parameters_string
            )
            return par_string.format(
                eq0="{:.3f}".format(converted_equ_val),
                fc0="{:.2f}".format(0),
                eq1="{:.3f}".format(converted_equ_val),
                fc1="{:.2f}".format(converted_fc),
            )

        # basic format of the Gromacs string
        master_string = "  {index} {func_type} {parameters}"

        def write_bond(key_list, equilibrium_values, force_constants):
            return master_string.format(
                index=format_index(key_list),
                func_type=6,
                parameters=format_bond(equilibrium_values, force_constants),
            )

        def write_angle(key_list, equilibrium_values, force_constants):
            return master_string.format(
                index=format_index(key_list),
                func_type=1,
                parameters=format_angle(
                    equilibrium_values, force_constants, restraint_lambda=False
                ),
            )

        def write_dihedral(
            key_list, equilibrium_values, force_constants, restraint_lambda
        ):
            if restraint_lambda:
                # In [ dihedral_restraints ], function type 1
                # means the dihedral is restrained harmonically.
                func_type = 1
            else:
                # In [ dihedrals ], function type 2
                # means the dihedral is restrained harmonically.
                func_type = 2
            return master_string.format(
                index=format_index(key_list),
                func_type=func_type,
                parameters=format_angle(
                    equilibrium_values, force_constants, restraint_lambda
                ),
            )

        # Writing the string
        output = [
            "[ intermolecular_interactions ]",
        ]

        output.append("[ bonds ]")
        output.append("; ai         aj      type bA         kA         bB         kB")
        # Bonds: r1-l1 (r0, kr)
        output.append(write_bond(("r1", "l1"), "r0", "kr"))

        output.append("[ angles ]")
        output.append(
            "; ai         aj         ak      type thA        fcA        thB        fcB"
        )
        # Angles: r2-r1-l1 (thetaA0, kthetaA)
        output.append(write_angle(("r2", "r1", "l1"), "thetaA0", "kthetaA"))
        # Angles: r1-l1-l2 (thetaB0, kthetaB)
        output.append(write_angle(("r1", "l1", "l2"), "thetaB0", "kthetaB"))

        if restraint_lambda:
            output.append("[ dihedral_restraints ]")
            output.append(
                "; ai         aj         ak         al      type phiA       dphiA fcA       phiB       dphiB  fcB"
            )
        else:
            output.append("[ dihedrals ]")
            output.append(
                "; ai         aj         ak         al      type phiA       fcA        phiB       fcB"
            )
        # Dihedrals: r3-r2-r1-l1 (phiA0, kphiA)
        output.append(
            write_dihedral(
                ("r3", "r2", "r1", "l1"),
                "phiA0",
                "kphiA",
                restraint_lambda=restraint_lambda,
            )
        )
        # Dihedrals: r2-r1-l1-l2 (phiB0, kphiB)
        output.append(
            write_dihedral(
                ("r2", "r1", "l1", "l2"),
                "phiB0",
                "kphiB",
                restraint_lambda=restraint_lambda,
            )
        )
        # Dihedrals: r1-l1-l2-l3 (phiC0, kphiC)
        output.append(
            write_dihedral(
                ("r1", "l1", "l2", "l3"),
                "phiC0",
                "kphiC",
                restraint_lambda=restraint_lambda,
            )
        )

        return "\n".join(output)

    def _gromacs_multiple_distance(self, perturbation_type=None):
        """
        Format the Gromacs string for multiple distance restraints.

        Parameters
        ----------
        perturbation_type : str, optional, default=None
            The type of perturbation to applied during the current stage of the free energy
            calculation. If the perturbation type is "release_restraint", the permanent distance
            restraint will be written as a distance restraint (topology file directive [ distance
            restraints ], not affected by any lambda values), while all other restraints will be
            written as restraint potentials (topology file directive [ bonds ], affected by bonded-lambda.
            This allows the bond restraints to be released while retaining the permanent distance restraint.
            For all other perturbation types, all restraints will be written as restraint potential bonds.

        Returns
        -------
        str
            The Gromacs string for the multiple distance restraints.
        """

        def _get_distance_restraint_restraint_str(r1, l1, r0, r_fb, kr):
            """
            Get the text line specifying a distance restraint restraint
            (unaffected by any lambdas).
            """
            # Calculate parameters.
            ai = self._system.getIndex(r1) + 1
            aj = self._system.getIndex(l1) + 1
            low = r0 - r_fb
            up1 = r0 + r_fb
            up2 = 100  # Set this unresonably high so that we never get a linear potential.

            # Format strings, remembering to convert units.
            restr_type_str = "2"  # No averaging
            restr_index_str = "0"  # Only one restraint
            low_str = "{:.3f}".format(low / _nanometer)
            up1_str = "{:.3f}".format(up1 / _nanometer)
            up2_str = "{:.3f}".format(up2)
            fac_str = (
                "1.0"  # Scale the force constant (specified in the mdp file) by 1.0.
            )

            # Format entire string.
            # ai, aj, type, index, type', low, up1, up2, fac
            distance_restraint_restraint_parameter_string = f"  {ai:<10} {aj:<10} {restr_type_str:<10} {restr_index_str:<10} {restr_type_str:<10} {low_str:<10} {up1_str:<10} {up2_str:<10} {fac_str:<10}"

            return distance_restraint_restraint_parameter_string

        def _get_restraint_potential_bond_str(r1, l1, r0, r_fb, kr):
            """
            Get the text line specifying a restraint potential bond
            (affected by bonded-lambda).
            """
            # Calculate parameters.
            ai = self._system.getIndex(r1) + 1
            aj = self._system.getIndex(l1) + 1
            low = r0 - r_fb
            up1 = r0 + r_fb
            up2 = 100  # Set this unresonably high so that we never get a linear potential.

            # Format strings, remembering to convert units.
            restr_type_str = "10"  # Restraint potential bond
            low_str = "{:.3f}".format(low / _nanometer)
            up1_str = "{:.3f}".format(up1 / _nanometer)
            up2_str = "{:.3f}".format(up2)
            kdr_0_str = "{:.2f}".format(
                0
            )  # Force constant is 0 when bonded-lambda = 0.
            kdr_str = "{:.2f}".format(kr / (_kj_per_mol / _nanometer**2))

            # Format entire string.
            # ai, aj, type, lowA, up1A, up2A, kdrA, lowB, up1B, up2B, kdrB
            restraint_potential_bond_parameter_string = f"  {ai:<10} {aj:<10} {restr_type_str:<10} {low_str:<10} {up1_str:<10} {up2_str:<10} {kdr_0_str:<10} {low_str:<10} {up1_str:<10} {up2_str:<10} {kdr_str:<10}"

            return restraint_potential_bond_parameter_string

        # Write the output string.
        output = [
            "[ intermolecular_interactions ]",
        ]

        output.append("[ bonds ]")
        output.append(
            "; ai         aj         type       lowA       up1A       up2A       kdrA       lowB       up1B       up2B       kdrB"
        )
        for restraint in self._restraint_dict["distance_restraints"]:
            output.append(
                _get_restraint_potential_bond_str(
                    restraint["r1"],
                    restraint["l1"],
                    restraint["r0"],
                    restraint["r_fb"],
                    restraint["kr"],
                )
            )
        if perturbation_type != "release_restraint":
            output.append(
                _get_restraint_potential_bond_str(
                    self._restraint_dict["permanent_distance_restraint"]["r1"],
                    self._restraint_dict["permanent_distance_restraint"]["l1"],
                    self._restraint_dict["permanent_distance_restraint"]["r0"],
                    self._restraint_dict["permanent_distance_restraint"]["r_fb"],
                    self._restraint_dict["permanent_distance_restraint"]["kr"],
                )
            )
        else:  # Perturbation type is release_restraint - write the permanent restraint as a distance restraint
            output.append("[ distance_restraints ]")
            output.append(
                "; ai         aj         type       index      type'      low        up1        up2        fac"
            )
            output.append(
                _get_distance_restraint_restraint_str(
                    self._restraint_dict["permanent_distance_restraint"]["r1"],
                    self._restraint_dict["permanent_distance_restraint"]["l1"],
                    self._restraint_dict["permanent_distance_restraint"]["r0"],
                    self._restraint_dict["permanent_distance_restraint"]["r_fb"],
                    self._restraint_dict["permanent_distance_restraint"]["kr"],
                )
            )

        return "\n".join(output)

    def _somd_boresch(self, perturbation_type=None):
        """Format the SOMD string for the Boresch restraints."""

        # Indices
        r1 = self._system.getIndex(self._restraint_dict["anchor_points"]["r1"])
        r2 = self._system.getIndex(self._restraint_dict["anchor_points"]["r2"])
        r3 = self._system.getIndex(self._restraint_dict["anchor_points"]["r3"])
        l1 = self._system.getIndex(self._restraint_dict["anchor_points"]["l1"])
        l2 = self._system.getIndex(self._restraint_dict["anchor_points"]["l2"])
        l3 = self._system.getIndex(self._restraint_dict["anchor_points"]["l3"])
        # Equilibrium values
        r0 = self._restraint_dict["equilibrium_values"]["r0"] / _angstrom
        thetaA0 = self._restraint_dict["equilibrium_values"]["thetaA0"] / _radian
        thetaB0 = self._restraint_dict["equilibrium_values"]["thetaB0"] / _radian
        phiA0 = self._restraint_dict["equilibrium_values"]["phiA0"] / _radian
        phiB0 = self._restraint_dict["equilibrium_values"]["phiB0"] / _radian
        phiC0 = self._restraint_dict["equilibrium_values"]["phiC0"] / _radian
        # Force constants. Halve these as SOMD defines force constants as E = kx**2
        kr = (
            self._restraint_dict["force_constants"]["kr"]
            / (_kcal_per_mol / (_angstrom * _angstrom))
        ) / 2
        kthetaA = (
            self._restraint_dict["force_constants"]["kthetaA"]
            / (_kcal_per_mol / (_radian * _radian))
        ) / 2
        kthetaB = (
            self._restraint_dict["force_constants"]["kthetaB"]
            / (_kcal_per_mol / (_radian * _radian))
        ) / 2
        kphiA = (
            self._restraint_dict["force_constants"]["kphiA"]
            / (_kcal_per_mol / (_radian * _radian))
        ) / 2
        kphiB = (
            self._restraint_dict["force_constants"]["kphiB"]
            / (_kcal_per_mol / (_radian * _radian))
        ) / 2
        kphiC = (
            self._restraint_dict["force_constants"]["kphiC"]
            / (_kcal_per_mol / (_radian * _radian))
        ) / 2

        restr_string = f'boresch restraints dictionary = {{"anchor_points":{{"r1":{r1}, "r2":{r2}, "r3":{r3}, "l1":{l1}, '
        restr_string += f'"l2":{l2}, "l3":{l3}}}, '
        restr_string += f'"equilibrium_values":{{"r0":{r0:.2f}, "thetaA0":{thetaA0:.2f}, "thetaB0":{thetaB0:.2f},"phiA0":{phiA0:.2f}, '
        restr_string += f'"phiB0":{phiB0:.2f}, "phiC0":{phiC0:.2f}}}, '
        restr_string += f'"force_constants":{{"kr":{kr:.2f}, "kthetaA":{kthetaA:.2f}, "kthetaB":{kthetaB:.2f}, "kphiA":{kphiA:.2f}, '
        restr_string += f'"kphiB":{kphiB:.2f}, "kphiC":{kphiC:.2f}}}}}'

        return restr_string

    def _somd_multiple_distance(self, perturbation_type=None):
        """Format the SOMD string for the multiple distance restraints."""

        def _add_restr_to_str(restr, restr_string):
            """Apend the information for a single restraint to the string."""
            # Indices
            r1 = self._system.getIndex(restr["r1"])
            l1 = self._system.getIndex(restr["l1"])
            # Equilibrium value
            r0 = restr["r0"] / _angstrom
            # Force constant. Halve these as SOMD defines force constants as E = kx**2
            kr = (restr["kr"] / (_kcal_per_mol / (_angstrom * _angstrom))) / 2
            # Flat bottomed radius
            r_fb = restr["r_fb"] / _angstrom

            restr_string += f"({r1}, {l1}): ({r0}, {kr}, {r_fb}), "
            return restr_string

        standard_restr_string = "distance restraints dictionary = {"

        for single_restraint in self._restraint_dict["distance_restraints"]:
            standard_restr_string = _add_restr_to_str(
                single_restraint, standard_restr_string
            )

        if perturbation_type == "restraint":
            # In this case, we want all restraints to be switched on, even "permanent" ones
            standard_restr_string = _add_restr_to_str(
                self._restraint_dict["permanent_distance_restraint"],
                standard_restr_string,
            )
            return standard_restr_string[:-2] + "}"
        else:  # Other perturbation types, we want the permanent restraints to be constantly on
            standard_restr_string = standard_restr_string[:-2] + "}\n"
            permanent_restr_string = "permanent distance restraints dictionary = {"
            permanent_restr_string = _add_restr_to_str(
                self._restraint_dict["permanent_distance_restraint"],
                permanent_restr_string,
            )
            return standard_restr_string + permanent_restr_string[:-2] + "}"

    def toString(self, engine, perturbation_type=None, restraint_lambda=False):
        """
        The method for convert the restraint to a format that could be used
        by MD Engines.

        Parameters
        ----------

        engine : str
            The molecular dynamics engine for which to generate the restraint.
            Available options are currently "GROMACS" or "SOMD" for Boresch restraints,
            or "SOMD" only for multiple distance restraints.
        perturbation_type : str, optional, default=None
            The type of perturbation to applied during the current stage of the free energy
            calculation. This is only used for multiple distance restraints, for which all
            restraints are converted to standard distance restraints to allow them to be
            turned on when the perturbation type is "restraint", but for which the permanent
            distance restraint is always active if the perturbation type is "release_restraint"
            (or any other perturbation type).
        restraint_lambda : str, optional, default=False
            Whether to use restraint_lambda in Gromacs, this would move the dihedral restraints
            from [ dihedrals ], which is controlled by the bonded-lambda to
            [ dihedral_restraints ], which is controlled by restraint-lambda.
        """
        engine = engine.strip().lower()
        match (self._restraint_type, engine):
            case "boresch", "gromacs":
                return self._gromacs_boresch(
                    perturbation_type, restraint_lambda=restraint_lambda
                )
            case "boresch", "somd":
                return self._somd_boresch(perturbation_type)
            case "multiple_distance", "gromacs":
                return self._gromacs_multiple_distance(perturbation_type)
            case "multiple_distance", "somd":
                return self._somd_multiple_distance(perturbation_type)
            case _:
                raise NotImplementedError(
                    f"Restraint type {self._restraint_type} not implemented "
                    f"yet for {engine}."
                )

    def getCorrection(self, method="analytical"):
        """
        Calculate the free energy of releasing the restraint
        to the standard state volume.

        Parameters
        ----------
        method : str
            The integration method to use for calculating the correction for
            releasing the restraint to the standard state concentration.
            "numerical" or "analytical". Note that the Boresch analytical
            correction can introduce errors when the restraints are weak,
            restrained angles are close to 0 or pi radians, or the restrained
            distance is close to 0.

        Returns
        ----------
        dG : float
            Free energy of releasing the restraint to the standard state volume,
            in kcal / mol.
        """
        # Constants. Take .value() to avoid issues with ** and log of GeneralUnit
        v0 = (
            ((_Sire_meter3 / 1000) / _Sire_mole) / _Sire_angstrom3
        ).value()  # standard state volume in A^3
        R = (
            _k_boltz.value() * _kcal_per_mol / _kelvin
        ).value()  # molar gas constant in kcal mol-1 K-1

        # Parameters
        T = self.T / _kelvin  # Temperature in Kelvin

        if self._restraint_type == "boresch":
            # Constants. Take .value() to avoid issues with ** and log of GeneralUnit
            v0 = (
                ((_Sire_meter3 / 1000) / _Sire_mole) / _Sire_angstrom3
            ).value()  # standard state volume in A^3
            R = (
                _k_boltz.value() * _kcal_per_mol / _kelvin
            ).value()  # molar gas constant in kcal mol-1 K-1

            # Parameters
            prefactor = (
                8 * (_np.pi**2) * v0
            )  # In A^3. Divide this to account for force constants of 0 in the
            # analytical correction

            if method == "numerical":
                # ========= Acknowledgement ===============
                # Calculation copied from restraints.py  in
                # Yank https://github.com/choderalab/yank
                # =========================================

                def numerical_distance_integrand(r, r0, kr):
                    """Integrand for harmonic distance restraint. Domain is on [0, infinity],
                    but this will be truncated to [0, 8 RT] for practicality.

                    Parameters
                    ----------
                        r : (float)
                            Distance to be integrated, in Angstrom
                        r0 : (float)
                            Equilibrium distance, in Angstrom
                        kr : (float)
                            Force constant, in kcal mol-1 A-2

                    Returns
                    ----------
                        float : Value of integrand
                    """
                    return (r**2) * _np.exp(-(kr * (r - r0) ** 2) / (2 * R * T))

                def numerical_angle_integrand(theta, theta0, ktheta):
                    """Integrand for harmonic angle restraints. Domain is on [0,pi].

                    Parameters
                    ----------
                        theta : (float)
                        Angle to be integrated, in radians
                        theta0 : (float)
                        Equilibrium angle, in radians
                        ktheta : (float)
                        Force constant, in kcal mol-1 rad-2

                    Returns
                    ----------
                        float: Value of integrand
                    """
                    return _np.sin(theta) * _np.exp(
                        -(ktheta * (theta - theta0) ** 2) / (2 * R * T)
                    )

                def numerical_dihedral_integrand(phi, phi0, kphi):
                    """Integrand for the harmonic dihedral restraints. Domain is on [-pi,pi].

                    Parameters
                    ----------
                        phi : (float)
                        Angle to be integrated, in radians
                        phi0 : (float)
                        Equilibrium angle, in radians
                        kphi : (float)
                        Force constant, in kcal mol-1 rad-2

                    Returns
                    ----------
                        float: Value of integrand
                    """
                    d_phi = abs(phi - phi0)
                    d_phi_corrected = min(
                        d_phi, 2 * _np.pi - d_phi
                    )  # correct for periodic boundaries
                    return _np.exp(-(kphi * d_phi_corrected**2) / (2 * R * T))

                # Radial
                r0 = self._restraint_dict["equilibrium_values"]["r0"] / _angstrom  # A
                kr = self._restraint_dict["force_constants"]["kr"] / (
                    _kcal_per_mol / _angstrom2
                )  # kcal mol-1 A-2
                dist_at_8RT = 4 * _np.sqrt(
                    (R * T) / kr
                )  # Dist. which gives restraint energy = 8 RT
                r_min = max(0, r0 - dist_at_8RT)
                r_max = r0 + dist_at_8RT
                integrand = lambda r: numerical_distance_integrand(r, r0, kr)
                z_r = _integrate.quad(integrand, r_min, r_max)[0]

                # Angular
                for angle in ["thetaA", "thetaB"]:
                    theta0 = (
                        self._restraint_dict["equilibrium_values"][f"{angle}0"]
                        / _radian
                    )  # rad
                    ktheta = self._restraint_dict["force_constants"][f"k{angle}"] / (
                        _kcal_per_mol / (_radian * _radian)
                    )  # kcal mol-1 rad-2
                    integrand = lambda theta: numerical_angle_integrand(
                        theta, theta0, ktheta
                    )
                    z_r *= _integrate.quad(integrand, 0, _np.pi)[0]

                # Dihedral
                for dihedral in ["phiA", "phiB", "phiC"]:
                    phi0 = (
                        self._restraint_dict["equilibrium_values"][f"{dihedral}0"]
                        / _radian
                    )  # rad
                    kphi = self._restraint_dict["force_constants"][f"k{dihedral}"] / (
                        _kcal_per_mol / (_radian * _radian)
                    )  # kcal mol-1 rad-2
                    integrand = lambda phi: numerical_dihedral_integrand(
                        phi, phi0, kphi
                    )
                    z_r *= _integrate.quad(integrand, -_np.pi, _np.pi)[0]

                # Compute dg and attach unit
                dg = -R * T * _np.log(prefactor / z_r)
                dg *= _kcal_per_mol

                return dg

            elif method == "analytical":
                # Only need three equilibrium values for the analytical correction
                r0 = (
                    self._restraint_dict["equilibrium_values"]["r0"] / _angstrom
                )  # Distance in A
                thetaA0 = (
                    self._restraint_dict["equilibrium_values"]["thetaA0"] / _radian
                )  # Angle in radians
                thetaB0 = (
                    self._restraint_dict["equilibrium_values"]["thetaB0"] / _radian
                )  # Angle in radians

                force_constants = []

                # Loop through and correct for force constants of zero,
                # which break the analytical correction. To account for this,
                # divide the prefactor accordingly. Note that setting
                # certain force constants to zero while others are non-zero
                # will result in unstable restraints, but this will be checked when
                # the restraint object is created
                for k, val in self._restraint_dict["force_constants"].items():
                    if val.value() == 0:
                        if k == "kr":
                            raise ValueError("The force constant kr must not be zero")
                        if k == "kthetaA":
                            prefactor /= 2 / _np.sin(thetaA0)
                        if k == "kthetaB":
                            prefactor /= 2 / _np.sin(thetaB0)
                        if k[:4] == "kphi":
                            prefactor /= 2 * _np.pi
                    else:
                        if k == "kr":
                            force_constants.append(val / (_kcal_per_mol / _angstrom2))
                        else:
                            force_constants.append(
                                val / (_kcal_per_mol / (_radian * _radian))
                            )

                # Calculation
                n_nonzero_k = len(force_constants)
                prod_force_constants = _np.prod(force_constants)
                numerator = prefactor * _np.sqrt(prod_force_constants)
                denominator = (
                    (r0**2)
                    * _np.sin(thetaA0)
                    * _np.sin(thetaB0)
                    * (2 * _np.pi * R * T) ** (n_nonzero_k / 2)
                )

                # Compute dg and attach unit
                dg = -R * T * _np.log(numerator / denominator)
                dg *= _kcal_per_mol

                return dg

            else:
                raise ValueError(
                    f"Correction method {method} is not supported. Please choose from 'numerical' or 'analytical'."
                )

        if self._restraint_type == "multiple_distance":
            if method == "analytical":
                raise NotImplementedError(
                    "The analytical correction is not supported for multiple distance restraints."
                )

            else:

                def _numerical_distance_integrand(r, r0, r_fb, kr):
                    """
                    Integrand for harmonic distance restraint.

                    Parameters
                    ----------
                    r : float
                        Distance to be integrated, in Angstrom
                    r0 : float
                        Equilibrium distance, in Angstrom
                    r_fb : float
                        Flat-bottomed radius, in Angstrom
                    kr : float
                        Force constant, in kcal mol-1 A-2

                    Returns
                    -------
                    float
                        Value of integrand

                    Notes
                    -----
                    The domain of the integrand is [0, infinity], but this will
                    be truncated to [0, 8 RT] for practicality.
                    """
                    r_eff = abs(r - r0) - r_fb
                    if r_eff < 0:
                        r_eff = 0
                    return (r**2) * _np.exp(-(kr * r_eff**2) / (2 * R * T))

                def _get_correction(r0, r_fb, kr):
                    """
                    Get the free energy of releasing the harmonic distance restraint.

                    Parameters
                    ----------
                    r0 : float
                        Equilibrium distance, in Angstrom
                    r_fb : float
                        Flat-bottomed radius, in Angstrom
                    kr : float
                        Force constant, in kcal mol-1 A-2

                    Returns
                    -------
                    float
                        Free energy of releasing the restraint

                    Notes
                    -----
                    The domain of the integrand is [0, infinity], but this will
                    be truncated to [0, 8 RT] for practicality.
                    """
                    dist_at_8RT = (
                        4 * _np.sqrt((R * T) / kr) + r_fb
                    )  # Dist. which gives restraint energy = 8 RT
                    r_min = max(0, r0 - dist_at_8RT)
                    r_max = r0 + dist_at_8RT
                    integrand = lambda r: _numerical_distance_integrand(r, r0, r_fb, kr)
                    z_r = _integrate.quad(integrand, r_min, r_max)[0]
                    dg = -R * T * _np.log(v0 / (4 * _np.pi * z_r))

                    # Attatch unit of kcal/mol
                    dg *= _kcal_per_mol

                    return dg

                # Get the parameters from the permanent distance restraint, which is not released
                r0 = (
                    self._restraint_dict["permanent_distance_restraint"]["r0"]
                    / _angstrom
                )
                r_fb = (
                    self._restraint_dict["permanent_distance_restraint"]["r_fb"]
                    / _angstrom
                )
                kr = self._restraint_dict["permanent_distance_restraint"]["kr"] / (
                    _kcal_per_mol / _angstrom2
                )
                _warnings.warn(
                    "The multiple distance restraint correction is assumes that only "
                    "the 'permanent_distance_restraint' is active."
                )
                return _get_correction(r0, r_fb, kr)

    @property
    def correction(self):
        """Give the free energy of removing the restraint."""
        method = "analytical" if self._restraint_type == "boresch" else "numerical"
        return self.getCorrection(method=method)
