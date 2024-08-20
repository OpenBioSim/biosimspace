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
__all__ = ["_AToMUtils"]
from .. import Protocol as _Protocol
from ..Types import Vector as _Vector
import math as _math
import warnings as _warnings


class _AToMUtils:
    # Internal class for creating openmm forces within an AToM process.
    def __init__(self, protocol):
        # Check for proper typing
        if not isinstance(
            protocol,
            (
                _Protocol.AToMMinimisation,
                _Protocol.AToMEquilibration,
                _Protocol.AToMAnnealing,
                _Protocol.AToMProduction,
            ),
        ):
            raise TypeError("Protocol must be an AToM protocol")
        self.protocol = protocol
        self.data = self.protocol.getData()

    def getAlignmentConstants(self):
        self.alignment_k_distance = self.protocol.getAlignKfSep().value()
        self.alignment_k_theta = self.protocol.getAlignKTheta().value()
        self.alignment_k_psi = self.protocol.getAlignKPsi().value()

    def getCMConstants(self):
        self.cm_kf = self.protocol.getCMKf().value()
        self.cm_tol = self.protocol.getCMTol().value()

    def findAbsoluteCoreIndices(self):
        import numpy as np

        self.lig1_first_atomnum = self.data["first_ligand1_atom_index"]
        self.lig1_rigid_atoms = list(
            np.add(self.lig1_first_atomnum, self.data["ligand1_rigid_core"])
        )
        self.lig2_first_atomnum = self.data["first_ligand2_atom_index"]
        self.lig2_rigid_atoms = list(
            np.add(self.lig2_first_atomnum, self.data["ligand2_rigid_core"])
        )

    def findAbsoluteCOMAtoms(self):
        import numpy as np

        self.protein_first_atomnum = self.data["first_protein_atom_index"]
        self.protein_com_atoms = list(
            np.add(self.protein_first_atomnum, self.data["protein_com_atoms"])
        )

        self.lig1_first_atomnum = self.data["first_ligand1_atom_index"]
        self.lig1_com_atoms = list(
            np.add(self.lig1_first_atomnum, self.data["ligand1_com_atoms"])
        )

        self.lig2_first_atomnum = self.data["first_ligand2_atom_index"]
        self.lig2_com_atoms = list(
            np.add(self.lig2_first_atomnum, self.data["ligand2_com_atoms"])
        )

    def getATMForceConstants(self, index=None):
        self.lig1_atoms = self.getLigand1AtomsAsList()
        self.lig2_atoms = self.getLigand2AtomsAsList()
        self.SCUmax = self.protocol.getSCUmax().value()
        self.SCU0 = self.protocol.getSCU0().value()
        self.SCa = self.protocol.getSCa()
        if isinstance(self.protocol, _Protocol.AToMProduction):
            if index is None:
                raise ValueError("Index must be set for AToMProduction protocol")
            self.lambda1 = self.protocol.getLambda1()[index]
            self.lambda2 = self.protocol.getLambda2()[index]
            self.alpha = self.protocol.getAlpha()[index].value()
            self.uh = self.protocol.getUh()[index].value()
            self.w0 = self.protocol.getW0()[index].value()
            self.direction = self.protocol.getDirection()[index]
            self.master_lambda = self.protocol.get_lambda_values()[index]
        elif isinstance(
            self.protocol, (_Protocol.AToMEquilibration, _Protocol.AToMAnnealing)
        ):
            self.lambda1 = self.protocol.getLambda1()
            self.lambda2 = self.protocol.getLambda2()
            self.alpha = self.protocol.getAlpha().value()
            self.uh = self.protocol.getUh().value()
            self.w0 = self.protocol.getW0().value()
            self.direction = self.protocol.getDirection()

    def _dump_atm_constants_to_dict(self):
        """Internal function to write all ATM window-dependent constants to a dictionary (string)
        to be used in sampling for analysis."""
        output = ""
        output += "atm_constants = {\n"
        output += "    'Lambda1': {},\n".format(self.protocol.getLambda1())
        output += "    'Lambda2': {},\n".format(self.protocol.getLambda2())
        output += "    'Alpha': {},\n".format(
            [i.value() for i in self.protocol.getAlpha()]
        )
        output += "    'Uh': {},\n".format([i.value() for i in self.protocol.getUh()])
        output += "    'W0': {},\n".format([i.value() for i in self.protocol.getW0()])
        output += "    'Direction': {}\n".format(self.protocol.getDirection())
        output += "}\n"

        output += "for key in atm_constants.keys():\n"
        output += "    if key in ['Alpha','Uh','W0']:\n"
        output += "        atm_constants[key] = [i for i in atm_constants[key] * kilocalories_per_mole]\n"

        return output

    def findDisplacement(self):
        d = self.data["displacement"]
        if isinstance(d, (list)):
            if not all(isinstance(x, float) for x in d):
                raise TypeError("Displacement must be a list of floats")
            self.displacement = d
        elif isinstance(d, _Vector):
            disp = [d.x(), d.y(), d.z()]
            self.displacement = disp
        else:
            raise TypeError("Displacement must be a list or BioSimSpace vector")

    def createDisplacement(self):
        self.findDisplacement()
        d = [round(x, 3) for x in self.displacement]
        output = ""
        output += "displacement = {}\n".format(d)
        output += "#BioSimSpace output is in angstrom, divide by 10 to convert to the expected units of nm\n"
        output += "displacement = [i/10.0 for i in displacement]\n"
        return output

    def createSoftcorePertE(self):
        """Create the softcorePertE function for the Gallachio lab analysis"""
        output = ""
        output += "def softCorePertE(u, umax, ub, a):\n"
        output += "    usc = u\n"
        output += "    if u > ub:\n"
        output += "        gu = (u-ub)/(a*(umax-ub))\n"
        output += "        zeta = 1. + 2.*gu*(gu + 1.)\n"
        output += "        zetap = np.power( zeta, a)\n"
        output += "        usc = (umax-ub)*(zetap - 1.)/(zetap + 1.) + ub\n"
        output += "    return usc\n"
        return output

    def createAlignmentForce(self, force_group=None):
        """
        Create the alignment force that keeps the ligands co-planar.

        parameters
        ----------
        force_group : None or list
            Group of the force to be added to the system. If none defined then no force group will be set
            (therefore it will default to 0). Only tested for single-point energies.
            If a list is given the groups will be assigned in the order [distance, angle, dihedral]
        """
        # This force is the same in every lambda window
        self.getAlignmentConstants()
        self.findAbsoluteCoreIndices()

        if force_group is not None and len(force_group) != 3:
            raise ValueError("Force group must be a list of three integers")
        output = "\n\n"
        output += "k_distance = {} * kilocalorie_per_mole / angstrom**2\n".format(
            self.alignment_k_distance
        )
        output += "k_theta = {} * kilocalorie_per_mole\n".format(self.alignment_k_theta)
        output += "k_psi = {} * kilocalorie_per_mole\n".format(self.alignment_k_psi)
        output += "idxs_a = {}\n".format(self.lig1_rigid_atoms)
        output += "idxs_b = {}\n".format(self.lig2_rigid_atoms)
        output += "\n\n"

        output += 'distance_energy_fn = "0.5 * k * ((x1 - x2 - dx)^2 + (y1 - y2 - dy)^2 + (z1 - z2 - dz)^2);"\n'
        output += "distance_force = CustomCompoundBondForce(2, distance_energy_fn)\n"
        output += "distance_force.addPerBondParameter('k')\n"
        output += "distance_force.addPerBondParameter('dx')\n"
        output += "distance_force.addPerBondParameter('dy')\n"
        output += "distance_force.addPerBondParameter('dz')\n"

        output += """distance_parameters = [
        k_distance.value_in_unit(kilojoules_per_mole / nanometer**2),
        displacement[0]*nanometer,
        displacement[1]*nanometer,
        displacement[2]*nanometer,
        ]\n"""

        output += (
            "distance_force.addBond((idxs_b[0], idxs_a[0]), distance_parameters)\n"
        )
        if force_group is not None:
            output += "distance_force.setForceGroup({})\n".format(force_group[0])
        output += "system.addForce(distance_force)\n"
        output += "\n\n"

        output += """angle_energy_fn = (
        "0.5 * k * (1 - cos_theta);"
        ""
        "cos_theta = (dx_1 * dx_2 + dy_1 * dy_2 + dz_1 * dz_2) / (norm_1 * norm_2);"
        ""
        "norm_1 = sqrt(dx_1^2 + dy_1^2 + dz_1^2);"
        "dx_1 = x2 - x1; dy_1 = y2 - y1; dz_1 = z2 - z1;"
        ""
        "norm_2 = sqrt(dx_2^2 + dy_2^2 + dz_2^2);"
        "dx_2 = x4 - x3; dy_2 = y4 - y3; dz_2 = z4 - z3;"
        )\n"""
        output += "angle_force = CustomCompoundBondForce(4, angle_energy_fn)\n"
        output += 'angle_force.addPerBondParameter("k")\n'
        output += """angle_force.addBond(
        (idxs_b[0], idxs_b[1], idxs_a[0], idxs_a[1]),
        [k_theta.value_in_unit(kilojoules_per_mole)],
        )\n"""
        if force_group is not None:
            output += "angle_force.setForceGroup({})\n".format(force_group[1])
        output += "system.addForce(angle_force)\n\n"

        # Femto dihedral form:
        # output += """dihedral_energy_fn = (
        # "0.5 * k * (1 - cos_phi);"
        # ""
        # "cos_phi = (v_x * w_x + v_y * w_y + v_z * w_z) / (norm_v * norm_w);"
        # ""
        # "norm_v = sqrt(v_x^2 + v_y^2 + v_z^2);"
        # "v_x = dx_31 - dot_31 * dx_21 / norm_21;"
        # "v_y = dy_31 - dot_31 * dy_21 / norm_21;"
        # "v_z = dz_31 - dot_31 * dz_21 / norm_21;"
        # ""
        # "dot_31 = (dx_31 * dx_21 + dy_31 * dy_21 + dz_31 * dz_21) / norm_21;"
        # "dx_31 = x3 - x1; dy_31 = y3 - y1; dz_31 = z3 - z1;"
        # ""
        # "norm_w = sqrt(w_x^2 + w_y^2 + w_z^2);"
        # "w_x = dx_54 - dot_54 * dx_21 / norm_21;"
        # "w_y = dy_54 - dot_54 * dy_21 / norm_21;"
        # "w_z = dz_54 - dot_54 * dz_21 / norm_21;"
        # ""
        # "dot_54 =(dx_54 * dx_21 + dy_54 * dy_21 + dz_54 * dz_21) / norm_21;"
        # "dx_54 = x5 - x4; dy_54 = y5 - y4; dz_54 = z5 - z4;"
        # ""
        # "norm_21 = sqrt(dx_21^2 + dy_21^2 + dz_21^2);"
        # "dx_21 = x2 - x1; dy_21 = y2 - y1; dz_21 = z2 - z1;"
        # )\n"""

        # Gallachio lab dihedral form:
        output += 'dihedral_energy_fn = "(k/2) * (1 - cosp) ; "\n'
        output += 'dihedral_energy_fn += "cosp = xvn*xwn + yvn*ywn + zvn*zwn ; "\n'
        output += 'dihedral_energy_fn += "xvn = xv/v ; yvn = yv/v; zvn = zv/v ;"\n'
        output += 'dihedral_energy_fn += "v = sqrt(xv^2 + yv^2 + zv^2 ) ;"\n'
        output += 'dihedral_energy_fn += "xv = xd0 - dot01*xdn1 ;"\n'
        output += 'dihedral_energy_fn += "yv = yd0 - dot01*ydn1 ;"\n'
        output += 'dihedral_energy_fn += "zv = zd0 - dot01*zdn1 ;"\n'
        output += 'dihedral_energy_fn += "dot01 = xd0*xdn1 + yd0*ydn1 + zd0*zdn1 ;"\n'
        output += 'dihedral_energy_fn += "xd0 = x3 - x1 ;"\n'
        output += 'dihedral_energy_fn += "yd0 = y3 - y1 ;"\n'
        output += 'dihedral_energy_fn += "zd0 = z3 - z1 ;"\n'
        output += 'dihedral_energy_fn += "xwn = xw/w ; ywn = yw/w; zwn = zw/w ;"\n'
        output += 'dihedral_energy_fn += "w = sqrt(xw^2 + yw^2 + zw^2) ;"\n'
        output += 'dihedral_energy_fn += "xw = xd3 - dot31*xdn1 ;"\n'
        output += 'dihedral_energy_fn += "yw = yd3 - dot31*ydn1 ;"\n'
        output += 'dihedral_energy_fn += "zw = zd3 - dot31*zdn1 ;"\n'
        output += 'dihedral_energy_fn += "dot31 = xd3*xdn1 + yd3*ydn1 + zd3*zdn1 ;"\n'
        output += 'dihedral_energy_fn += "xd3 = x5 - x4 ;"\n'
        output += 'dihedral_energy_fn += "yd3 = y5 - y4 ;"\n'
        output += 'dihedral_energy_fn += "zd3 = z5 - z4 ;"\n'
        output += 'dihedral_energy_fn += "xdn1 = xd1/dn1 ; ydn1 = yd1/dn1 ; zdn1 = zd1/dn1 ;"\n'
        output += 'dihedral_energy_fn += "dn1 = sqrt(xd1^2 + yd1^2 + zd1^2) ;"\n'
        output += 'dihedral_energy_fn += "xd1 = x2 - x1 ;"\n'
        output += 'dihedral_energy_fn += "yd1 = y2 - y1 ;"\n'
        output += 'dihedral_energy_fn += "zd1 = z2 - z1 ;"\n'

        output += "dihedral_force = CustomCompoundBondForce(5, dihedral_energy_fn)\n"
        output += 'dihedral_force.addPerBondParameter("k")\n'
        output += """dihedral_force.addBond(
        (idxs_b[0], idxs_b[1], idxs_b[2], idxs_a[0], idxs_a[2]),
        [0.5 * k_psi.value_in_unit(kilojoules_per_mole)],
        )\n"""
        output += """dihedral_force.addBond(
        (idxs_a[0], idxs_a[1], idxs_a[2], idxs_b[0], idxs_b[2]),
        [0.5 * k_psi.value_in_unit(kilojoules_per_mole)],
        )\n"""
        if force_group is not None:
            output += "dihedral_force.setForceGroup({})\n".format(force_group[2])
        output += "system.addForce(dihedral_force)\n\n"
        return output

    def getLigand1AtomsAsList(self):
        import numpy as np

        return list(
            np.arange(
                self.data["first_ligand1_atom_index"],
                self.data["last_ligand1_atom_index"] + 1,
            )
        )

    def getLigand2AtomsAsList(self):
        import numpy as np

        return list(
            np.arange(
                self.data["first_ligand2_atom_index"],
                self.data["last_ligand2_atom_index"] + 1,
            )
        )

    def createATMForce(
        self,
        index,
        force_group=None,
    ):
        """
        Create a string which can be added directly to an openmm script to add an ATM force to the system.

        Parameters
        ----------
        index : int
            Index of current window - used to set window-dependent variables.
        force_group : int
            Group of the force to be added to the system. Shuld only be needed when testing single-point energies.
        """
        self.findDisplacement()
        self.getATMForceConstants(index)
        output = ""
        output += "#Parameters for ATM force in  original units\n"
        output += "lig1_atoms = {}\n".format(self.lig1_atoms)
        output += "lig2_atoms = {}\n".format(self.lig2_atoms)
        if isinstance(self.protocol, _Protocol.AToMProduction):
            output += "window_index = {}\n".format(index)
        output += "lambda1 = {}\n".format(self.lambda1)
        output += "lambda2 = {}\n".format(self.lambda2)
        output += "alpha = {} * kilocalories_per_mole\n".format(self.alpha)
        output += "uh = {} * kilocalories_per_mole\n".format(self.uh)
        output += "w0 = {} * kilocalories_per_mole\n".format(self.w0)
        output += "direction = {}\n".format(self.direction)
        output += "sc_Umax = {} * kilocalories_per_mole\n".format(self.SCUmax)
        output += "sc_U0 = {} * kilocalories_per_mole\n".format(self.SCU0)
        output += "sc_a = {}\n".format(self.SCa)

        if isinstance(self.protocol, _Protocol.AToMProduction):
            output += self._dump_atm_constants_to_dict()

        output += "\n\n #Define ATM force\n"
        output += """atm_force = ATMForce(
        lambda1,
        lambda2,
        alpha.value_in_unit(kilojoules_per_mole),
        uh.value_in_unit(kilojoules_per_mole),
        w0.value_in_unit(kilojoules_per_mole),
        sc_Umax.value_in_unit(kilojoules_per_mole),
        sc_U0.value_in_unit(kilojoules_per_mole),
        sc_a,
        direction,
        )"""

        output += "\n\n #Add ATM force to system\n"
        output += "for _ in prm.topology.atoms():\n"
        output += "    atm_force.addParticle(Vec3(0.0,0.0,0.0))"
        output += "\n"
        # TODO: add offset - check convesion of a list to a Vec3
        # Assuming that offset is the 3-vector which deifnes the ligand displacement
        # need to convert displacement to nm
        output += "for i in lig1_atoms:\n"
        output += "    atm_force.setParticleParameters(i, Vec3(*displacement))\n"
        output += "for i in lig2_atoms:\n"
        output += "    atm_force.setParticleParameters(i, -Vec3(*displacement))\n"
        output += "\n"
        output += "nonbonded_force_id = [i for i, force in enumerate(system.getForces()) if isinstance(force, NonbondedForce)][0]\n"
        output += "nonbonded = copy.deepcopy(system.getForce(nonbonded_force_id))\n"
        output += "system.removeForce(nonbonded_force_id)\n"
        output += "atm_force.addForce(nonbonded)\n"
        if force_group is not None:
            output += "atm_force.setForceGroup({})\n".format(force_group)
        output += "system.addForce(atm_force)\n"
        return output

    def createCOMRestraint(self, force_group=None):
        """
        Create a string containing the CM-CM restriants for two groups of atoms.
        In most cases these will be some combination of protein and ligand atoms.
        Constants for the force are set in the protocol.

        parameters
        ----------
        force_group : None or int
            Group of the force to be added to the system. If none defined then no force group will be set
            (therefore it will default to 0). Only tested for single-point energies."""

        self.findAbsoluteCOMAtoms()
        # Groups contained within the constraint
        protein_com = self.protein_com_atoms
        lig1_com = self.lig1_com_atoms
        lig2_com = self.lig2_com_atoms
        self.getCMConstants()
        # Constants for the force
        kf_cm = self.cm_kf
        tol_cm = self.cm_tol
        output = ""
        output += "protein_com = {}\n".format(protein_com)
        output += "lig1_com = {}\n".format(lig1_com)
        output += "lig2_com = {}\n".format(lig2_com)
        output += "# Constants for the CM-CM force in their input units\n"
        output += "kfcm = {} * kilocalorie_per_mole / angstrom**2\n".format(kf_cm)
        output += "tolcm = {} * angstrom \n".format(tol_cm)

        # Add expression for cm restraint
        output += 'expr = "0.5 * kfcm * step(dist - tolcm) * (dist - tolcm)^2;""dist = sqrt((x1 - x2 - offx)^2 + (y1 - y2 - offy)^2 + (z1 - z2 - offz)^2);"\n'
        output += "force_CMCM = CustomCentroidBondForce(2, expr)\n"
        output += "force_CMCM.addPerBondParameter('kfcm')\n"
        output += "force_CMCM.addPerBondParameter('tolcm')\n"
        output += "force_CMCM.addPerBondParameter('offx')\n"
        output += "force_CMCM.addPerBondParameter('offy')\n"
        output += "force_CMCM.addPerBondParameter('offz')\n"
        output += "force_CMCM.addGroup(protein_com)\n"
        output += "force_CMCM.addGroup(lig1_com)\n"

        output += """parameters_bound = (
        kfcm.value_in_unit(kilojoules_per_mole / nanometer**2),
        tolcm.value_in_unit(nanometer),
        0.0 * nanometer,
        0.0 * nanometer,
        0.0 * nanometer,
        )\n"""
        output += "force_CMCM.addBond((1,0), parameters_bound)\n"
        output += "numgroups = force_CMCM.getNumGroups()\n"

        output += "force_CMCM.addGroup(protein_com)\n"
        output += "force_CMCM.addGroup(lig2_com)\n"
        output += """parameters_free = (
        kfcm.value_in_unit(kilojoules_per_mole / nanometer**2),
        tolcm.value_in_unit(nanometer),
        displacement[0] * nanometer,
        displacement[1] * nanometer,
        displacement[2] * nanometer,
        )\n"""

        output += "force_CMCM.addBond((numgroups+1,numgroups+0), parameters_free)\n"
        if force_group is not None:
            if not isinstance(force_group, int):
                raise TypeError("Force group must be an integer")
            output += "force_CMCM.setForceGroup({})\n".format(force_group)
        output += "system.addForce(force_CMCM)\n"
        output += "#End of CM-CM force\n\n"

        return output

    def create_flat_bottom_restraint(self, restrained_atoms, force_group=None):
        """Flat bottom restraint for atom-compatible position restraints

        Parameters
        ----------
        restrained_atoms : list
            List of atom indices to be restrained. Need to be explicitly given due to the ability to parse strings in the protocol.

        force_group : None or int
            Group of the force to be added to the system. If none defined then no force group will be set
            (therefore it will default to 0). Only tested for single-point energies.
        """
        # Still using the position restraint mixin, get the values of the relevant constants
        pos_const = self.protocol.getForceConstant().value()
        pos_width = self.protocol.getPosRestWidth().value()
        output = ""
        output += "fc = {} * kilocalorie_per_mole / angstrom**2\n".format(pos_const)
        output += "tol = {} * angstrom\n".format(pos_width)
        output += "restrained_atoms = {}\n".format(restrained_atoms)
        output += "positions = prm.positions\n"
        output += 'posrestforce = CustomExternalForce("0.5*fc*select(step(dist-tol), (dist-tol)^2, 0); dist = periodicdistance(x,y,z,x0,y0,z0)")\n'

        output += 'posrestforce.addPerParticleParameter("x0")\n'
        output += 'posrestforce.addPerParticleParameter("y0")\n'
        output += 'posrestforce.addPerParticleParameter("z0")\n'
        output += 'posrestforce.addPerParticleParameter("fc")\n'
        output += 'posrestforce.addPerParticleParameter("tol")\n'

        output += "for i in restrained_atoms:\n"
        output += "    x1 = positions[i][0].value_in_unit_system(openmm.unit.md_unit_system)\n"
        output += "    y1 = positions[i][1].value_in_unit_system(openmm.unit.md_unit_system)\n"
        output += "    z1 = positions[i][2].value_in_unit_system(openmm.unit.md_unit_system)\n"
        output += "    fc1 = fc.value_in_unit(kilojoules_per_mole / nanometer**2)\n"
        output += "    tol1 = tol.value_in_unit(nanometer)\n"
        output += "    posrestforce.addParticle(i, [x1, y1, z1, fc1, tol1])\n"
        if force_group is not None:
            if not isinstance(force_group, int):
                raise TypeError("Force group must be an integer")
            output += "posrestforce.setForceGroup({})\n".format(force_group)
        output += "system.addForce(posrestforce)\n"
        return output

    def createAnnealingProtocol(self):
        """
        Create a string which can be added directly to an openmm script to add an annealing protocol to the system.
        """
        anneal_runtime = self.protocol.getRunTime()
        num_cycles = self.protocol.getAnnealNumCycles()
        cycle_numsteps = int(
            (anneal_runtime / num_cycles) / self.protocol.getTimeStep()
        )

        prot = self.protocol.getAnnealValues()
        # Find all entries whose keys contain "start" and create a dictionary of these entries
        # Also remove the word "start" from the key
        start = {k.replace("_start", ""): v for k, v in prot.items() if "start" in k}
        # Same for "end"
        end = {k.replace("_end", ""): v for k, v in prot.items() if "end" in k}
        # write protocol to output in dictionary format
        output = ""
        output += f"values_start = {start}\n"
        output += f"values_end = {end}\n"
        output += "increments = {\n"
        output += f"    key: (values_end[key] - values_start[key]) / {num_cycles}\n"
        output += "    for key in values_start.keys()\n"
        output += "}\n"
        # First set all values using the start values
        output += "for key in values_start.keys():\n"
        output += "    simulation.context.setParameter(key, values_start[key])\n"
        # Now perform the annealing in cycles
        output += f"for i in range({int(num_cycles)}):\n"
        output += f"    simulation.step({cycle_numsteps})\n"
        output += "    print(f'Cycle {i+1}')\n"
        output += "    state = simulation.context.getState(getPositions=True, getVelocities=True)\n"
        output += "    for key in values_start.keys():\n"
        output += "        simulation.context.setParameter(key, simulation.context.getParameter(key) + increments[key])\n"
        output += "simulation.saveState('openmm.xml')"
        return output

    def createRestartLogic(self, total_cycles, steps_per_cycle):
        # Creates the logic to calculate, at run time, the number of cycles that need to be run
        # based on the number of steps that have already been run
        output = ""
        output += f"total_required_cycles = {total_cycles}\n"
        output += "if not is_restart:\n"
        output += "    steps_so_far = 0\n"
        output += "    numcycles = total_required_cycles\n"
        output += "else:\n"
        output += "    steps_so_far = step\n"
        output += f"    cycles_so_far = steps_so_far / {steps_per_cycle}\n"
        output += "    numcycles = int(total_required_cycles - cycles_so_far)\n"
        return output

    def createLoopWithReporting(
        self,
        name,
        steps_per_cycle,
        report_interval,
        timestep,
        inflex_point,
    ):
        """Creates the loop in which simulations are run, stopping each cycle
        to report the potential energies required for MBAR analysis.

        Parameters
        ----------
        cycles : int
            Number of cycles to run the simulation for.
        steps_per_cycle : int
            Number of steps to run the simulation for in each cycle.
        report_interval : int (in ps)
            Interval at which to report the potential energies.
        timestep : float (in ps)
            Timestep used in the simulation.
        steps : int
            Total number of steps that have been performed so far (Default 0).
        inflex_point : int
            The index at which the protocol changes direction. Potentials only need to be calculated for each half of the protocol.
        """
        output = ""
        output += "# Reporting for MBAR:\n"
        # round master lambda to 4 d.p. to avoid floating point errors
        output += f"master_lambda_list = {[round(i,4) for i in self.protocol.get_lambda_values()]}\n"
        output += f"master_lambda = master_lambda_list[window_index]\n"

        output += "if is_restart:\n"
        output += "    try:\n"
        output += "        MBAR_df = pd.read_csv(f'energies_{master_lambda}.csv')\n"
        output += "        energies = MBAR_df.to_dict('list')\n"
        output += "        energies = {float(k) if k.replace('.', '').isdigit() else k: v for k, v in energies.items()}\n"
        output += "    except FileNotFoundError:\n"
        output += "        raise FileNotFoundError('MBAR data not found, unable to restart')\n"
        output += "else:\n"
        output += "    energies = {}\n"
        output += "    energies['time'] = []\n"
        output += "    energies['fep-lambda'] = []\n"
        output += "    energies['temperature'] = []\n"
        output += "    for i in master_lambda_list:\n"
        output += "        energies[i] = []\n"
        output += f"\n# Run the simulation in cycles, with each cycle having {report_interval} steps.\n"
        output += "# Timestep in ps\n"
        output += f"timestep = {timestep}\n"
        output += f"inflex_point = {inflex_point}\n"
        output += f"for x in range(0, numcycles):\n"
        output += f"    simulation.step({steps_per_cycle})\n"
        output += f"    steps_so_far += {steps_per_cycle}\n"
        output += "    time = steps_so_far * timestep\n"
        output += "    energies['time'].append(time)\n"
        output += "    energies['fep-lambda'].append(master_lambda)\n"
        output += "    energies['temperature'].append(integrator.getTemperature().value_in_unit(kelvin))\n"
        output += "    #now loop over all simulate lambda values, set the values in the context, and calculate potential energy\n"
        output += "    # do the first half of master lambda if direction == 1\n"
        output += "    if direction == 1:\n"
        output += (
            "        for ind, lam in enumerate(master_lambda_list[:inflex_point]):\n"
        )
        output += "            for key in atm_constants.keys():\n"
        output += "                if key in ['Alpha','Uh','W0']:\n"
        output += "                    simulation.context.setParameter(key, atm_constants[key][ind].value_in_unit(kilojoules_per_mole))\n"
        output += "                else:\n"
        output += "                    simulation.context.setParameter(key, atm_constants[key][ind])\n"
        output += "            state = simulation.context.getState(getEnergy=True)\n"
        output += "            energies[lam].append(state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))\n"
        output += "    #fill the rest of the dictionary with NaNs\n"
        output += "        for lam in master_lambda_list[inflex_point:]:\n"
        output += "            energies[lam].append(float('nan'))\n"
        output += "    # do the second half of master lambda if direction == -1\n"
        output += "    else:\n"
        output += "    #fill the first half of the dictionary with NaNs\n"
        output += "        for lam in master_lambda_list[:inflex_point]:\n"
        output += "            energies[lam].append(float('nan'))\n"
        output += (
            "        for ind, lam in enumerate(master_lambda_list[inflex_point:]):\n"
        )
        output += "            for key in atm_constants.keys():\n"
        output += "                if key in ['Alpha','Uh','W0']:\n"
        output += "                    simulation.context.setParameter(key, atm_constants[key][ind+inflex_point].value_in_unit(kilojoules_per_mole))\n"
        output += "                else:\n"
        output += "                    simulation.context.setParameter(key, atm_constants[key][ind+inflex_point])\n"
        output += "            state = simulation.context.getState(getEnergy=True)\n"
        output += "            energies[lam].append(state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))\n"
        output += (
            "    #Now reset lambda-dependent values back to their original state\n"
        )
        output += "    simulation.context.setParameter('Lambda1',lambda1)\n"
        output += "    simulation.context.setParameter('Lambda2',lambda2)\n"
        output += "    simulation.context.setParameter('Alpha',alpha)\n"
        output += "    simulation.context.setParameter('Uh',uh)\n"
        output += "    simulation.context.setParameter('W0',w0)\n"
        output += "    simulation.context.setParameter('Direction',direction)\n"
        output += f"    simulation.saveState('{name}.xml')\n"
        output += "    #now dump data to a csv file\n"
        output += f"    df = pd.DataFrame(energies)\n"
        output += "    df.set_index(['time', 'fep-lambda'], inplace=True)\n"
        output += "    df.to_csv(f'energies_{master_lambda}.csv')\n"
        output += "#Dump final data to csv file\n"
        output += "df = pd.DataFrame(energies)\n"
        output += "df.set_index(['time', 'fep-lambda'], inplace=True)\n"
        output += "df.to_csv(f'energies_{master_lambda}.csv')\n"
        return output

    def createSoftcorePertELoop(
        self,
        name,
        steps_per_cycle,
        report_interval,
        timestep,
    ):
        """Recreation of Gallachio lab analysis - currently uses {cycles} to define sampling frequency"""
        output = ""
        output += f"\n# Run the simulation in cycles, with each cycle having {report_interval} steps.\n"
        output += "# Timestep in ps\n"
        output += f"timestep = {timestep}\n"
        output += "\n"
        output += "#Create dictionary for storing results in the same manner as the Gallachio lab code\n"
        # Logic for restarting simulations
        output += "#Reporting for UWHAM:\n"
        output += "if is_restart:\n"
        # first UWHAM
        output += "    try:\n"
        output += f"        UWHAM_df = pd.read_csv('{name}.csv')\n"
        output += "        result = UWHAM_df.to_dict('list')\n"
        output += "    except FileNotFoundError:\n"
        output += "        raise FileNotFoundError('UWHAM data not found, unable to restart')\n"
        output += "else:\n"
        output += "    result = {}\n"
        output += "    result['window'] = []\n"
        output += "    result['temperature'] = []\n"
        output += "    result['direction'] = []\n"
        output += "    result['lambda1'] = []\n"
        output += "    result['lambda2'] = []\n"
        output += "    result['alpha'] = []\n"
        output += "    result['uh'] = []\n"
        output += "    result['w0'] = []\n"
        output += "    result['pot_en'] = []\n"
        output += "    result['pert_en'] = []\n"
        output += "    result['metad_offset'] = []\n"

        output += f"for x in range(0, numcycles):\n"
        output += f"    simulation.step({steps_per_cycle})\n"
        output += (
            "    state = simulation.context.getState(getEnergy = True, groups = -1)\n"
        )
        output += "    pot_energy = state.getPotentialEnergy()\n"
        output += "    (u1, u0, alchemicalEBias) = atm_force.getPerturbationEnergy(simulation.context)\n"
        output += "    umcore = simulation.context.getParameter(atm_force.Umax())* kilojoules_per_mole\n"
        output += "    ubcore = simulation.context.getParameter(atm_force.Ubcore())* kilojoules_per_mole\n"
        output += "    acore = simulation.context.getParameter(atm_force.Acore())\n"
        output += "    uoffset = 0.0 * kilojoules_per_mole\n"
        output += (
            "    direction = simulation.context.getParameter(atm_force.Direction())\n"
        )
        output += "    if direction > 0:\n"
        output += (
            "        pert_e = softCorePertE(u1-(u0+uoffset), umcore, ubcore, acore)\n"
        )
        output += "    else:\n"
        output += (
            "        pert_e = softCorePertE(u0-(u1+uoffset), umcore, ubcore, acore)\n"
        )
        output += "    result['window'].append(window_index)\n"
        output += "    result['temperature'].append(integrator.getTemperature().value_in_unit(kelvin))\n"
        output += "    result['direction'].append(direction)\n"
        output += "    result['lambda1'].append(lambda1)\n"
        output += "    result['lambda2'].append(lambda2)\n"
        output += (
            "    result['alpha'].append(alpha.value_in_unit(kilocalories_per_mole))\n"
        )
        output += "    result['uh'].append(uh.value_in_unit(kilocalories_per_mole))\n"
        output += "    result['w0'].append(w0.value_in_unit(kilocalories_per_mole))\n"
        output += "    result['pot_en'].append(pot_energy.value_in_unit(kilocalories_per_mole))\n"
        output += "    result['pert_en'].append(pert_e.value_in_unit(kilocalories_per_mole))\n"
        output += "    result['metad_offset'].append(0.0)\n"
        output += "    #save the state of the simulation\n"
        output += f"    simulation.saveState('{name}.xml')\n"
        output += "    #now dump data to a csv file\n"
        output += "    df = pd.DataFrame(result)\n"
        output += "    df.set_index('window', inplace=True)\n"
        output += f"    df.to_csv(f'{name}.csv')\n"

        output += "#now convert the final dictionary to a pandas dataframe\n"
        output += "df = pd.DataFrame(result)\n"
        output += "df.set_index('window', inplace= True)\n"
        output += f"df.to_csv('{name}.csv')\n"
        return output

    def createReportingBoth(
        self,
        name,
        steps_per_cycle,
        timestep,
        inflex_point,
    ):
        output = ""
        output += "# Timestep in ps\n"
        output += f"timestep = {timestep}\n"
        output += "\n"
        # Logic for restarting simulations
        output += "#Reporting for UWHAM:\n"
        output += "if is_restart:\n"
        # first UWHAM
        output += "    try:\n"
        output += f"        UWHAM_df = pd.read_csv('{name}.csv')\n"
        output += "        result = UWHAM_df.to_dict('list')\n"
        output += "    except FileNotFoundError:\n"
        output += "        raise FileNotFoundError('UWHAM data not found, unable to restart')\n"
        output += "else:\n"
        output += "    result = {}\n"
        output += "    result['window'] = []\n"
        output += "    result['temperature'] = []\n"
        output += "    result['direction'] = []\n"
        output += "    result['lambda1'] = []\n"
        output += "    result['lambda2'] = []\n"
        output += "    result['alpha'] = []\n"
        output += "    result['uh'] = []\n"
        output += "    result['w0'] = []\n"
        output += "    result['pot_en'] = []\n"
        output += "    result['pert_en'] = []\n"
        output += "    result['metad_offset'] = []\n"

        output += "# Reporting for MBAR:\n"
        # round master lambda to 4 d.p. to avoid floating point errors
        output += f"master_lambda_list = {[round(i,4) for i in self.protocol.get_lambda_values()]}\n"
        output += f"master_lambda = master_lambda_list[window_index]\n"
        output += "if is_restart:\n"
        output += "    try:\n"
        output += "        MBAR_df = pd.read_csv(f'energies_{master_lambda}.csv')\n"
        output += "        energies = MBAR_df.to_dict('list')\n"
        output += "        energies = {float(k) if k.replace('.', '').isdigit() else k: v for k, v in energies.items()}\n"
        output += "    except FileNotFoundError:\n"
        output += "        raise FileNotFoundError('MBAR data not found, unable to restart')\n"
        output += "else:\n"
        output += "    energies = {}\n"
        output += "    energies['time'] = []\n"
        output += "    energies['fep-lambda'] = []\n"
        output += "    energies['temperature'] = []\n"
        output += "    for i in master_lambda_list:\n"
        output += "        energies[i] = []\n"

        output += f"inflex_point = {inflex_point}\n"

        output += "# Now run the simulation.\n"
        output += f"for x in range(0, numcycles):\n"
        output += f"    simulation.step({steps_per_cycle})\n"
        output += f"    steps_so_far += {steps_per_cycle}\n"
        output += "    time = steps_so_far * timestep\n"
        output += (
            "    state = simulation.context.getState(getEnergy = True, groups = -1)\n"
        )
        output += "    pot_energy = state.getPotentialEnergy()\n"
        output += "    (u1, u0, alchemicalEBias) = atm_force.getPerturbationEnergy(simulation.context)\n"
        output += "    umcore = simulation.context.getParameter(atm_force.Umax())* kilojoules_per_mole\n"
        output += "    ubcore = simulation.context.getParameter(atm_force.Ubcore())* kilojoules_per_mole\n"
        output += "    acore = simulation.context.getParameter(atm_force.Acore())\n"
        output += "    uoffset = 0.0 * kilojoules_per_mole\n"
        output += (
            "    direction = simulation.context.getParameter(atm_force.Direction())\n"
        )
        output += "    if direction > 0:\n"
        output += (
            "        pert_e = softCorePertE(u1-(u0+uoffset), umcore, ubcore, acore)\n"
        )
        output += "    else:\n"
        output += (
            "        pert_e = softCorePertE(u0-(u1+uoffset), umcore, ubcore, acore)\n"
        )
        output += "    result['window'].append(window_index)\n"
        output += "    result['temperature'].append(integrator.getTemperature().value_in_unit(kelvin))\n"
        output += "    result['direction'].append(direction)\n"
        output += "    result['lambda1'].append(lambda1)\n"
        output += "    result['lambda2'].append(lambda2)\n"
        output += (
            "    result['alpha'].append(alpha.value_in_unit(kilocalories_per_mole))\n"
        )
        output += "    result['uh'].append(uh.value_in_unit(kilocalories_per_mole))\n"
        output += "    result['w0'].append(w0.value_in_unit(kilocalories_per_mole))\n"
        output += "    result['pot_en'].append(pot_energy.value_in_unit(kilocalories_per_mole))\n"
        output += "    result['pert_en'].append(pert_e.value_in_unit(kilocalories_per_mole))\n"
        output += "    result['metad_offset'].append(0.0)\n"
        output += "    energies['time'].append(time)\n"
        output += "    energies['fep-lambda'].append(master_lambda)\n"
        output += "    energies['temperature'].append(integrator.getTemperature().value_in_unit(kelvin))\n"
        output += "    #now loop over all simulate lambda values, set the values in the context, and calculate potential energy\n"
        output += "    # do the first half of master lambda if direction == 1\n"
        output += "    if direction == 1:\n"
        output += (
            "        for ind, lam in enumerate(master_lambda_list[:inflex_point]):\n"
        )
        output += "            for key in atm_constants.keys():\n"
        output += "                if key in ['Alpha','Uh','W0']:\n"
        output += "                    simulation.context.setParameter(key, atm_constants[key][ind].value_in_unit(kilojoules_per_mole))\n"
        output += "                else:\n"
        output += "                    simulation.context.setParameter(key, atm_constants[key][ind])\n"
        output += "            state = simulation.context.getState(getEnergy=True)\n"
        output += "            energies[lam].append(state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))\n"
        output += "    #fill the rest of the dictionary with NaNs\n"
        output += "        for lam in master_lambda_list[inflex_point:]:\n"
        output += "            energies[lam].append(float('nan'))\n"
        output += "    # do the second half of master lambda if direction == -1\n"
        output += "    else:\n"
        output += "    #fill the first half of the dictionary with NaNs\n"
        output += "        for lam in master_lambda_list[:inflex_point]:\n"
        output += "            energies[lam].append(float('nan'))\n"
        output += (
            "        for ind, lam in enumerate(master_lambda_list[inflex_point:]):\n"
        )
        output += "            for key in atm_constants.keys():\n"
        output += "                if key in ['Alpha','Uh','W0']:\n"
        output += "                    simulation.context.setParameter(key, atm_constants[key][ind+inflex_point].value_in_unit(kilojoules_per_mole))\n"
        output += "                else:\n"
        output += "                    simulation.context.setParameter(key, atm_constants[key][ind+inflex_point])\n"
        output += "            state = simulation.context.getState(getEnergy=True)\n"
        output += "            energies[lam].append(state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))\n"
        output += (
            "    #Now reset lambda-dependent values back to their original state\n"
        )
        output += "    simulation.context.setParameter('Lambda1',lambda1)\n"
        output += "    simulation.context.setParameter('Lambda2',lambda2)\n"
        output += "    simulation.context.setParameter('Alpha',alpha)\n"
        output += "    simulation.context.setParameter('Uh',uh)\n"
        output += "    simulation.context.setParameter('W0',w0)\n"
        output += "    simulation.context.setParameter('Direction',direction)\n"
        output += "    #save the state of the simulation\n"
        output += f"    simulation.saveState('{name}.xml')\n"
        output += "    #now dump UWHAM data to a csv file\n"
        output += "    df = pd.DataFrame(result)\n"
        output += "    df.set_index('window', inplace=True)\n"
        output += f"    df.to_csv(f'{name}.csv')\n"
        output += "    #now dump MBAR data to a csv file\n"
        output += "    df = pd.DataFrame(energies)\n"
        output += "    df.set_index(['time', 'fep-lambda'], inplace=True)\n"
        output += "    df.to_csv(f'energies_{master_lambda}.csv')\n"

        output += "#now convert the UWHAM dictionary to a pandas dataframe\n"
        output += "df = pd.DataFrame(result)\n"
        output += "df.set_index('window', inplace= True)\n"
        output += f"df.to_csv('{name}.csv')\n"

        output += "# same for MBAR\n"
        output += "df = pd.DataFrame(energies)\n"
        output += "df.set_index(['time', 'fep-lambda'], inplace=True)\n"
        output += "df.to_csv(f'energies_{master_lambda}.csv')\n"

        return output

    def createSinglePointTest(
        self,
        inflex_point,
        name,
        atm_force_group=None,
        position_restraint_force_group=None,
        alignment_force_groups=None,
        com_force_group=None,
    ):
        """Create a single point test for the ATM force"""
        output = ""
        output += "# Create the dictionary which will hold the energies\n"
        output += f"master_lambda_list = {[round(i,4) for i in self.protocol.get_lambda_values()]}\n"
        output += "energies = {}\n"
        output += f"for i in master_lambda_list[:{inflex_point}]:\n"
        output += "    energies[i] = []\n"
        # First we can check the potential of forces that are not lambda-dependent, this will only work if the ATMforce is in its own group

        if (
            (position_restraint_force_group is not None)
            and (alignment_force_groups is not None)
            and (com_force_group is not None)
        ):
            output += "non_lambda_forces = {}\n"
            output += f"pos_state = simulation.context.getState(getEnergy=True, groups={{{position_restraint_force_group}}})\n"
            output += "non_lambda_forces['position_restraint'] = pos_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)\n"
            output += f"alignment_force_groups = {alignment_force_groups}\n"
            output += "for counter,group in enumerate(alignment_force_groups):\n"
            output += "    if counter == 0:\n"
            output += "        name='distance'\n"
            output += "    elif counter == 1:\n"
            output += "        name='angle'\n"
            output += "    elif counter == 2:\n"
            output += "        name='dihedral'\n"
            output += "    align_state = simulation.context.getState(getEnergy=True, groups={group})\n"
            output += "    non_lambda_forces[name] = align_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)\n"
            output += f"com_state = simulation.context.getState(getEnergy=True, groups={{{com_force_group}}})\n"
            output += "non_lambda_forces['com'] = com_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)\n"
            # now save as a dataframe
            output += "df = pd.DataFrame(non_lambda_forces,index=[0])\n"
            output += "df.to_csv(f'non_lambda_forces.csv')\n"

        output += "#now loop over all simulate lambda values, set the values in the context, and calculate potential energy\n"
        output += f"for ind, lam in enumerate(master_lambda_list[:{inflex_point}]):\n"
        output += "    for key in atm_constants.keys():\n"
        output += "        if key in ['Alpha','Uh','W0']:\n"
        output += "            simulation.context.setParameter(key, atm_constants[key][ind].value_in_unit(kilojoules_per_mole))\n"
        output += "        else:\n"
        output += "            simulation.context.setParameter(key, atm_constants[key][ind])\n"
        if atm_force_group is None:
            output += "    state = simulation.context.getState(getEnergy=True)\n"
        else:
            group_placeholder = f"{atm_force_group}"
            output += f"    state = simulation.context.getState(getEnergy=True, groups={{{group_placeholder}}})\n"
        output += "    energies[lam].append(state.getPotentialEnergy().value_in_unit(kilojoules_per_mole))\n"
        output += (
            "    #Now reset lambda-dependent values back to their original state\n"
        )
        output += "    simulation.context.setParameter('Lambda1',lambda1)\n"
        output += "    simulation.context.setParameter('Lambda2',lambda2)\n"
        output += "    simulation.context.setParameter('Alpha',alpha)\n"
        output += "    simulation.context.setParameter('Uh',uh)\n"
        output += "    simulation.context.setParameter('W0',w0)\n"
        output += "    simulation.context.setParameter('Direction',direction)\n"
        output += "#now convert the dictionary to a pandas dataframe, with both time and fep-lambda as index columns\n"
        output += "df = pd.DataFrame(energies)\n"
        output += "df.to_csv(f'energies_singlepoint.csv')\n"
        output += "simulation.step(1)\n"
        output += f"simulation.saveState('{name}.xml')\n"
        return output
