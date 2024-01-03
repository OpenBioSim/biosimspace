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
from ..Protocol import Protocol as _Protocol


class _AToMUtils:
    # Internal class for creating openmm forces within an AToM process.
    def __init__(self, protocol):
        # Check for proper typing
        if not isinstance(protocol, _Protocol.AToM):
            raise TypeError("Protocol must be a BioSimSpace.Protocol object")
        self.protocol = protocol
        self.data = self.protocol.getData()

    def getAlignmentConstants(self):
        self.alignment_k_distance = self.protocol.getAlignKfSep()
        self.alignment_k_theta = self.protocol.getAlignKTheta()
        self.alignment_k_psi = self.protocol.getAlignKPsi()

    def findAbsoluteCoreIndices(self):
        import numpy as np

        self.lig1_first_atomnum = self.data["first_lig1_atom_index"]
        self.lig1_rigid_atoms = list(
            np.add(self.lig1_first_atomnum, self.data["lig1_rigid_core"])
        )
        self.lig2_first_atomnum = self.data["first_lig2_atom_index"]
        self.lig2_rigid_atoms = list(
            np.add(self.lig1_first_atomnum, self.data["lig2_rigid_core"])
        )

    def findDisplacement(self):
        if not isinstance(self.data["displacement"], list):
            raise TypeError("Displacement must be a list")
        elif not all(isinstance(x, float) for x in self.data["displacement"]):
            raise TypeError("Displacement must be a list of floats")
        self.displacement = self.data["displacement"]

    def createAlignmentForce(self):
        # This force is the same in every lambda window
        self.get_alignment_constants()
        self.find_absolute_core_indices()
        self.find_displacement()

        output = "\n\n"
        d = [round(x, 3) for x in self.displacement]
        output += "displacment = {}\n".format(d)
        output += "k_distance = {}\n".format(self.alignment_k_distance)
        output += "k_theta = {}\n".format(self.alignment_k_theta)
        output += "k_psi = {}\n".format(self.alignment_k_psi)
        output += "idxs_a = {}\n".format(self.lig1_rigid_atoms)
        output += "idxs_b = {}\n".format(self.lig2_rigid_atoms)
        output += "\n\n"

        output += 'distance_energy_fn = "0.5 * k * ((x1 - x2 - dx)^2 + (y1 - y2 - dy)^2 + (z1 - z2 - dz)^2);"\n'
        output += (
            "distance_force = openmm.CustomCompoundBondForce(2, distance_energy_fn)\n"
        )
        output += "distance_force.addPerBondParameter('k')\n"
        output += "distance_force.addPerBondParameter('dx')\n"
        output += "distance_force.addPerBondParameter('dy')\n"
        output += "distance_force.addPerBondParameter('dz')\n"

        output += """distance_parameters = [
        k_distance*openmm.unit(
            openmm.unit.kilojoules_per_mole / openmm.unit.angstrom**2
        ),
        displacment[0]*openmm.unit.angstrom,
        displacment[1]*openmm.unit.angstrom,
        displacment[2]*openmm.unit.angstrom,
        ]\n"""

        output += "distance_force.addBond((idxs_b[0], idxs_a[0]), distance_parameters)"
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
        output += "angle_force = openmm.CustomCompoundBondForce(4, angle_energy_fn)\n"
        output += 'angle_force.addPerBondParameter("k")\n'
        output += """angle_force.addBond(
        (idxs_b[0], idxs_b[1], idxs_a[0], idxs_a[1]),
        [k_theta*openmm.unit.kilojoules_per_mole],
        )\n"""
        output += "system.addForce(angle_force)\n\n"

        output += """dihedral_energy_fn = (
        "0.5 * k * (1 - cos_phi);"
        ""
        "cos_phi = (v_x * w_x + v_y * w_y + v_z * w_z) / (norm_v * norm_w);"
        ""
        "norm_v = sqrt(v_x^2 + v_y^2 + v_z^2);"
        "v_x = dx_31 - dot_31 * dx_21 / norm_21;"
        "v_y = dy_31 - dot_31 * dy_21 / norm_21;"
        "v_z = dz_31 - dot_31 * dz_21 / norm_21;"
        ""
        "dot_31 = (dx_31 * dx_21 + dy_31 * dy_21 + dz_31 * dz_21) / norm_21;"
        "dx_31 = x3 - x1; dy_31 = y3 - y1; dz_31 = z3 - z1;"
        ""
        "norm_w = sqrt(w_x^2 + w_y^2 + w_z^2);"
        "w_x = dx_54 - dot_54 * dx_21 / norm_21;"
        "w_y = dy_54 - dot_54 * dy_21 / norm_21;"
        "w_z = dz_54 - dot_54 * dz_21 / norm_21;"
        ""
        "dot_54 =(dx_54 * dx_21 + dy_54 * dy_21 + dz_54 * dz_21) / norm_21;"
        "dx_54 = x5 - x4; dy_54 = y5 - y4; dz_54 = z5 - z4;"
        ""
        "norm_21 = sqrt(dx_21^2 + dy_21^2 + dz_21^2);"
        "dx_21 = x2 - x1; dy_21 = y2 - y1; dz_21 = z2 - z1;"
        )\n"""

        output += (
            "dihedral_force = openmm.CustomCompoundBondForce(5, dihedral_energy_fn)\n"
        )
        output += 'dihedral_force.addPerBondParameter("k")\n'
        output += """dihedral_force.addBond(
        (idxs_b[0], idxs_b[1], idxs_b[2], idxs_a[0], idxs_a[2]),
        [0.5 * k_psi * unit(openmm.unit.kilojoules_per_mole)],
        )\n"""
        output += """dihedral_force.addBond(
        (idxs_a[0], idxs_a[1], idxs_a[2], idxs_b[0], idxs_b[2]),
        [0.5 * k_dihedral*unit(openmm.unit.kilojoules_per_mole)],
        )\n"""
        output += "system.addForce(dihedral_force)\n\n"
        return output

    def getLigand1AtomsAsList(self):
        import numpy as np

        return list(
            np.arange(
                self.data["first_lig1_atom_index"],
                self.data["last_ligand1_atom_index"] + 1,
            )
        )

    def getLigand2AtomsAsList(self):
        import numpy as np

        return list(
            np.arange(
                self.data["first_lig2_atom_index"],
                self.data["last_ligand2_atom_index"] + 1,
            )
        )

    def createATMForce(
        self,
        lambda1,
        lambda2,
        displacement,
        alpha,
        u0,
        w0,
        direction,
        sc_Umax,
        sc_U0,
        sc_a,
    ):
        output = ""
        output += "#Parameters for ATM force\n"
        output += "lig1_atoms = {}\n".format(self.getLigand1AtomsAsList())
        output += "lig2_atoms = {}\n".format(self.getLigand2AtomsAsList())
        output += "displacement = {}\n".format(displacement)
        output += "lambda1 = {}\n".format(lambda1)
        output += "lambda2 = {}\n".format(lambda2)
        output += "alpha = {}\n".format(alpha)
        output += "u0 = {}\n".format(u0)
        output += "w0 = {}\n".format(w0)
        output += "direction = {}\n".format(direction)
        output += "sc_Umax = {}\n".format(sc_Umax)
        output += "sc_U0 = {}\n".format(sc_U0)
        output += "sc_a = {}\n".format(sc_a)

        output += "\n\n #Define ATM force\n"
        output += """atm_force = openmm.ATMForce(
        lambda1,
        lambda2,
        alpha * openmm.unit.kilojoules_per_mole,
        u0 * openmm.unit.kilojoules_per_mole,   
        w0 * openmm.unit.kilojoules_per_mole,
        sc_Umax * openmm.unit.kilojoules_per_mole,
        sc_u0 * openmm.unit.kilojoules_per_mole,
        sc_a,
        direction,
        )"""

        output += "\n\n #Add ATM force to system\n"
        output += "for _ in range(len(self._topology.atoms)):\n"
        output += "    atm_force.addParticle(openmm.Vec3(0.0,0.0,0.0))\n"
        # TODO: add offset - check convesion of a list to a Vec3
