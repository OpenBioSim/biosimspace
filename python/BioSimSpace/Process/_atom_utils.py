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
            self.direction = self.protocol.getDirections()[index]
        elif isinstance(
            self.protocol, (_Protocol.AToMEquilibration, _Protocol.AToMAnnealing)
        ):
            self.lambda1 = self.protocol.getLambda1()
            self.lambda2 = self.protocol.getLambda2()
            self.alpha = self.protocol.getAlpha().value()
            self.uh = self.protocol.getUh().value()
            self.w0 = self.protocol.getW0().value()
            self.direction = self.protocol.getDirection()

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
        output += "displacement = {}\n".format(self.displacement)
        output += "#BioSimSpace output is in angstrom, divide by 10 to convert to the expected units of nm\n"
        output += "displacement = [i/10.0 for i in displacement]\n"
        return output

    def createAlignmentForce(self):
        # This force is the same in every lambda window
        self.getAlignmentConstants()
        self.findAbsoluteCoreIndices()

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
    ):
        """
        Create a string which can be added directly to an openmm script to add an ATM force to the system.

        Parameters
        ----------
        index : int
            Index of current window - used to set window-dependent variables.
        """
        self.findDisplacement()
        self.getATMForceConstants(index)
        output = ""
        output += "#Parameters for ATM force in  original units\n"
        output += "lig1_atoms = {}\n".format(self.lig1_atoms)
        output += "lig2_atoms = {}\n".format(self.lig2_atoms)
        output += "lambda1 = {}\n".format(self.lambda1)
        output += "lambda2 = {}\n".format(self.lambda2)
        output += "alpha = {} * kilocalories_per_mole\n".format(self.alpha)
        output += "uh = {} * kilocalories_per_mole\n".format(self.uh)
        output += "w0 = {} * kilocalories_per_mole\n".format(self.w0)
        output += "direction = {}\n".format(self.direction)
        output += "sc_Umax = {} * kilocalories_per_mole\n".format(self.SCUmax)
        output += "sc_U0 = {} * kilocalories_per_mole\n".format(self.SCU0)
        output += "sc_a = {}\n".format(self.SCa)

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
        output += "system.addForce(atm_force)\n"
        return output

    def createCOMRestraint(self):
        """
        Create a string containing the CM-CM restriants for two groups of atoms.
        In most cases these will be some combination of protein and ligand atoms.
        Constants for the force are set in the protocol."""
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
        output += 'expr = "(kfcm/2)*step(d12-tolcm)*(d12-tolcm)^2 "\n'
        output += 'expr += " ; d12 = sqrt((x1 - offx - x2)^2 + (y1 - offy - y2)^2 + (z1 - offz - z2)^2 ) ; "\n'
        output += "force_CMCM = CustomCentroidBondForce(2, expr)\n"
        output += "force_CMCM.addPerBondParameter('kfcm')\n"
        output += "force_CMCM.addPerBondParameter('tolcm')\n"
        output += "force_CMCM.addPerBondParameter('offx')\n"
        output += "force_CMCM.addPerBondParameter('offy')\n"
        output += "force_CMCM.addPerBondParameter('offz')\n"
        output += "system.addForce(force_CMCM)\n"
        output += "force_CMCM.addGroup(protein_com)\n"
        output += "force_CMCM.addGroup(lig1_com)\n"
        output += "force_CMCM.addGroup(lig2_com)\n"

        output += """parameters_free = (
        kfcm.value_in_unit(kilojoules_per_mole / nanometer**2),
        tolcm.value_in_unit(nanometer),
        displacement[0] * nanometer,
        displacement[1] * nanometer,
        displacement[2] * nanometer,
        )\n"""

        output += """parameters_bound = (
        kfcm.value_in_unit(kilojoules_per_mole / nanometer**2),
        tolcm.value_in_unit(nanometer),
        0.0 * nanometer,
        0.0 * nanometer,
        0.0 * nanometer,
        )\n"""

        output += "force_CMCM.addBond((0,1), parameters_bound)\n"
        output += "force_CMCM.addBond((0,2), parameters_free)\n"
        output += "#End of CM-CM force\n\n"

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
