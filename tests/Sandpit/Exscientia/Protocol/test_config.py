import pandas as pd
import pytest

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align._decouple import decouple
from BioSimSpace.Sandpit.Exscientia.Protocol import (
    ConfigFactory,
    Equilibration,
    FreeEnergy,
    FreeEnergyEquilibration,
    FreeEnergyMinimisation,
    Production,
)

from BioSimSpace.Sandpit.Exscientia.Align._decouple import decouple
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
from BioSimSpace.Sandpit.Exscientia.Units.Angle import radian, degree
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol
from BioSimSpace.Sandpit.Exscientia.Units.Temperature import kelvin
from BioSimSpace.Sandpit.Exscientia.FreeEnergy import Restraint
from BioSimSpace.Sandpit.Exscientia._SireWrappers import Molecule
from BioSimSpace.Sandpit.Exscientia._Utils import _try_import, _have_imported


from tests.Sandpit.Exscientia.conftest import (
    url,
    has_gromacs,
    has_antechamber,
    has_openff,
)


class TestAmber:
    @pytest.fixture(scope="class")
    def system(self):
        m0 = BSS.IO.readMolecules(
            [f"{url}/CAT-13c.prm7.bz2", f"{url}/CAT-13c.rst7.bz2"]
        ).getMolecule(0)
        return m0.toSystem()

    def test_NVT(self, system):
        config = ConfigFactory(system, Equilibration(pressure=None))
        res = [x.strip().strip(",") for x in config.generateAmberConfig()]
        assert "ntb=1" in res


class TestAmberRBFE:
    @pytest.fixture(scope="class")
    def system(self):
        m0 = BSS.IO.readMolecules(
            [f"{url}/CAT-13c.prm7.bz2", f"{url}/CAT-13c.rst7.bz2"]
        ).getMolecule(0)
        m1 = BSS.IO.readMolecules(
            [f"{url}/CAT-13a.prm7.bz2", f"{url}/CAT-13a.rst7.bz2"]
        ).getMolecule(0)
        atom_mapping = BSS.Align.matchAtoms(m0, m1)
        m0 = BSS.Align.rmsdAlign(m0, m1, atom_mapping)
        merged = BSS.Align.merge(m0, m1)
        return merged.toSystem()

    def test_generate_fep_masks(self, system):
        config = ConfigFactory(system, FreeEnergyMinimisation())
        res = config._generate_amber_fep_masks(0.004)
        expected_res = {
            "noshakemask": '""',
            "scmask1": '"@25-27,29-31"',
            "scmask2": '""',
            "timask1": '"@1-45"',
            "timask2": '"@46-84"',
        }
        assert res == expected_res

    def test_generate_restraint_masks(self, system):
        protocol = FreeEnergyEquilibration(restraint=[0, 24], force_constant=4.3)
        config = ConfigFactory(system, protocol)
        res = [x.strip().strip(",") for x in config.generateAmberConfig()]
        expected_res = {"ntr=1", 'restraintmask="@1,25,46"', "restraint_wt=4.3"}
        assert expected_res.issubset(res)

    @pytest.mark.parametrize(
        "protocol", [Equilibration, Production, FreeEnergyEquilibration, FreeEnergy]
    )
    def test_tau_t(self, system, protocol):
        config = ConfigFactory(system, protocol(tau_t=BSS.Types.Time(2, "picosecond")))
        res = [x.strip().strip(",") for x in config.generateAmberConfig()]
        expected_res = {"gamma_ln=0.50000"}
        assert expected_res.issubset(res)


class TestGromacsRBFE:
    @staticmethod
    @pytest.fixture(scope="class")
    def system():
        m0 = BSS.IO.readMolecules(
            [f"{url}/CAT-13a.prm7.bz2", f"{url}/CAT-13a.rst7.bz2"]
        ).getMolecule(0)
        m1 = BSS.IO.readMolecules(
            [f"{url}/CAT-13c.prm7.bz2", f"{url}/CAT-13c.rst7.bz2"]
        ).getMolecule(0)
        atom_mapping = BSS.Align.matchAtoms(m0, m1)
        m0 = BSS.Align.rmsdAlign(m0, m1, atom_mapping)
        merged = BSS.Align.merge(m0, m1)
        return merged.toSystem()

    @pytest.mark.parametrize(
        "protocol", [Equilibration, Production, FreeEnergyEquilibration, FreeEnergy]
    )
    def test_tau_t(self, system, protocol):
        config = ConfigFactory(system, protocol(tau_t=BSS.Types.Time(2, "picosecond")))
        res = config.generateGromacsConfig()
        expected_res = {"tau-t = 2.00000"}
        assert expected_res.issubset(res)

    @pytest.mark.skipif(
        has_gromacs is False, reason="Requires GROMACS to be installed."
    )
    def test_fep(self, system):
        """Test if the default config writer will write the expected thing."""
        protocol = FreeEnergyMinimisation(
            lam=0.0,
            lam_vals=None,
            min_lam=0.0,
            max_lam=1.0,
            num_lam=11,
            steps=10000,
            perturbation_type="full",
        )
        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            system,
            protocol,
            engine="GROMACS",
        )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", "r") as f:
            mdp_text = f.read()
            assert (
                "fep-lambdas          = 0.00000 0.10000 0.20000 0.30000 0.40000 0.50000 0.60000 0.70000 0.80000 0.90000 1.00000"
                in mdp_text
            )
            assert "init-lambda-state = 6" in mdp_text

    @pytest.mark.skipif(
        has_gromacs is False, reason="Requires GROMACS to be installed."
    )
    def test_fep_df(self, system):
        """Test if the default config writer configured with the pd.DataFrame
        will write the expected thing.
        """
        protocol = FreeEnergyMinimisation(
            lam=pd.Series(data={"fep": 0.0}),
            lam_vals=None,
            min_lam=pd.Series(data={"fep": 0.0}),
            max_lam=pd.Series(data={"fep": 1.0}),
            num_lam=11,
            steps=10000,
            perturbation_type="full",
        )
        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            system,
            protocol,
            engine="GROMACS",
        )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", "r") as f:
            mdp_text = f.read()
            assert (
                "fep-lambdas          = 0.00000 0.10000 0.20000 0.30000 0.40000 0.50000 0.60000 0.70000 0.80000 0.90000 1.00000"
                in mdp_text
            )
            assert "init-lambda-state = 6" in mdp_text

    @pytest.mark.skipif(
        has_gromacs is False, reason="Requires GROMACS to be installed."
    )
    def test_staged_fep_df(self, system):
        """Test if the multi-stage lambda will be written correctly."""
        protocol = FreeEnergyMinimisation(
            lam=pd.Series(data={"bonded": 0.0, "coul": 0.0, "vdw": 0.0}),
            lam_vals=pd.DataFrame(
                data={
                    "bonded": [0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0],
                    "coul": [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
                    "vdw": [0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0],
                }
            ),
            steps=10000,
            perturbation_type="full",
        )
        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            system,
            protocol,
            engine="GROMACS",
        )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", "r") as f:
            mdp_text = f.read()
            assert (
                "bonded-lambdas       = 0.00000 0.25000 0.50000 0.75000 1.00000 1.00000 1.00000"
                in mdp_text
            )
            assert (
                "coul-lambdas         = 0.00000 0.00000 0.00000 0.50000 1.00000 1.00000 1.00000"
                in mdp_text
            )
            assert (
                "vdw-lambdas          = 0.00000 0.00000 0.00000 0.00000 0.00000 0.50000 1.00000"
                in mdp_text
            )
            assert "init-lambda-state = 6" in mdp_text


@pytest.mark.skipif(
    has_antechamber is False or has_openff is False,
    reason="Requires ambertools/antechamber and openff to be installed",
)
@pytest.fixture(scope="module")
def system_and_mdr_restraint():
    # Benzene.
    m = BSS.Parameters.openff_unconstrained_2_0_0("c1ccccc1").getMolecule()
    m._sire_object = m._sire_object.edit().rename("LIG").molecule().commit()

    # Assign atoms for restraint
    atom_1 = m.getAtoms()[0]
    atom_2 = m.getAtoms()[1]
    atom_3 = m.getAtoms()[2]
    atom_4 = m.getAtoms()[3]
    atom_5 = m.getAtoms()[4]
    atom_6 = m.getAtoms()[5]

    mol = decouple(m)
    system = mol.toSystem()

    # Create random restraint dictionary
    restraint_dict = {
        "distance_restraints": [
            {
                "l1": atom_1,
                "r1": atom_2,
                "r0": 3 * angstrom,
                "kr": 10 * kcal_per_mol / angstrom**2,
                "r_fb": 1 * angstrom,
            },
            {
                "l1": atom_3,
                "r1": atom_4,
                "r0": 3 * angstrom,
                "kr": 10 * kcal_per_mol / angstrom**2,
                "r_fb": 1 * angstrom,
            },
        ],
        "permanent_distance_restraint": {
            "l1": atom_5,
            "r1": atom_6,
            "r0": 3 * angstrom,
            "kr": 10 * kcal_per_mol / angstrom**2,
            "r_fb": 1 * angstrom,
        },
    }

    restraint = Restraint(
        system, restraint_dict, 298 * kelvin, restraint_type="multiple_distance"
    )

    return system, restraint


@pytest.mark.skipif(
    has_antechamber is False or has_openff is False,
    reason="Requires ambertools/antechamber and openff to be installed",
)
class TestGromacsABFE:
    @staticmethod
    @pytest.fixture(scope="class")
    def system():
        # Benzene.
        m0 = BSS.Parameters.openff_unconstrained_2_0_0("c1ccccc1").getMolecule()
        m0._sire_object = m0._sire_object.edit().rename("LIG").molecule().commit()
        protocol = FreeEnergyMinimisation(
            lam=0.0,
            lam_vals=None,
            min_lam=0.0,
            max_lam=1.0,
            num_lam=11,
            steps=10000,
            perturbation_type="full",
        )
        return m0, protocol

    @pytest.mark.skipif(
        has_gromacs is False, reason="Requires GROMACS to be installed."
    )
    def test_decouple_vdw_q(self, system):
        m, protocol = system
        """Test the decoupling where lambda0 = vdw-q and lambda1=none."""
        mol = decouple(m)
        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            mol.toSystem(),
            protocol,
            engine="GROMACS",
        )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", "r") as f:
            mdp_text = f.read()
            assert "couple-moltype = LIG" in mdp_text
            assert "couple-lambda0 = vdw-q" in mdp_text
            assert "couple-lambda1 = none" in mdp_text
            assert "couple-intramol = yes" in mdp_text


    def test_decouple_perturbable(self, system):
        m, protocol = system
        mol = decouple(m)
        sire_mol = mol._sire_object
        c = sire_mol.cursor()
        for key in [
            "charge",
            "LJ",
            "bond",
            "angle",
            "dihedral",
            "improper",
            "forcefield",
            "intrascale",
            "mass",
            "element",
            "atomtype",
            "coordinates",
            "velocity",
            "ambertype",
        ]:
            if f"{key}1" not in c and key in c:
                c[f"{key}0"] = c[key]
                c[f"{key}1"] = c[key]

        c["is_perturbable"] = True
        sire_mol = c.commit()
        mol = Molecule(sire_mol)

        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            mol.toSystem(),
            protocol,
            engine="GROMACS",
        )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", "r") as f:
            mdp_text = f.read()
            assert "couple-moltype" not in mdp_text
            assert "couple-lambda0" not in mdp_text
            assert "couple-lambda1" not in mdp_text
            assert "couple-intramol" not in mdp_text


    @pytest.mark.skipif(
        has_gromacs is False, reason="Requires GROMACS to be installed."
    )
    def test_annihilate_vdw2q(self, system):
        """Test the annihilation where lambda0 = vdw-q and lambda1=none."""
        m, protocol = system
        mol = decouple(m, charge=(False, True), LJ=(True, False), intramol=False)
        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            mol.toSystem(),
            protocol,
            engine="GROMACS",
        )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", "r") as f:
            mdp_text = f.read()
            assert "couple-moltype = LIG" in mdp_text
            assert "couple-lambda0 = vdw" in mdp_text
            assert "couple-lambda1 = q" in mdp_text
            assert "couple-intramol = no" in mdp_text

    @pytest.mark.skipif(
        has_gromacs is False, reason="Requires GROMACS to be installed."
    )
    def test_mdr_force_constant(self, system, system_and_mdr_restraint):
        """
        Check that when the restraint type == 'release_restraint', the
        force constant of the permanent distance restraint (not affected by
        restraint-lambda) is written to the MDP file.
        """
        system, restraint = system_and_mdr_restraint
        protocol = FreeEnergy(
            perturbation_type="release_restraint",
        )

        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            system,
            protocol,
            engine="GROMACS",
            restraint=restraint,
        )
        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", "r") as f:
            mdp_text = f.read()
            assert "disre-fc = 4184.0" in mdp_text

    @pytest.mark.skipif(
        has_gromacs is False, reason="Requires GROMACS to be installed."
    )
    def test_sc_parameters(self, system):
        """Test if the soft core parameters have been written.
        The default sc-alpha is 0, which means the soft-core of the vdw is not
        turned on by default. This checks if the this value has been changed to
        0.5.
        """
        m, protocol = system
        mol = decouple(m)
        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            mol.toSystem(),
            protocol,
            engine="GROMACS",
        )

        with open(f"{freenrg._work_dir}/lambda_6/gromacs.mdp", "r") as f:
            mdp_text = f.read()
            assert "sc-alpha = 0.5" in mdp_text


@pytest.mark.skipif(
    has_antechamber is False or has_openff is False,
    reason="Requires ambertools/antechamber and openff to be installed",
)
class TestSomdABFE:
    @staticmethod
    @pytest.fixture(scope="class")
    def system_and_boresch_restraint():
        # Benzene.
        m = BSS.Parameters.openff_unconstrained_2_0_0("c1ccccc1").getMolecule()

        # Assign atoms for restraint
        atom_1 = m.getAtoms()[0]
        atom_2 = m.getAtoms()[1]
        atom_3 = m.getAtoms()[2]
        atom_4 = m.getAtoms()[3]
        atom_5 = m.getAtoms()[4]
        atom_6 = m.getAtoms()[5]

        mol = decouple(m)
        system = mol.toSystem()

        # Create random restraint dictionary
        restraint_dict = {
            "anchor_points": {
                "r1": atom_1,
                "r2": atom_2,
                "r3": atom_3,
                "l1": atom_4,
                "l2": atom_5,
                "l3": atom_6,
            },
            "equilibrium_values": {
                "r0": 5.08 * angstrom,
                "thetaA0": 64.051 * degree,
                "thetaB0": 39.618 * degree,
                "phiA0": 2.59 * radian,
                "phiB0": -1.20 * radian,
                "phiC0": 2.63 * radian,
            },
            "force_constants": {
                "kr": 10 * kcal_per_mol / angstrom**2,
                "kthetaA": 10 * kcal_per_mol / (radian * radian),
                "kthetaB": 10 * kcal_per_mol / (radian * radian),
                "kphiA": 10 * kcal_per_mol / (radian * radian),
                "kphiB": 10 * kcal_per_mol / (radian * radian),
                "kphiC": 10 * kcal_per_mol / (radian * radian),
            },
        }

        restraint = Restraint(
            system, restraint_dict, 298 * kelvin, restraint_type="Boresch"
        )

        return system, restraint

    def test_turn_on_restraint_boresch(self, system_and_boresch_restraint):
        """Test for turning on multiple distance restraints"""
        system, restraint = system_and_boresch_restraint
        protocol = FreeEnergy(perturbation_type="restraint")
        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            system, protocol, engine="SOMD", restraint=restraint
        )

        # Test .cfg file
        with open(f"{freenrg._work_dir}/lambda_0.0000/somd.cfg", "r") as f:
            cfg_text = f.read()
            assert "use boresch restraints = True" in cfg_text
            assert 'boresch restraints dictionary = {"anchor_points":{"r1":1,'
            ' "r2":2, "r3":3, "l1":4, "l2":5, "l3":6}, "equilibrium_values"'
            ':{"r0":5.08, "thetaA0":1.12, "thetaB0":0.69,"phiA0":2.59, "phiB0"'
            ':-1.20, "phiC0":2.63}, "force_constants":{"kr":5.00, "kthetaA":5.00,'
            ' "kthetaB":5.00, "kphiA":5.00, "kphiB":5.00, "kphiC":5.00}}' in cfg_text
            assert "turn on receptor-ligand restraints mode = True" in cfg_text

        # Test .pert file
        with open(f"{freenrg._work_dir}/lambda_0.0000/somd.pert", "r") as f:
            pert_text = f.read()
            # Check all atoms present
            carbons = [f"C{i}" for i in range(1, 7)]
            hydrogens = [f"H{i}" for i in range(1, 7)]
            atoms = carbons + hydrogens
            for atom in atoms:
                assert atom in pert_text
            # Check perturbations are correct
            lines = [
                "molecule LIG",
                "atom",
                "initial_type   C1",
                "final_type     C1",
                "initial_LJ     3.48065 0.08688",
                "final_LJ       3.48065 0.08688",
                "initial_charge -0.13000",
                "final_charge   -0.13000",
                "endatom",
                "endmolecule",
            ]
            for line in lines:
                assert line in pert_text

    def test_turn_on_restraint_mdr(self, system_and_mdr_restraint):
        """Test for turning on the mdr restraint"""
        system, restraint = system_and_mdr_restraint
        protocol = FreeEnergy(perturbation_type="restraint")
        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            system, protocol, engine="SOMD", restraint=restraint
        )

        # Test .cfg file
        with open(f"{freenrg._work_dir}/lambda_0.0000/somd.cfg", "r") as f:
            cfg_text = f.read()
            assert "use distance restraints = True" in cfg_text
            assert (
                "distance restraints dictionary = {(1, 0): (3.0, 5.0, 1.0), (3, 2): "
                "(3.0, 5.0, 1.0), (5, 4): (3.0, 5.0, 1.0)}"
            ) in cfg_text
            assert "turn on receptor-ligand restraints mode = True" in cfg_text

        # Test .pert file
        with open(f"{freenrg._work_dir}/lambda_0.0000/somd.pert", "r") as f:
            pert_text = f.read()
            # Check all atoms present
            carbons = [f"C{i}" for i in range(1, 7)]
            hydrogens = [f"H{i}" for i in range(1, 7)]
            atoms = carbons + hydrogens
            for atom in atoms:
                assert atom in pert_text
            # Check perturbations are correct
            lines = [
                "molecule LIG",
                "atom",
                "initial_type   C1",
                "final_type     C1",
                "initial_LJ     3.48065 0.08688",
                "final_LJ       3.48065 0.08688",
                "initial_charge -0.13000",
                "final_charge   -0.13000",
                "endatom",
                "endmolecule",
            ]
            for line in lines:
                assert line in pert_text

    def test_release_restraint_mdr(self, system_and_mdr_restraint):
        """Test for releasing the non-permanent mdr restraints"""
        system, restraint = system_and_mdr_restraint
        protocol = FreeEnergy(perturbation_type="release_restraint")
        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            system, protocol, engine="SOMD", restraint=restraint
        )

        # Test .cfg file
        with open(f"{freenrg._work_dir}/lambda_0.0000/somd.cfg", "r") as f:
            cfg_text = f.read()
            assert "use distance restraints = True" in cfg_text
            assert "use permanent distance restraints = True" in cfg_text
            assert "turn on receptor-ligand restraints mode = True" in cfg_text
            assert (
                "distance restraints dictionary = {(1, 0): (3.0, 5.0, 1.0), "
                "(3, 2): (3.0, 5.0, 1.0)}"
            ) in cfg_text
            assert (
                "permanent distance restraints dictionary = {(5, 4): (3.0, 5.0, 1.0)}"
                in cfg_text
            )

        # Test .pert file
        with open(f"{freenrg._work_dir}/lambda_0.0000/somd.pert", "r") as f:
            pert_text = f.read()
            # Don't check atoms present as they are given randomised names
            # and the atom types are all du
            # Check perturbations are correct
            lines = [
                "molecule LIG",
                "atom",
                "initial_type   du",
                "final_type     du",
                "initial_LJ     0.00000 0.00000",
                "final_LJ       0.00000 0.00000",
                "initial_charge 0.00000",
                "final_charge   0.00000",
                "endatom",
                "endmolecule",
            ]
            for line in lines:
                assert line in pert_text

    def test_discharge(self, system_and_boresch_restraint):
        """Test for discharging the ligand"""
        system, restraint = system_and_boresch_restraint
        protocol = FreeEnergy(perturbation_type="discharge_soft")
        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            system, protocol, engine="SOMD", restraint=restraint
        )

        # Test .cfg file
        with open(f"{freenrg._work_dir}/lambda_0.0000/somd.cfg", "r") as f:
            cfg_text = f.read()
            assert "use boresch restraints = True" in cfg_text
            assert 'boresch restraints dictionary = {"anchor_points":{"r1":1,'
            ' "r2":2, "r3":3, "l1":4, "l2":5, "l3":6}, "equilibrium_values"'
            ':{"r0":5.08, "thetaA0":1.12, "thetaB0":0.69,"phiA0":2.59, "phiB0"'
            ':-1.20, "phiC0":2.63}, "force_constants":{"kr":5.00, "kthetaA":5.00,'
            ' "kthetaB":5.00, "kphiA":5.00, "kphiB":5.00, "kphiC":5.00}}' in cfg_text
            assert "turn on receptor-ligand restraints mode = True" not in cfg_text

        # Test .pert file
        with open(f"{freenrg._work_dir}/lambda_0.0000/somd.pert", "r") as f:
            pert_text = f.read()
            # Check all atoms present
            carbons = [f"C{i}" for i in range(1, 7)]
            hydrogens = [f"H{i}" for i in range(1, 7)]
            atoms = carbons + hydrogens
            for atom in atoms:
                assert atom in pert_text
            # Check perturbations are correct
            lines = [
                "molecule LIG",
                "atom",
                "initial_type   C1",
                "final_type     C1",
                "initial_LJ     3.48065 0.08688",
                "final_LJ       3.48065 0.08688",
                "initial_charge -0.13000",
                "final_charge   -0.00000",
                "endatom",
                "endmolecule",
            ]
            for line in lines:
                assert line in pert_text

    def test_vanish(self, system_and_boresch_restraint):
        """Test for vanishing the ligand"""
        system, restraint = system_and_boresch_restraint
        protocol = FreeEnergy(perturbation_type="vanish_soft")
        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            system, protocol, engine="SOMD", restraint=restraint
        )

        # Test .cfg file
        with open(f"{freenrg._work_dir}/lambda_0.0000/somd.cfg", "r") as f:
            cfg_text = f.read()
            assert "use boresch restraints = True" in cfg_text
            assert 'boresch restraints dictionary = {"anchor_points":{"r1":1,'
            ' "r2":2, "r3":3, "l1":4, "l2":5, "l3":6}, "equilibrium_values"'
            ':{"r0":5.08, "thetaA0":1.12, "thetaB0":0.69,"phiA0":2.59, "phiB0"'
            ':-1.20, "phiC0":2.63}, "force_constants":{"kr":5.00, "kthetaA":5.00,'
            ' "kthetaB":5.00, "kphiA":5.00, "kphiB":5.00, "kphiC":5.00}}' in cfg_text
            assert "turn on receptor-ligand restraints mode = True" not in cfg_text

        # Test .pert file
        with open(f"{freenrg._work_dir}/lambda_0.0000/somd.pert", "r") as f:
            pert_text = f.read()
            # Check all atoms present
            carbons = [f"C{i}" for i in range(1, 7)]
            hydrogens = [f"H{i}" for i in range(1, 7)]
            atoms = carbons + hydrogens
            for atom in atoms:
                assert atom in pert_text
            # Check perturbations are correct
            lines = [
                "molecule LIG",
                "atom",
                "initial_type   C1",
                "final_type     du",
                "initial_LJ     3.48065 0.08688",
                "final_LJ       0.00000 0.00000",
                "initial_charge -0.00000",
                "final_charge   -0.00000",
                "endatom",
                "endmolecule",
            ]
            for line in lines:
                assert line in pert_text

    def test_discharge_and_vanish(self, system_and_boresch_restraint):
        """Test for simultaneously discharging and vanishing the ligand"""
        system, restraint = system_and_boresch_restraint
        protocol = FreeEnergy(perturbation_type="full")
        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            system, protocol, engine="SOMD", restraint=restraint
        )

        # Test .cfg file
        with open(f"{freenrg._work_dir}/lambda_0.0000/somd.cfg", "r") as f:
            cfg_text = f.read()
            assert "use boresch restraints = True" in cfg_text
            assert 'boresch restraints dictionary = {"anchor_points":{"r1":1,'
            ' "r2":2, "r3":3, "l1":4, "l2":5, "l3":6}, "equilibrium_values"'
            ':{"r0":5.08, "thetaA0":1.12, "thetaB0":0.69,"phiA0":2.59, "phiB0"'
            ':-1.20, "phiC0":2.63}, "force_constants":{"kr":5.00, "kthetaA":5.00,'
            ' "kthetaB":5.00, "kphiA":5.00, "kphiB":5.00, "kphiC":5.00}}' in cfg_text
            assert "turn on receptor-ligand restraints mode = True" not in cfg_text

        # Test .pert file
        with open(f"{freenrg._work_dir}/lambda_0.0000/somd.pert", "r") as f:
            pert_text = f.read()
            # Check all atoms present
            carbons = [f"C{i}" for i in range(1, 7)]
            hydrogens = [f"H{i}" for i in range(1, 7)]
            atoms = carbons + hydrogens
            for atom in atoms:
                assert atom in pert_text
            # Check perturbations are correct
            lines = [
                "molecule LIG",
                "atom",
                "initial_type   C1",
                "final_type     du",
                "initial_LJ     3.48065 0.08688",
                "final_LJ       0.00000 0.00000",
                "initial_charge -0.13000",
                "final_charge   0.00000",
                "endatom",
                "endmolecule",
            ]
            for line in lines:
                assert line in pert_text


class TestAmberASFE:
    @staticmethod
    @pytest.fixture(scope="class")
    def system():
        mol = BSS.IO.readMolecules(
            [f"{url}/CAT-13c.prm7.bz2", f"{url}/CAT-13c.rst7.bz2"]
        ).getMolecule(0)
        mol = decouple(mol)
        return mol.toSystem()

    def test_generate_fep_masks(self, system):
        config = ConfigFactory(system, FreeEnergyMinimisation())
        res = config._generate_amber_fep_masks(0.001)
        expected_res = {
            "noshakemask": '"@1-45"',
            "scmask1": '"@1-45"',
            "scmask2": '""',
            "timask1": '"@1-45"',
            "timask2": '""',
        }
        for key in expected_res:
            assert expected_res[key] == res[key]
