import bz2
from math import exp
import pandas as pd
import pathlib
import pytest
import numpy as np
import os
import shutil
import time

from tests.Sandpit.Exscientia.conftest import (
    has_alchemlyb,
    has_alchemlyb_parquet,
    has_alchemtest,
    has_gromacs,
    url,
)
from tests.conftest import root_fp

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Protocol import FreeEnergyEquilibration
from BioSimSpace.Sandpit.Exscientia.Align._decouple import decouple
from BioSimSpace.Sandpit.Exscientia import Types as _Types
from BioSimSpace.Sandpit.Exscientia.FreeEnergy import Restraint
from BioSimSpace.Sandpit.Exscientia.Units.Angle import radian, degree
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
from BioSimSpace.Sandpit.Exscientia.Units.Temperature import kelvin


@pytest.mark.skipif(
    has_alchemtest is False, reason="Requires alchemtest and alchemlyb to be installed."
)
class TestRelativeAnalysis:
    @staticmethod
    @pytest.fixture(scope="class")
    def gmx_data(tmp_path_factory):
        from alchemtest.gmx import load_ABFE

        outdir = tmp_path_factory.mktemp("gromacs")
        for leg in ["complex", "ligand"]:
            for index, filename in enumerate(load_ABFE().data[leg]):
                dir = outdir / leg / f"lambda_{index}"
                dir.mkdir(parents=True)
                eng = dir / "gromacs.xvg"
                eng.symlink_to(filename)
        return outdir

    @staticmethod
    @pytest.fixture(scope="class")
    def amber_data(tmp_path_factory):
        from alchemtest.amber import load_bace_example

        outdir = tmp_path_factory.mktemp("amber")
        for leg, stage in [("complex", "decharge"), ("solvated", "vdw")]:
            for index, filename in enumerate(load_bace_example().data[leg][stage]):
                dir = outdir / f"{leg}_{stage}" / f"lambda_{index}"
                dir.mkdir(parents=True)
                with open(dir / "amber.out", "w") as f:
                    with bz2.open(filename, "rt") as bz_file:
                        f.write(bz_file.read())
        return outdir

    @staticmethod
    @pytest.fixture(scope="class")
    def gmx_complex(gmx_data):
        outdir = gmx_data
        complex, _ = BSS.FreeEnergy.AlchemicalFreeEnergy.analyse(
            work_dir=str(outdir / "complex"),
            temperature=310 * BSS.Units.Temperature.kelvin,
        )
        return complex

    @staticmethod
    @pytest.fixture(scope="class")
    def gmx_ligand(gmx_data):
        outdir = gmx_data
        ligand, _ = BSS.FreeEnergy.AlchemicalFreeEnergy.analyse(
            work_dir=str(outdir / "ligand"),
            temperature=310 * BSS.Units.Temperature.kelvin,
        )
        return ligand

    @staticmethod
    @pytest.fixture(scope="class")
    def amber_complex_decharge(amber_data):
        outdir = amber_data
        complex, _ = BSS.FreeEnergy.AlchemicalFreeEnergy.analyse(
            work_dir=str(outdir / "complex_decharge"),
            temperature=298 * BSS.Units.Temperature.kelvin,
        )
        return complex

    @staticmethod
    @pytest.fixture(scope="class")
    def amber_solvated_vdw(amber_data):
        outdir = amber_data
        complex, _ = BSS.FreeEnergy.AlchemicalFreeEnergy.analyse(
            work_dir=str(outdir / "solvated_vdw"),
            temperature=298 * BSS.Units.Temperature.kelvin,
        )
        return complex

    @pytest.mark.parametrize(
        "fixture,length,energy",
        [
            ("gmx_ligand", 20, 7.654472744451637),
            ("gmx_complex", 30, 21.819752),
            ("amber_complex_decharge", 5, -5.25352),
            ("amber_solvated_vdw", 12, 2.261816),
        ],
    )
    def test_pmf(self, fixture, length, energy, request):
        pmf = request.getfixturevalue(fixture)
        assert len(pmf) == length
        assert len(pmf[0]) == 3
        np.testing.assert_allclose(
            pmf[-1][1] / BSS.Units.Energy.kcal_per_mol, energy, atol=0.1
        )

    def test_difference(self, gmx_complex, gmx_ligand):
        dG, error = BSS.FreeEnergy.AlchemicalFreeEnergy.difference(
            gmx_complex, gmx_ligand
        )
        np.testing.assert_allclose(
            dG / BSS.Units.Energy.kcal_per_mol, 14.216101, atol=0.1
        )


@pytest.mark.skipif(
    has_alchemlyb_parquet is False, reason="Requires alchemlyb > 2.1.0."
)
class TestAnalysePARQUET:
    @staticmethod
    @pytest.fixture(scope="class")
    def data(tmp_path_factory):
        outdir = tmp_path_factory.mktemp("out")
        shutil.copytree(
            f"{root_fp}/Sandpit/Exscientia/input/parquet", outdir / "parquet"
        )
        return str(outdir / "parquet")

    def test_analyse(self, data):
        result = BSS.FreeEnergy.AlchemicalFreeEnergy.analyse(
            data, temperature=300 * BSS.Units.Temperature.kelvin, estimator="MBAR"
        )
        assert np.isclose(
            result[0][-1][-1] / BSS.Units.Energy.kcal_per_mol, 20.87341050030068, atol=1
        )


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
@pytest.mark.skipif(
    has_alchemlyb is False, reason="Requires alchemlyb to be installed."
)
class Test_gmx_ABFE:
    @staticmethod
    @pytest.fixture(scope="class")
    def freenrg():
        m = BSS.IO.readMolecules(
            [f"{root_fp}/input/ala.top", f"{root_fp}/input/ala.crd"]
        )[0]
        decouple_m = decouple(m)
        solvated = BSS.Solvent.tip3p(
            molecule=decouple_m, box=3 * [3 * BSS.Units.Length.nanometer]
        )
        protocol = FreeEnergyEquilibration(
            lam=pd.Series(data={"coul": 0.5, "vdw": 0.0}),
            lam_vals=pd.DataFrame(
                data={
                    "coul": [0.0, 0.5, 1.0, 1.0, 1.0],
                    "vdw": [0.0, 0.0, 0.0, 0.5, 1.0],
                }
            ),
            runtime=_Types.Time(0, "nanoseconds"),
        )
        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            solvated,
            protocol,
            engine="GROMACS",
        )
        freenrg.run()
        freenrg.wait()
        return freenrg

    def test_file_exist(self, freenrg):
        """Test if all the files are there."""
        path = pathlib.Path(freenrg.workDir())
        for i in range(5):
            assert (path / f"lambda_{i}" / "gromacs.xvg").is_file()

    def test_lambda(self, freenrg):
        """Test if the xvg files contain the correct lambda."""
        path = pathlib.Path(freenrg.workDir())
        for i, (coul, vdw) in enumerate(
            zip([0.0, 0.5, 1.0, 1.0, 1.0], [0.0, 0.0, 0.0, 0.5, 1.0])
        ):
            from alchemlyb.parsing.gmx import extract_u_nk

            u_nk = extract_u_nk(path / f"lambda_{i}" / "gromacs.xvg", 300)
            assert u_nk.index.names == ["time", "coul-lambda", "vdw-lambda"]
            assert np.isclose(u_nk.index.values[0][1], coul)
            assert np.isclose(u_nk.index.values[0][2], vdw)


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
class Test_Somd_ABFE:
    @staticmethod
    @pytest.fixture(scope="class")
    def freenrg():
        # Just use a single ligand with anchors based all on atoms in the ligand, for simplicity
        m = BSS.IO.readMolecules(
            [f"{url}/ligand01.prm7.bz2", f"{url}/ligand01.rst7.bz2"]
        ).getMolecule(0)

        # Assign atoms for restraint
        atom_1 = m.getAtoms()[0]
        atom_2 = m.getAtoms()[1]
        atom_3 = m.getAtoms()[2]
        atom_4 = m.getAtoms()[3]
        atom_5 = m.getAtoms()[4]
        atom_6 = m.getAtoms()[5]

        decoupled_m = decouple(m)
        solvated = BSS.Solvent.tip3p(
            molecule=decoupled_m, box=3 * [3 * BSS.Units.Length.nanometer]
        )

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
            solvated, restraint_dict, 298 * kelvin, restraint_type="Boresch"
        )
        protocol = BSS.Protocol.FreeEnergy(
            lam_vals=[0.0, 0.5, 1.0], runtime=0.0001 * BSS.Units.Time.nanosecond
        )
        freenrg = BSS.FreeEnergy.AlchemicalFreeEnergy(
            solvated, protocol, engine="SOMD", restraint=restraint
        )

        freenrg.run()
        freenrg.wait()
        # Sleep to allow files to be written
        time.sleep(25)
        return freenrg

    def test_files_exist(self, freenrg):
        """
        Test if the files have been created. Note that e.g. gradients.dat
        are not created until later in the simulation, so their presence is
        not tested for.
        """
        path = pathlib.Path(freenrg.workDir())
        for lam in ["0.0000", "0.5000", "1.0000"]:
            assert (path / f"lambda_{lam}" / "simfile.dat").is_file()
            assert (path / f"lambda_{lam}" / "SYSTEM.s3").is_file()
            assert (path / f"lambda_{lam}" / "somd.cfg").is_file()
            assert (path / f"lambda_{lam}" / "somd.rst7").is_file()
            assert (path / f"lambda_{lam}" / "somd.prm7").is_file()
            assert (path / f"lambda_{lam}" / "somd.err").is_file()
            assert (path / f"lambda_{lam}" / "somd.out").is_file()
            assert (path / f"lambda_{lam}" / "somd.pert").is_file()

    def test_correct_pert_file(self, freenrg):
        """Check that pert file is correct."""
        path = pathlib.Path(freenrg.workDir()) / "lambda_0.0000"
        with open(os.path.join(path, "somd.pert"), "rt") as f:
            lines = f.readlines()

            for i, line in enumerate(lines):
                # Check that the end-state properties are correct.
                if "final_type" in line:
                    assert "final_type     du" in line
                if "final_LJ" in line:
                    assert "final_LJ       0.00000 0.00000" in line
                if "final_charge" in line:
                    assert "final_charge   0.00000" in line
                # Check that the initial state properties are correct for the first and last atoms.
                if line == "        name           C1\n":
                    assert "initial_type   c" in lines[i + 1]
                    assert "initial_LJ     3.39967 0.08600" in lines[i + 3]
                    assert "initial_charge 0.67120" in lines[i + 5]
                if "name           O3" in line:
                    assert "initial_type   o" in lines[i + 1]
                    assert "initial_LJ     2.95992 0.21000" in lines[i + 3]
                    assert "initial_charge -0.52110" in lines[i + 5]

    def test_correct_conf_file(self, freenrg):
        """Check that lambda data is correct in somd.cfg"""
        path = pathlib.Path(freenrg.workDir())
        for lam in ["0.0000", "0.5000", "1.0000"]:
            with open(os.path.join(path, f"lambda_{lam}", "somd.cfg"), "rt") as f:
                lines = f.readlines()
                assert "lambda array = 0.0, 0.5, 1.0\n" in lines
                assert f"lambda_val = {float(lam):.1f}\n" in lines
                assert f"perturbed residue number = 1\n" in lines
