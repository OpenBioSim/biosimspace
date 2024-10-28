import pytest

import numpy as np

from functools import partial

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align import decouple
from BioSimSpace.Sandpit.Exscientia._Exceptions import AnalysisError
from BioSimSpace.Sandpit.Exscientia.FreeEnergy import (
    RestraintSearch,
    Restraint,
)
from BioSimSpace.Sandpit.Exscientia.Trajectory import Trajectory
from BioSimSpace.Sandpit.Exscientia.Units.Angle import degree
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
from BioSimSpace.Sandpit.Exscientia.Units.Length import nanometer

from tests.Sandpit.Exscientia.conftest import (
    url,
    has_gromacs,
    has_mdanalysis,
    has_mdrestraints_generator,
)


@pytest.fixture(scope="module")
def setup_system():
    "Setup the system for the tests"
    ligand = BSS.IO.readMolecules(
        [
            f"{url}/ligand04.prm7.bz2",
            f"{url}/ligand04.rst7.bz2",
        ]
    ).getMolecule(0)
    decouple_system = decouple(ligand).toSystem()
    protocol = BSS.Protocol.Production(runtime=BSS.Types.Time(8, "FEMTOSECOND"))
    return ligand, decouple_system, protocol


# Make sure GROMACS is installed.
has_gromacs = BSS._gmx_exe is not None


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_run_Gromacs(setup_system):
    "Test if the normal run works on Gromacs"
    ligand, decouple_system, protocol = setup_system
    restraint_search = RestraintSearch(
        decouple_system, protocol=protocol, engine="GROMACS"
    )
    restraint_search.start()
    restraint_search.wait()
    assert not restraint_search._process.isError()


def test_run_Somd(setup_system):
    "Test if the normal run works with SOMD"
    ligand, decouple_system, protocol = setup_system
    restraint_search = RestraintSearch(
        decouple_system, protocol=protocol, engine="SOMD"
    )
    restraint_search.start()
    restraint_search.wait()
    assert not restraint_search._process.isError()


@pytest.mark.skipif(
    (
        has_mdrestraints_generator is False
        or has_gromacs is False
        or has_mdanalysis is False
    ),
    reason="Requires MDRestraintsGenerator, MDAnalysis and Gromacs to be installed.",
)
class TestMDRestraintsGenerator_analysis:
    @staticmethod
    @pytest.fixture(scope="class")
    def restraint_search(tmp_path_factory):
        outdir = tmp_path_factory.mktemp("out")
        system = BSS.IO.readMolecules(
            [
                f"{url}/crd.gro.bz2",
                f"{url}/complex.top.bz2",
            ]
        )
        ligand = system.getMolecule(1)
        decoupled_ligand = decouple(ligand)
        protein = system.getMolecule(0)
        new_system = (protein + decoupled_ligand).toSystem()

        protocol = BSS.Protocol.Production()
        restraint_search = BSS.FreeEnergy.RestraintSearch(
            new_system,
            protocol=protocol,
            engine="GROMACS",
            work_dir=str(outdir),
        )
        traj, top = BSS.IO.expand(url, ["traj.xtc", "complex.tpr"], ".bz2")
        restraint_search._process.getTrajectory = lambda: Trajectory(
            trajectory=traj, topology=top
        )
        restraint = restraint_search.analyse(
            method="MDRestraintsGenerator",
            restraint_type="Boresch",
            block=False,
        )
        return restraint, outdir

    def test_sanity(self, restraint_search):
        restraint, outdir = restraint_search
        assert isinstance(restraint, Restraint)

    def test_plots(self, restraint_search):
        """Test if all the plots has been generated correctly."""
        restraint, outdir = restraint_search
        assert (outdir / "bond_1.png").is_file()
        assert (outdir / "angle_1.png").is_file()
        assert (outdir / "angle_2.png").is_file()
        assert (outdir / "dihedral_1.png").is_file()
        assert (outdir / "dihedral_2.png").is_file()
        assert (outdir / "dihedral_3.png").is_file()

    def test_dG_off(self, restraint_search):
        """Test if the restraint generated is valid."""
        restraint, outdir = restraint_search
        assert (outdir / "dG_off.dat").is_file()
        dG = np.loadtxt(outdir / "dG_off.dat")
        assert isinstance(dG, np.ndarray)

    def test_top(self, restraint_search):
        """Test if the restraint generated has the same energy."""
        restraint, outdir = restraint_search
        assert (outdir / "BoreschRestraint.top").is_file()
        with open(outdir / "BoreschRestraint.top", "r") as f:
            assert "intermolecular_interactions" in f.read()

@pytest.mark.skipif(
    (has_gromacs is False or has_mdanalysis is False),
    reason="Requires MDAnalysis and Gromacs to be installed.",
)
class TestBSS_analysis:
    """Test selection of restraints using the inbuilt BSS method."""

    @staticmethod
    def _restraint_search_general(tmp_path_factory, ligand_idx=1):
        outdir = tmp_path_factory.mktemp("out")
        system = BSS.IO.readMolecules(
            [
                f"{url}/crd.gro.bz2",
                f"{url}/complex.top.bz2",
            ]
        )
        protein = system.getMolecule(0)
        if ligand_idx == 0:
            ligand = protein.copy()
        else:
            ligand = system.getMolecule(ligand_idx)
        decoupled_ligand = decouple(ligand)
        new_system = (protein + decoupled_ligand).toSystem()

        protocol = BSS.Protocol.Production()
        restraint_search = BSS.FreeEnergy.RestraintSearch(
            new_system, protocol=protocol, engine="GROMACS", work_dir=str(outdir)
        )
        traj, top = BSS.IO.expand(url, ["traj.xtc", "complex.tpr"], ".bz2")
        restraint_search._process.getTrajectory = lambda: Trajectory(
            trajectory=traj, topology=top
        )
        return restraint_search, outdir

    @pytest.fixture(scope="class")
    def _restraint_search(self, tmp_path_factory):
        return partial(self._restraint_search_general, tmp_path_factory, ligand_idx=1)()

    @pytest.fixture(scope="class")
    def _restraint_search_protein(self, tmp_path_factory):
        # Select the protein as the ligand.
        return partial(self._restraint_search_general, tmp_path_factory,ligand_idx=0)()

    @staticmethod
    @pytest.fixture(scope="class")
    def boresch_restraint(_restraint_search):
        restraint_search, outdir = _restraint_search
        restraint = restraint_search.analyse(
            method="BSS", restraint_type="Boresch", block=False
        )
        return restraint, outdir

    @staticmethod
    @pytest.fixture(scope="class")
    def boresch_restraint_k20(_restraint_search):
        restraint_search, outdir = _restraint_search
        restraint = restraint_search.analyse(
            method="BSS",
            restraint_type="Boresch",
            block=False,
            force_constant=20 * kcal_per_mol / angstrom**2,
        )
        return restraint, outdir

    @staticmethod
    @pytest.fixture(scope="class")
    def multiple_distance_restraint(_restraint_search):
        restraint_search, outdir = _restraint_search
        restraint = restraint_search.analyse(
            method="BSS", restraint_type="multiple_distance", block=False
        )
        return restraint, outdir

    def test_sanity(self, boresch_restraint):
        restraint, _ = boresch_restraint
        assert isinstance(restraint, Restraint)

    def test_plots_boresch(self, boresch_restraint):
        """Test if all the plots have been generated correctly"""
        restraint, outdir = boresch_restraint
        assert (outdir / "restraint_idx0_dof_time.png").is_file()
        assert (outdir / "restraint_idx0_dof_hist.png").is_file()

    def test_dG_off_boresch(self, boresch_restraint):
        """Test if the restraint generated has the same energy"""
        restraint, _ = boresch_restraint
        assert np.isclose(-9.2133, restraint.correction.value(), atol=0.01)

    def test_bond_boresch(self, boresch_restraint):
        restraint, _ = boresch_restraint
        equilibrium_values_r0 = (
            restraint._restraint_dict["equilibrium_values"]["r0"] / nanometer
        )
        assert np.isclose(0.4987, equilibrium_values_r0, atol=0.001)

    def test_angles_boresch(self, boresch_restraint):
        restraint, _ = boresch_restraint
        equilibrium_values_thetaA0 = (
            restraint._restraint_dict["equilibrium_values"]["thetaA0"] / degree
        )
        assert np.isclose(119.3375, equilibrium_values_thetaA0, atol=0.001)
        equilibrium_values_thetaB0 = (
            restraint._restraint_dict["equilibrium_values"]["thetaB0"] / degree
        )
        assert np.isclose(100.8382, equilibrium_values_thetaB0, atol=0.001)

    def test_dihedrals_boresch(self, boresch_restraint):
        restraint, _ = boresch_restraint
        equilibrium_values_phiA0 = (
            restraint._restraint_dict["equilibrium_values"]["phiA0"] / degree
        )
        assert np.isclose(29.1947, equilibrium_values_phiA0, atol=0.001)
        equilibrium_values_phiB0 = (
            restraint._restraint_dict["equilibrium_values"]["phiB0"] / degree
        )
        assert np.isclose(-136.9009, equilibrium_values_phiB0, atol=0.001)
        equilibrium_values_phiC0 = (
            restraint._restraint_dict["equilibrium_values"]["phiC0"] / degree
        )
        assert np.isclose(-102.3908, equilibrium_values_phiC0, atol=0.001)

    def test_index_boresch(self, boresch_restraint):
        restraint, _ = boresch_restraint
        idxs = {
            k: restraint._restraint_dict["anchor_points"][k].index()
            for k in restraint._restraint_dict["anchor_points"]
        }
        assert idxs == {"r1": 1560, "r2": 1558, "r3": 1562, "l1": 8, "l2": 9, "l3": 13}

    def test_force_constant_boresch(self, boresch_restraint_k20):
        restraint, _ = boresch_restraint_k20
        for force_constant in restraint._restraint_dict["force_constants"].values():
            assert np.isclose(20, force_constant.value(), atol=0.01)

    def test_analysis_failure_boresch(self, _restraint_search):
        restraint_search, _ = _restraint_search
        with pytest.raises(AnalysisError):
            restraint = restraint_search.analyse(
                method="BSS",
                restraint_type="Boresch",
                block=False,
                cutoff=0.1 * angstrom,
            )

    def test_non_overlapping_anchor_selections_forbidden(self, _restraint_search_protein):
        """
        Ensure that anchor atoms cannot be in the same molecule
        (the protein has been selected as the ligand, so we are
        guaranteed that anchor atoms will be in the same molecule).
        """
        restraint_search, _ = _restraint_search_protein
        with pytest.raises(AnalysisError):
            restraint = restraint_search.analyse(
                method="BSS",
                restraint_type="Boresch",
                block=False,
            )

    def test_plots_mdr(self, multiple_distance_restraint):
        """Test if all the plots have been generated correctly"""
        restraint, outdir = multiple_distance_restraint
        assert (outdir / "distance_restraints_time.png").is_file()
        assert (outdir / "distance_restraints_hist.png").is_file()

    def test_dG_off_mdr(self, multiple_distance_restraint):
        """Test if the restraint generated has the same energy"""
        restraint, _ = multiple_distance_restraint
        assert np.isclose(-0.0790, restraint.correction.value(), atol=0.001)

    def test_dict_mdr(self, multiple_distance_restraint):
        """Check that the MDR restraint dictionary is correct"""
        restraint, _ = multiple_distance_restraint
        restr_dict = restraint._restraint_dict
        # Check the permanent distance restraint - this should always be the same.
        assert restr_dict["permanent_distance_restraint"]["r1"].index() == 1539
        assert restr_dict["permanent_distance_restraint"]["l1"].index() == 10
        assert restr_dict["permanent_distance_restraint"]["r0"].unit() == "ANGSTROM"
        # Check close
        assert restr_dict["permanent_distance_restraint"][
            "r0"
        ].value() == pytest.approx(8.9019, abs=1e-4)
        assert restr_dict["permanent_distance_restraint"]["kr"].dimensions() == (
            1,
            0,
            -2,
            0,
            0,
            -1,
            0,
        )
        assert restr_dict["permanent_distance_restraint"]["kr"].value() == 40.0
        assert restr_dict["permanent_distance_restraint"]["r_fb"].unit() == "ANGSTROM"
        assert restr_dict["permanent_distance_restraint"][
            "r_fb"
        ].value() == pytest.approx(0.5756, abs=1e-4)
        # Check the distance restraints - the parameters will vary slightly due to random initialisation
        # of the clustering algorithm during restraint searching.
        assert len(restr_dict["distance_restraints"]) == 9
        idxs_receptor = [
            pair["r1"].index()
            for pair in restraint._restraint_dict["distance_restraints"]
        ]
        idxs_ligand = [
            pair["l1"].index()
            for pair in restraint._restraint_dict["distance_restraints"]
        ]
        # The search algorithm is not deterministic, so the receptor anchor points
        # are not always the same. The ligand anchor points are always the same as
        # all heavy atoms are used.
        expected_idxs_ligand = [13, 9, 16, 12, 11, 14, 15, 8, 17].sort()
        assert len(idxs_receptor) == len(idxs_ligand)
        assert idxs_ligand.sort() == expected_idxs_ligand

    def test_parameters_mdr(self, multiple_distance_restraint):
        """
        Check the parameters selected for the permanent distance restraint.
        Other restraints are not tested as the multiple distance restraint
        selection algorithm is not deterministic.
        """
        restraint, _ = multiple_distance_restraint
        permanent_restraint = restraint._restraint_dict["permanent_distance_restraint"]
        assert np.isclose(8.9019, permanent_restraint["r0"] / angstrom, atol=0.001)
        assert np.isclose(
            40,
            permanent_restraint["kr"] / (kcal_per_mol / angstrom**2),
        )
        assert np.isclose(0.5756, permanent_restraint["r_fb"] / angstrom, atol=0.001)

    def test_analysis_failure_mdr(self, _restraint_search):
        restraint_search, _ = _restraint_search
        with pytest.raises(AnalysisError):
            restraint = restraint_search.analyse(
                method="BSS",
                restraint_type="multiple_distance",
                block=False,
                cutoff=0.1 * angstrom,
            )
