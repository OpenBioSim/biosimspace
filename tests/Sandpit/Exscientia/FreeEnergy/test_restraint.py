import pytest

import numpy as np

from sire.legacy.Units import angstrom3 as _Sire_angstrom3
from sire.legacy.Units import k_boltz as _k_boltz
from sire.legacy.Units import meter3 as _Sire_meter3
from sire.legacy.Units import mole as _Sire_mole

import BioSimSpace.Sandpit.Exscientia as BSS
from BioSimSpace.Sandpit.Exscientia.Align import decouple
from BioSimSpace.Sandpit.Exscientia.FreeEnergy import Restraint
from BioSimSpace.Sandpit.Exscientia.Units.Length import angstrom
from BioSimSpace.Sandpit.Exscientia.Units.Angle import radian, degree
from BioSimSpace.Sandpit.Exscientia.Units.Energy import kcal_per_mol
from BioSimSpace.Sandpit.Exscientia.Units.Temperature import kelvin

# Store the tutorial URL.
url = BSS.tutorialUrl()

################### Test Borech Restraint ################


@pytest.fixture(scope="session")
def boresch_restraint_component():
    """Generate a the components required to create a restraint."""
    ligand = BSS.IO.readMolecules(
        [f"{url}/ligand01.prm7.bz2", f"{url}/ligand01.rst7.bz2"]
    ).getMolecule(0)
    decoupled_ligand = decouple(ligand)

    protein = BSS.IO.readMolecules(
        [f"{url}/1jr5.crd.bz2", f"{url}/1jr5.top.bz2"]
    ).getMolecule(0)

    system = (protein + decoupled_ligand).toSystem()

    # Assign three atoms from the ligand
    ligand_2 = decoupled_ligand.getAtoms()[2]
    ligand_1 = decoupled_ligand.getAtoms()[1]
    ligand_0 = decoupled_ligand.getAtoms()[0]

    # And three atoms from the protein
    protein_0 = protein.getAtoms()[0]
    protein_1 = protein.getAtoms()[1]
    protein_2 = protein.getAtoms()[2]

    restraint_dict = {
        "anchor_points": {
            "r1": protein_0,
            "r2": protein_1,
            "r3": protein_2,
            "l1": ligand_0,
            "l2": ligand_1,
            "l3": ligand_2,
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

    return system, restraint_dict


@pytest.fixture(scope="session")
def boresch_restraint(boresch_restraint_component):
    """Generate the Boresch restraint object."""
    system, restraint_dict = boresch_restraint_component

    restraint = Restraint(
        system, restraint_dict, 300 * kelvin, restraint_type="Boresch"
    )
    return restraint


def test_sanity_boresch(boresch_restraint):
    """Sanity check."""
    assert isinstance(boresch_restraint, Restraint)


def test_numerical_correction_boresch(boresch_restraint):
    dG = boresch_restraint.getCorrection(method="numerical") / kcal_per_mol
    assert np.isclose(-7.2, dG, atol=0.1)


def test_analytical_correction_boresch(boresch_restraint):
    dG = boresch_restraint.getCorrection(method="analytical") / kcal_per_mol
    assert np.isclose(-7.2, dG, atol=0.1)
    assert isinstance(boresch_restraint, Restraint)


test_force_constants_boresch = [
    ({"kr": 0}, ValueError),
    ({"kthetaA": 0}, ValueError),
    ({"kthetaB": 0}, ValueError),
    (
        {
            "kthetaA": 0,
            "kphiA": 0,
            "kphiB": 0,
        },
        None,
    ),
]


@pytest.mark.parametrize("force_constants, expected", test_force_constants_boresch)
def test_input_force_constants_boresch(
    boresch_restraint_component, force_constants, expected
):
    print(force_constants)
    system, restraint_dict = boresch_restraint_component
    dict_copy = restraint_dict.copy()
    force_constants_copy = restraint_dict["force_constants"].copy()
    force_constants_copy.update(force_constants)
    dict_copy["force_constants"] = force_constants_copy
    if expected is None:
        Restraint(system, dict_copy, 300 * kelvin, restraint_type="Boresch")
    else:
        with pytest.raises(expected):
            Restraint(system, dict_copy, 300 * kelvin, restraint_type="Boresch")


class TestGromacsOutputBoresch:
    @staticmethod
    @pytest.fixture(scope="class")
    def Topology(boresch_restraint):
        return boresch_restraint.toString(engine="Gromacs").split("\n")

    def test_sanity(self, Topology):
        """Sanity check."""
        assert "intermolecular_interactions" in Topology[0]

    def test_bond(self, Topology):
        ai, aj, type, bA, kA, bB, kB = Topology[3].split()
        assert ai == "1"
        assert aj == "1496"
        assert bA == "0.508"
        assert bB == "0.508"
        assert kB == "4184.00"

    def test_angle(self, Topology):
        ai, aj, ak, type, thA, kA, thB, kB = Topology[6].split()
        assert ai == "2"
        assert aj == "1"
        assert ak == "1496"
        assert thA == "64.051"
        assert thB == "64.051"
        assert kB == "41.84"
        ai, aj, ak, type, thA, kA, thB, kB = Topology[7].split()
        assert ai == "1"
        assert aj == "1496"
        assert ak == "1497"

    def test_dihedral(self, Topology):
        ai, aj, ak, al, type, phiA, kA, phiB, kB = Topology[10].split()
        assert ai == "3"
        assert aj == "2"
        assert ak == "1"
        assert al == "1496"
        assert type == "2"
        assert phiA == "148.396"
        assert kA == "0.00"
        assert phiB == "148.396"
        assert kB == "41.84"
        ai, aj, ak, al, type, phiA, kA, phiB, kB = Topology[11].split()
        assert ai == "2"
        assert aj == "1"
        assert ak == "1496"
        assert al == "1497"
        ai, aj, ak, al, type, phiA, kA, phiB, kB = Topology[12].split()
        assert ai == "1"
        assert aj == "1496"
        assert ak == "1497"
        assert al == "1498"


class TestGromacsOutputBoreschRestraintLambda(TestGromacsOutputBoresch):
    @staticmethod
    @pytest.fixture(scope="class")
    def Topology(boresch_restraint):
        return boresch_restraint.toString(
            engine="Gromacs", restraint_lambda=True
        ).split("\n")

    def test_dihedral(self, Topology):
        assert "dihedral_restraints" in Topology[8]
        ai, aj, ak, al, type, phiA, dphiA, kA, phiB, dphiB, kB = Topology[10].split()
        assert ai == "3"
        assert aj == "2"
        assert ak == "1"
        assert al == "1496"
        assert type == "1"
        assert phiA == "148.396"
        assert dphiA == "0.00"
        assert dphiB == "0.00"
        assert kA == "0.00"
        assert phiB == "148.396"
        assert kB == "41.84"


class TestSomdOutputBoresch:
    @staticmethod
    @pytest.fixture(scope="class")
    def getRestraintSomd(boresch_restraint):
        boresch_str = boresch_restraint.toString(engine="SOMD").split("=")[1].strip()
        boresch_dict = eval(boresch_str)
        return boresch_dict

    def test_sanity(self, getRestraintSomd):
        "Sanity check"
        boresch_dict = getRestraintSomd
        assert type(boresch_dict) == dict

    def test_indices(self, getRestraintSomd):
        anchor_points = getRestraintSomd["anchor_points"]
        assert anchor_points["r1"] == 0
        assert anchor_points["r2"] == 1
        assert anchor_points["r3"] == 2
        assert anchor_points["l1"] == 1495
        assert anchor_points["l2"] == 1496
        assert anchor_points["l3"] == 1497

    def test_equil_vals(self, getRestraintSomd):
        equil_vals = getRestraintSomd["equilibrium_values"]
        assert equil_vals["r0"] == 5.08
        assert equil_vals["thetaA0"] == 1.12
        assert equil_vals["thetaB0"] == 0.69
        assert equil_vals["phiA0"] == 2.59
        assert equil_vals["phiB0"] == -1.20
        assert equil_vals["phiC0"] == 2.63


################### Test Multiple Distance Restraint ###################


@pytest.fixture(scope="session")
def mdr_restraint_component():
    """
    Generate a the components required to create a
    multiple distance restraints restraint object.
    """
    ligand = BSS.IO.readMolecules(
        [f"{url}/ligand01.prm7.bz2", f"{url}/ligand01.rst7.bz2"]
    ).getMolecule(0)
    decoupled_ligand = decouple(ligand)

    protein = BSS.IO.readMolecules(
        [f"{url}/1jr5.crd.bz2", f"{url}/1jr5.top.bz2"]
    ).getMolecule(0)

    system = (protein + decoupled_ligand).toSystem()

    # Create three distance restraints
    restraint_dict = {
        "distance_restraints": [
            {
                "l1": decoupled_ligand.getAtoms()[0],
                "r1": protein.getAtoms()[0],
                "r0": 3 * angstrom,
                "kr": 10 * kcal_per_mol / angstrom**2,
                "r_fb": 1 * angstrom,
            },
            {
                "l1": decoupled_ligand.getAtoms()[1],
                "r1": protein.getAtoms()[1],
                "r0": 3 * angstrom,
                "kr": 10 * kcal_per_mol / angstrom**2,
                "r_fb": 1 * angstrom,
            },
        ],
        "permanent_distance_restraint": {
            "l1": decoupled_ligand.getAtoms()[2],
            "r1": protein.getAtoms()[2],
            "r0": 3 * angstrom,
            "kr": 10 * kcal_per_mol / angstrom**2,
            "r_fb": 1 * angstrom,
        },
    }

    return system, restraint_dict


@pytest.fixture(scope="session")
def mdr_restraint(mdr_restraint_component):
    """Generate the multiple distance restraints restraint object."""
    system, restraint_dict = mdr_restraint_component

    restraint = Restraint(
        system, restraint_dict, 300 * kelvin, restraint_type="multiple_distance"
    )
    return restraint


@pytest.fixture(scope="session")
def mdr_restraint_component_fb_r0():
    """
    Generate a the components required to create a
    multiple distance restraints restraint object with
    the values for r0 and r_fb set to 0 for the permanent restraint.
    """
    ligand = BSS.IO.readMolecules(
        [f"{url}/ligand01.prm7.bz2", f"{url}/ligand01.rst7.bz2"]
    ).getMolecule(0)
    decoupled_ligand = decouple(ligand)

    protein = BSS.IO.readMolecules(
        [f"{url}/1jr5.crd.bz2", f"{url}/1jr5.top.bz2"]
    ).getMolecule(0)

    system = (protein + decoupled_ligand).toSystem()

    # Create three distance restraints
    restraint_dict = {
        "distance_restraints": [
            {
                "l1": decoupled_ligand.getAtoms()[0],
                "r1": protein.getAtoms()[0],
                "r0": 3 * angstrom,
                "kr": 10 * kcal_per_mol / angstrom**2,
                "r_fb": 1 * angstrom,
            },
            {
                "l1": decoupled_ligand.getAtoms()[1],
                "r1": protein.getAtoms()[1],
                "r0": 3 * angstrom,
                "kr": 10 * kcal_per_mol / angstrom**2,
                "r_fb": 1 * angstrom,
            },
        ],
        "permanent_distance_restraint": {
            "l1": decoupled_ligand.getAtoms()[2],
            "r1": protein.getAtoms()[2],
            "r0": 0 * angstrom,
            "kr": 10 * kcal_per_mol / angstrom**2,
            "r_fb": 0 * angstrom,
        },
    }

    return system, restraint_dict


@pytest.fixture(scope="session")
def mdr_restraint_fb_r0(mdr_restraint_component_fb_r0):
    """
    Generate the multiple distance restraints restraint object with
    the values for r0 and r_fb set to 0 for the permanent restraint.
    """
    system, restraint_dict = mdr_restraint_component_fb_r0

    restraint = Restraint(
        system, restraint_dict, 300 * kelvin, restraint_type="multiple_distance"
    )
    return restraint


def test_sanity_mdr(mdr_restraint):
    """Sanity check."""
    assert isinstance(mdr_restraint, Restraint)


def test_numerical_correction_mdr(mdr_restraint):
    dG = mdr_restraint.getCorrection(method="numerical") / kcal_per_mol
    assert np.isclose(-0.991, dG, atol=0.001)


def test_numerical_correction_mdr_fb0(mdr_restraint_fb_r0):
    """
    Check that the numerical multiple distance restraints correction is
    consistent with the analytical result when the flat-bottom radius is 0
    and r0 = 0.
    """
    # Get the numerical correction
    dG_numerical = mdr_restraint_fb_r0.getCorrection(method="numerical") / kcal_per_mol

    # Calculate the analytical correction
    # Constants
    v0 = (
        ((_Sire_meter3 / 1000) / _Sire_mole) / _Sire_angstrom3
    ).value()  # standard state volume in A^3
    R = (
        _k_boltz.value() * kcal_per_mol / kelvin
    ).value()  # molar gas constant in kcal mol-1 K-1

    T = mdr_restraint_fb_r0.T / kelvin  # Temperature in Kelvin

    kr = mdr_restraint_fb_r0._restraint_dict["permanent_distance_restraint"]["kr"] / (
        kcal_per_mol / angstrom**2
    )
    dG_analytical = (
        R * T * np.log((2 * np.pi * R * T / kr) ** (3 / 2) / v0)
    )  # kcal mol-1

    assert np.isclose(dG_analytical, dG_numerical, atol=0.001)


def test_analytical_correction_mdr(mdr_restraint):
    with pytest.raises(NotImplementedError):
        dG = mdr_restraint.getCorrection(method="analytical") / kcal_per_mol


class TestGromacsOutputMDR:
    @staticmethod
    @pytest.fixture(scope="class")
    def getRestraintGromacsFull(mdr_restraint):
        """The form of the restraints when the perturbation type is `full`."""
        mdr_str = mdr_restraint.toString(
            engine="GROMACS", perturbation_type="full"
        ).split("\n")
        return mdr_str

    @staticmethod
    @pytest.fixture(scope="class")
    def getRestraintGromacsRelease(mdr_restraint):
        """The form of the restraints when the perturbation type is `release_restraint`."""
        mdr_str = mdr_restraint.toString(
            engine="GROMACS", perturbation_type="release_restraint"
        ).split("\n")
        return mdr_str

    def test_n_lines(self, getRestraintGromacsFull, getRestraintGromacsRelease):
        """Test that the correct number of lines have been written."""
        assert len(getRestraintGromacsFull) == 6
        assert len(getRestraintGromacsRelease) == 8

    def test_section_names(self, getRestraintGromacsFull, getRestraintGromacsRelease):
        """Check that the section headings are correct."""
        assert (
            getRestraintGromacsFull[0].strip()
            == getRestraintGromacsRelease[0].strip()
            == "[ intermolecular_interactions ]"
        )
        assert (
            getRestraintGromacsFull[1].strip()
            == getRestraintGromacsRelease[1].strip()
            == "[ bonds ]"
        )
        assert getRestraintGromacsRelease[5].strip() == "[ distance_restraints ]"

    def test_comments(self, getRestraintGromacsFull, getRestraintGromacsRelease):
        """Check that the value label comments are as expected."""
        bond_restr_str = "; ai         aj         type       lowA       up1A       up2A       kdrA       lowB       up1B       up2B       kdrB"
        distance_restr_str = "; ai         aj         type       index      type'      low        up1        up2        fac"
        assert (
            getRestraintGromacsFull[2].strip()
            == getRestraintGromacsRelease[2].strip()
            == bond_restr_str
        )
        assert getRestraintGromacsRelease[6].strip() == distance_restr_str

    def test_at_nums(self, getRestraintGromacsFull, getRestraintGromacsRelease):
        """Check that the anchor point atom numbers are correct."""
        full_restr_lines = getRestraintGromacsFull[3:]
        release_restr_lines = getRestraintGromacsRelease[3:5] + [
            getRestraintGromacsRelease[7]
        ]
        for lines in [full_restr_lines, release_restr_lines]:
            ais = [x.split()[0] for x in lines]
            ajs = [x.split()[1] for x in lines]
            assert ais == ["1", "2", "3"]
            assert ajs == ["1496", "1497", "1498"]

    def test_restraint_vals_bond_restraint(
        self, getRestraintGromacsFull, getRestraintGromacsRelease
    ):
        """Check that the bond restraint values are correct."""
        full_restr_lines = getRestraintGromacsFull[3:]
        release_restr_lines = getRestraintGromacsRelease[3:5]

        # Strip off the atom numbers.
        full_restr_lines = [x.split()[2:] for x in full_restr_lines]
        release_restr_lines = [x.split()[2:] for x in release_restr_lines]

        # Check that all lines are the same.
        assert all([x == full_restr_lines[0] for x in full_restr_lines])
        assert all([x == release_restr_lines[0] for x in release_restr_lines])
        assert full_restr_lines[0] == release_restr_lines[0]

        # Pick one line (equivalence of all lines checked above) and check the values.
        line = [x.strip() for x in full_restr_lines[0]]
        assert line == [
            "10",
            "0.200",
            "0.400",
            "100.000",
            "0.00",
            "0.200",
            "0.400",
            "100.000",
            "4184.00",
        ]

    def test_restraint_vals_distance_restraint(self, getRestraintGromacsRelease):
        """Check that the distance restraint values are correct."""
        restr_vals = [x.strip() for x in getRestraintGromacsRelease[7].split()[2:]]
        assert restr_vals == ["2", "0", "2", "0.200", "0.400", "100.000", "1.0"]


class TestSomdOutputMDR:
    @staticmethod
    @pytest.fixture(scope="class")
    def getRestraintSomd(mdr_restraint):
        """The standard form of the restraints for any perturbation type other than `restraint`."""
        mdr_str = (
            mdr_restraint.toString(engine="SOMD", perturbation_type="release_restraint")
            .split("\n")[0]
            .split("=")[1]
            .strip()
        )
        permanent_mdr_str = (
            mdr_restraint.toString(engine="SOMD", perturbation_type="release_restraint")
            .split("\n")[1]
            .split("=")[1]
            .strip()
        )

        mdr_dict = eval(mdr_str)
        permanent_mdr_str = eval(permanent_mdr_str)
        return mdr_dict, permanent_mdr_str

    @staticmethod
    @pytest.fixture(scope="class")
    def getRestraintSomdRestraintPert(mdr_restraint):
        mdr_str = (
            mdr_restraint.toString(engine="SOMD", perturbation_type="restraint")
            .split("=")[1]
            .strip()
        )
        mdr_dict = eval(mdr_str)
        return mdr_dict

    def test_sanity(self, getRestraintSomd):
        "Sanity check"
        mdr_dict_std = getRestraintSomd[0]
        mdr_dict_permanent = getRestraintSomd[1]
        assert isinstance(mdr_dict_std, dict)
        assert isinstance(mdr_dict_permanent, dict)

    def test_dict_vals(self, getRestraintSomd):
        std_mdr_dict = getRestraintSomd[0]
        expected_anchors_std = [{0, 1495}, {1, 1496}]
        permanent_mdr_dict = getRestraintSomd[1]
        expected_anchors_permanent = [{2, 1497}]

        # Test std restraints
        for i, restraint in enumerate(std_mdr_dict):
            assert set(restraint) == expected_anchors_std[i]
            assert std_mdr_dict[restraint][0] == 3.0  # r0
            # kr is halved because of the definition of the force constant in SOMD
            assert std_mdr_dict[restraint][1] == 5.0  # kr
            assert std_mdr_dict[restraint][2] == 1.0  # r_fb

        # Test permanent restraint
        for i, restraint in enumerate(permanent_mdr_dict):
            assert set(restraint) == expected_anchors_permanent[i]
            assert permanent_mdr_dict[restraint][0] == 3.0  # r0
            # kr is halved because of the definition of the force constant in SOMD
            assert permanent_mdr_dict[restraint][1] == 5.0  # kr
            assert permanent_mdr_dict[restraint][2] == 1.0  # r_fb

    def test_dict_vals_restraint_pert(self, getRestraintSomdRestraintPert):
        expected_anchors = [{0, 1495}, {1, 1496}, {2, 1497}]

        for i, restraint in enumerate(getRestraintSomdRestraintPert):
            assert set(restraint) == expected_anchors[i]
            assert getRestraintSomdRestraintPert[restraint][0] == 3.0  # r0
            # kr is halved because of the definition of the force constant in SOMD
            assert getRestraintSomdRestraintPert[restraint][1] == 5.0  # kr
            assert getRestraintSomdRestraintPert[restraint][2] == 1.0  # r_fb
