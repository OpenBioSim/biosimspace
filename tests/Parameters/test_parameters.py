import pytest
import tempfile

import BioSimSpace as BSS

from tests.conftest import url, has_openff, has_tleap, has_antechamber


@pytest.fixture(scope="session")
def molecule0():
    return BSS.IO.readMolecules(f"{url}/4LYT_Fixed.pdb.bz2")[0]


@pytest.fixture(scope="session")
def molecule1():
    return BSS.IO.readMolecules(f"{url}/3G8K_Fixed.pdb.bz2")[0]


@pytest.mark.skipif(has_tleap is False, reason="Requires tLEaP to be installed.")
@pytest.mark.parametrize("ff", BSS.Parameters.amberProteinForceFields())
def test_disulphide(molecule0, ff):
    """Test parameterisation in the presence of disulphide bridges."""

    # Try to parameterise with the named force field. If working, this should
    # auto-detect disulphide bonds and add the appropriate bond records to the
    # tLEaP input script.
    molecule = getattr(BSS.Parameters, ff)(molecule0).getMolecule()

    # Check that we actually generate records for four disulphide bonds.
    bonds = BSS.Parameters._Protocol.AmberProtein._get_disulphide_bonds(
        molecule0._sire_object
    )
    assert len(bonds) == 4

    # Check that the bond parameters are present in the molecule.
    bonds = molecule.search("bonds from element S to element S")
    assert len(bonds) == 4


@pytest.mark.skipif(has_tleap is False, reason="Requires tLEaP to be installed.")
def test_disulphide_renumber(molecule1, ff="ff14SB"):
    """
    Test parameterisation in the presence of disulphide bridges using a
    multi-chain PDB with duplicate residue numbering
    """

    # Try to parameterise with the named force field. If working, this should
    # auto-detect disulphide bonds and add the appropriate bond records to the
    # tLEaP input script.
    molecule = getattr(BSS.Parameters, ff)(molecule1).getMolecule()

    # Check that we actually generate records for eight disulphide bonds.
    bonds = BSS.Parameters._Protocol.AmberProtein._get_disulphide_bonds(
        molecule1._sire_object
    )
    assert len(bonds) == 8

    # Check that the bond parameters are present in the molecule.
    bonds = molecule.search("bonds from element S to element S")
    assert len(bonds) == 8


@pytest.mark.skipif(
    has_antechamber is False or has_openff is False,
    reason="Requires AmberTools/antechamber and OpenFF to be installed.",
)
def test_molecule_rename():
    """
    Test that a parameterised molecule generated from a SMILES string
    starting with the "[" character is renamed with a "smiles:" prefix
    so that it can be parsed by gmx when used in a GROMACS topology file.
    """

    # Create the parameterised molecule.
    mol = BSS.Parameters.openff_unconstrained_2_0_0(
        "[C@@H](C(F)(F)F)(OC(F)F)Cl"
    ).getMolecule()

    # Make sure the name is correct.
    assert mol._sire_object.name().value() == "smiles:[C@@H](C(F)(F)F)(OC(F)F)Cl"


@pytest.mark.skipif(
    has_antechamber is False or has_tleap is False,
    reason="Requires AmberTools/antechamber and tLEaP to be installed.",
)
def test_leap_commands(molecule0):
    """
    Test that custom commands are correctly inserted into the LEaP input
    script.
    """

    # Create lists of pre- and post-commands.
    pre_mol_commands = ["command1", "command2"]
    post_mol_commands = ["command3", "command4"]

    with tempfile.TemporaryDirectory() as tmp:
        # This will fail, but we only want to check the LEaP script.
        try:
            mol = BSS.Parameters.ff14SB(
                molecule0,
                work_dir=tmp,
                pre_mol_commands=pre_mol_commands,
                post_mol_commands=post_mol_commands,
            ).getMolecule()
        except:
            pass

        # Load the script and check that the custom parameters are present and
        # in the correct order.
        with open(f"{tmp}/leap.txt", "r") as f:
            script = f.readlines()

            # Create lists to store the indices of the custom commands.
            line_pre = [-1 for _ in range(len(pre_mol_commands))]
            line_post = [-1 for _ in range(len(post_mol_commands))]

            # Loop over the lines in the script and store the line numbers
            # where the custom commands are found.
            for x, line in enumerate(script):
                for y, command in enumerate(pre_mol_commands):
                    if command in line:
                        line_pre[y] = x
                for y, command in enumerate(post_mol_commands):
                    if command in line:
                        line_post[y] = x

            # Make sure the lines are found.
            for line in line_pre:
                assert line != -1
            for line in line_post:
                assert line != -1

            # Make sure the lines are in the correct order.
            for x in range(len(line_pre) - 1):
                assert line_pre[x] < line_pre[x + 1]
            assert line_pre[-1] < line_post[0]
            for x in range(len(line_post) - 1):
                assert line_post[x] < line_post[x + 1]


@pytest.mark.skipif(
    has_antechamber is False or has_tleap is False,
    reason="Requires AmberTools/antechamber and tLEaP to be installed.",
)
def test_smiles_stereo():
    """
    Test that SMILES string stereochemistry is correctly preserved when
    parameterising a molecule.
    """

    from rdkit import Chem

    # Define the SMILES string.
    smiles = "CC[C@@H](O)C"

    # Create the parameterised molecule.
    mol = BSS.Parameters.gaff(smiles).getMolecule()

    # Convert to RDKit format.
    rdmol0 = Chem.MolFromSmiles(smiles)
    rdmol1 = BSS.Convert.toRDKit(mol)

    # Get the SMILES string.
    rdmol0_smiles = Chem.MolToSmiles(rdmol0)
    rdmol1_smiles = Chem.MolToSmiles(Chem.RemoveHs(rdmol1))

    # Make sure the SMILES strings are the same.
    assert rdmol0_smiles == rdmol1_smiles


@pytest.mark.skipif(
    has_antechamber is False or has_tleap is False,
    reason="Requires AmberTools/antechamber and tLEaP to be installed.",
)
def test_acdoctor():
    """
    Test that parameterising negatively charged molecules works when acdoctor
    is disabled.
    """

    # Load the molecule.
    mol = BSS.IO.readMolecules(f"{url}/negative_charge.sdf")[0]

    # Make sure parameterisation fails when acdoctor is enabled.
    with pytest.raises(BSS._Exceptions.ParameterisationError):
        BSS.Parameters.gaff(mol).getMolecule()

    # Make sure parameterisation works when acdoctor is disabled.
    mol = BSS.Parameters.gaff(mol, acdoctor=False).getMolecule()
