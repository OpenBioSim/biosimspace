import BioSimSpace as BSS

from sire.legacy import Mol as SireMol

import pytest


def test_system(system):
    """
    Check that system conversions work as expected.
    """

    # Store the number of molecules, residues, and atoms.
    num_mol = system.nMolecules()
    num_res = system.nResidues()
    num_atm = system.nAtoms()

    # Convert to Sire format.
    sire_system = BSS.Convert.toSire(system)

    assert sire_system.num_molecules() == num_mol
    assert sire_system.num_residues() == num_res
    assert sire_system.num_atoms() == num_atm

    # Convert to RDKit format. The result will be a container of
    # RDKit molecules.
    rdmols = BSS.Convert.toRDKit(sire_system)

    assert len(rdmols) == num_mol

    # Work out atom total. (No concept of a residue.)
    nats = 0
    for rdmol in rdmols:
        nats += rdmol.GetNumAtoms()

    assert nats == num_atm

    # Convert back to BioSimSpace.
    bss_system = BSS.Convert.toBioSimSpace(rdmols).toSystem()

    assert bss_system.nMolecules() == num_mol
    assert bss_system.nAtoms() == num_atm


def test_residue(system):
    """
    Check that we can convert single residues.
    """

    # Convert to Sire format.
    sire_residue = BSS.Convert.toSire(system[0].getResidues()[0])

    assert isinstance(sire_residue, SireMol.Residue)

    # Now convert back.
    bss_residue = BSS.Convert.toBioSimSpace(sire_residue)

    assert isinstance(bss_residue, BSS._SireWrappers.Residue)


def test_atom(system):
    """
    Check that we can convert single atoms.
    """

    # Convert to Sire format.
    sire_atom = BSS.Convert.toSire(system[0].getAtoms()[0])

    assert isinstance(sire_atom, SireMol.Atom)

    # Now convert back.
    bss_atom = BSS.Convert.toBioSimSpace(sire_atom)

    assert isinstance(bss_atom, BSS._SireWrappers.Atom)


def test_molecule_rename():
    """
    Test that a parameterised molecule generated from a SMILES string
    starting with the "[" character is renamed with a "smiles:" prefix
    so that it can be parsed by gmx when used in a GROMACS topology file.
    """

    # Create the parameterised molecule.
    mol = BSS.Convert.smiles("[C@@H](C(F)(F)F)(OC(F)F)Cl")

    # Make sure the name is correct.
    assert mol._sire_object.name().value() == "smiles:[C@@H](C(F)(F)F)(OC(F)F)Cl"


@pytest.mark.parametrize("lig", ["2dd", "2hh", "2ii", "2s", "2x"])
def test_sdf_stereo(lig):
    """
    Test that SDF stereochemistry is correctly preserved when converting
    to RDKit format via Sire.
    """

    import tempfile

    from rdkit import Chem
    from tests.conftest import url

    # Load the molecule.
    mol = BSS.IO.readMolecules(f"{url}/lig_{lig}.sdf")[0]

    with tempfile.NamedTemporaryFile() as f:
        # Write the molecule to a temporary file.
        BSS.IO.saveMolecules(f.name, mol, "sdf")

        # Load the molecule using RDKit.
        rdmol = Chem.SDMolSupplier(f"{f.name}.sdf")[0]

        # Convert the molecule using Sire.
        rdmol_sire = BSS.Convert.toRDKit(mol)

        # Get the SMILES strings.
        smiles = Chem.CanonSmiles(Chem.MolToSmiles(rdmol))
        smiles_sire = Chem.CanonSmiles(Chem.MolToSmiles(rdmol_sire))

        # Check that the SMILES strings are the same.
        assert smiles == smiles_sire
