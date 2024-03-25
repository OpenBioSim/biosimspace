import glob
import os
import pytest
import tempfile

import BioSimSpace.Sandpit.Exscientia as BSS

from tests.Sandpit.Exscientia.conftest import has_amber, has_openff
from tests.conftest import root_fp


def test_file_cache():
    """
    Simple test to see if cached files are used when repeatedly writing
    the same system to the same file format.
    """

    # Clear the file cache.
    BSS.IO.clearCache()

    # Load the molecular system.
    s = BSS.IO.readMolecules([f"{root_fp}/input/ala.crd", f"{root_fp}/input/ala.top"])

    # Create a temporary working directory.
    tmp_dir = tempfile.TemporaryDirectory()
    tmp_path = tmp_dir.name

    # Write the system to PDB and PRM7 format.
    BSS.IO.saveMolecules(f"{tmp_path}/tmp", s, ["pdb", "prm7"])

    # Check that the file cache has two entries.
    assert len(BSS.IO._file_cache._cache) == 2

    # Write to PDB and GroTop format. The PDB from the cache should be reused.
    BSS.IO.saveMolecules(f"{tmp_path}/tmp2", s, ["pdb", "grotop"])

    # Check that the file cache has three entries, i.e. only the GroTop was added.
    assert len(BSS.IO._file_cache._cache) == 3

    # The directory should now contain 4 files.
    assert len(glob.glob(f"{tmp_path}/*")) == 4

    # Now delete one of the files on disk and re-write. This should create a new
    # entry in the cache.
    os.remove(f"{tmp_path}/tmp.pdb")
    BSS.IO.saveMolecules(f"{tmp_path}/tmp3", s, ["pdb", "grotop"])

    # Check that the file cache still has three entries.
    assert len(BSS.IO._file_cache._cache) == 3

    # The directory should now contain 5 files.
    assert len(glob.glob(f"{tmp_path}/*")) == 5

    # Now "corrupt" a file on disk so that its MD5 checksum is no longer
    # valid.
    with open(f"{tmp_path}/tmp2.pdb", "w") as f:
        pass

    # Write back to PDB and Gro87 format. The PDB file is now invalid, so
    # a new one will be written and added to the cache.
    BSS.IO.saveMolecules(f"{tmp_path}/tmp4", s, ["pdb", "gro87"])

    # Check that the file cache still has three entries.
    assert len(BSS.IO._file_cache._cache) == 4

    # The directory should now contain 7 files.
    assert len(glob.glob(f"{tmp_path}/*")) == 7

    # Now delete a key and check that the number of atoms is decremented.

    # Get the first key.
    key = list(BSS.IO._file_cache._cache.keys())[0]

    # Store the number of atoms for this key.
    num_atoms = BSS.IO._file_cache._cache[key][0].nAtoms()

    # Store the total number of atoms in the cache.
    total_atoms = BSS.IO._file_cache._cache._num_atoms

    # Delete the key.
    del BSS.IO._file_cache._cache[key]

    # Make sure the number of atoms in the cache was decremented.
    assert BSS.IO._file_cache._cache._num_atoms == total_atoms - num_atoms

    # Clear the file cache.
    BSS.IO.clearCache()

    # The cache should now be empty.
    assert len(BSS.IO._file_cache._cache) == 0

    # Disable the cache.
    BSS.IO.disableCache()

    # Write to PDB and GroTop format. The PDB from the cache should not be reused.
    BSS.IO.saveMolecules(f"{tmp_path}/tmp5", s, ["pdb", "grotop"])

    # The cache should still be empty.
    assert len(BSS.IO._file_cache._cache) == 0


@pytest.mark.skipif(
    has_amber is False or has_openff is False,
    reason="Requires AMBER and OpenFF to be installed",
)
def test_file_cache_mol_nums():
    """
    Make sure that systems can be cached if they have the same UID, but
    contain different MolNUms.
    """

    # Enable the cache.
    BSS.IO.enableCache()

    # Clear the file cache.
    BSS.IO.clearCache()

    # Create an initial system.
    system = BSS.Parameters.openff_unconstrained_2_0_0("CO").getMolecule().toSystem()

    # Create two different 5 atom molecules.
    mol0 = BSS.Parameters.openff_unconstrained_2_0_0("C").getMolecule()
    mol1 = BSS.Parameters.openff_unconstrained_2_0_0("CF").getMolecule()

    # Create two new systems by adding the different molecules to the original
    # system. These will have the same UID, but different molecule numbers.
    system0 = system + mol0
    system1 = system + mol1

    # Create a temporary working directory.
    tmp_dir = tempfile.TemporaryDirectory()
    tmp_path = tmp_dir.name

    # Write the two systems to PDB format.
    BSS.IO.saveMolecules(f"{tmp_path}/tmp0", system0, "pdb")
    BSS.IO.saveMolecules(f"{tmp_path}/tmp1", system1, "pdb")

    # The cache shold have two entries.
    assert len(BSS.IO._file_cache._cache) == 2
