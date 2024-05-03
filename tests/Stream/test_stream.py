import os
import pytest

import BioSimSpace as BSS


@pytest.fixture(autouse=True)
def run_around_tests():
    yield
    os.remove("test.bss")


def test_system(system):
    """Test streaming of a Sire system."""

    # Stream to file.
    BSS.Stream.save(system, "test")

    # Stream from file.
    s = BSS.Stream.load("test.bss")

    # Check that both systems contain the same number of molecules.
    assert system.nMolecules() == s.nMolecules()

    # Check that the molecules in both systems contain the same number
    # of residues and atoms.
    for m0, m1 in zip(system, s):
        assert m0.nResidues() == m1.nResidues()
        assert m0.nAtoms() == m1.nAtoms()


def test_molecule(system):
    """Test streaming of a Sire molecule."""

    # Extract the first molecule.
    m0 = system[0]

    # Stream to file.
    BSS.Stream.save(m0, "test")

    # Stream from file.
    m1 = BSS.Stream.load("test.bss")

    # Check that the molecules contain the same number of residues and atoms.
    assert m0.nResidues() == m1.nResidues()
    assert m0.nAtoms() == m1.nAtoms()


def test_molecules(system):
    """Test streaming of a Sire molecule group."""

    # Stream to file.
    BSS.Stream.save(system.getMolecules(), "test")

    # Stream from file.
    m = BSS.Stream.load("test.bss")

    # Check that the system contain the same number of molecules as the
    # molecule group.
    assert system.nMolecules() == len(m)

    # Check that the molecules in the system and molecule group contain the
    # same number of residues and atoms.
    for m0, m1 in zip(system, m):
        assert m0.nResidues() == m1.nResidues()
        assert m0.nAtoms() == m1.nAtoms()


def test_residue(system):
    """Test streaming of a Sire residue."""

    # Extract the first residue of the first molecule.
    r0 = system[0].getResidues()[0]

    # Stream to file.
    BSS.Stream.save(r0, "test")

    # Stream from file.
    r1 = BSS.Stream.load("test.bss")

    # Check that the residues contain the same number of atoms.
    assert r0.nAtoms() == r1.nAtoms()


def test_atom(system):
    """Test streaming of a Sire atom."""

    # Extract the first atom of the first molecule.
    a0 = system[0].getAtoms()[0]

    # Stream to file.
    BSS.Stream.save(a0, "test")

    # Stream from file.
    a1 = BSS.Stream.load("test.bss")

    # Check that the atom elements are the same.
    assert a0.element() == a1.element()


@pytest.mark.xfail(reason="Sire streaming for Selector<T> types is not yet exposed.")
def test_select_result(system):
    """Test streaming of a Sire select result."""

    # Search for all carbon atoms.
    s0 = system.search("element C")

    # Stream to file.
    BSS.Stream.save(s0, "test")

    # Stream from file.
    s1 = BSS.Stream.load("test.bss")

    # Check that the number of search results is the same.
    assert len(s0.atoms()) == len(s1.atoms())

    # Check that the results have the same element and index.
    for a0, a1 in zip(s0.atoms(), s1.atoms()):
        assert a0.index() == a1.index()
        assert a0.element() == a1.element()


def test_metadata(system):
    """Test that streamed metadata is saved and recovered correctly."""

    # Add the metadata for the system.
    original_metadata = BSS.Stream._stream._add_metadata(system)

    # Stream to file.
    BSS.Stream.save(system, "test")

    # Query the metadata.
    metadata = BSS.Stream.getMetadata("test.bss")

    # Make sure that the metadata is the same.
    for k, v in metadata.items():
        assert original_metadata[k] == v
