import warnings

from .._SireWrappers import Molecule as _Molecule, System as _System


def _mark_alchemical_ion(molecule):
    """
    Mark the ion molecule as being alchemical ion.

    This enables one to use

        * :meth:`~BioSimSpace.Sandpit.Exscientia._SireWrappers._system.System.getAlchemicalIon` to get the alchemical ion.
        * :meth:`~BioSimSpace.Sandpit.Exscientia._SireWrappers._system.System.getAlchemicalIonIdx` to get the index of alchemical ion.
        * :meth:`~BioSimSpace.Sandpit.Exscientia._SireWrappers._molecule.Molecule.isAlchemicalIon` to check if a molecule is an alchemical ion.


    Parameters
    ----------

    molecule : BioSimSpace._SireWrappers.Molecule
        The molecule to be marked as alchemical ion.

    Returns
    -------

    alchemical_ion : BSS._SireWrappers.Molecule
        The molecule marked as being alchemical ion.
    """
    # Validate input.

    if not isinstance(molecule, _Molecule):
        raise TypeError(
            "'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'"
        )

    # Cannot decouple a perturbable molecule.
    if molecule.isAlchemicalIon():
        warnings.warn("'molecule' has already been marked as alchemical ion!")

    # Create a copy of this molecule.
    mol = _Molecule(molecule)
    mol_sire = mol._sire_object

    # Edit the molecule
    mol_edit = mol_sire.edit()

    mol_edit.setProperty("AlchemicalIon", True)

    # Update the Sire molecule object of the new molecule.
    mol._sire_object = mol_edit.commit()

    return mol


def _get_protein_com_idx(system: _System) -> int:
    """return the index of the atom that is closest to the center of
    mass of the biggest molecule in the system.

    Args:
        system: The input system.

    Returns:
        atom_index
    """
    biggest_mol_idx = max(range(system.nMolecules()), key=lambda x: system[x].nAtoms())

    atom_offset = 0
    for i, mol in enumerate(system):
        if biggest_mol_idx == i:
            return atom_offset + mol.getCOMIdx()
        else:
            atom_offset += mol.nAtoms()
