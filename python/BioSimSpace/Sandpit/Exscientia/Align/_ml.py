__all__ = ["make_ml"]


def make_ml(molecule):
    """
    Mark the molecule as being represented as a ML ligand.

    This enables one to use

        * :meth:`~BioSimSpace.Sandpit.Exscientia._SireWrappers._system.System.getMLMolecules` to get the ML molecule.
        * :meth:`~BioSimSpace.Sandpit.Exscientia._SireWrappers._system.System.nMLMolecules` to get the number of ML molecule.
        * :meth:`~BioSimSpace.Sandpit.Exscientia._SireWrappers._molecule.Molecule.isML` to check if a molecule is a ML molecule.


    Parameters
    ----------

    molecule : BioSimSpace._SireWrappers.Molecule
        The molecule to be marked as ML ligand.

    Returns
    -------

    ML_molecule : BSS._SireWrappers.Molecule
        The molecule marked as being ML.
    """
    from .._SireWrappers import Molecule as _Molecule
    from .._Exceptions import IncompatibleError as _IncompatibleError

    # Validate input.

    if not isinstance(molecule, _Molecule):
        raise TypeError(
            "'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'"
        )

    # Cannot decouple a perturbable molecule.
    if molecule.isML():
        raise _IncompatibleError("'molecule' has already been marked as ML Molecule!")

    # Create a copy of this molecule.
    mol = _Molecule(molecule)
    mol_sire = mol._sire_object

    # Edit the molecule
    mol_edit = mol_sire.edit()

    mol_edit.set_property("ML", True)

    # Update the Sire molecule object of the new molecule.
    mol._sire_object = mol_edit.commit()

    return mol
