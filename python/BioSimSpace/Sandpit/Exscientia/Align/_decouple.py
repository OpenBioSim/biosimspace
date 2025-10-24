import warnings

from sire.legacy import Base as _SireBase
from sire.legacy import Mol as _SireMol
from sire.legacy import MM as _SireMM
from sire.legacy import Units as _SireUnits

from .._Exceptions import IncompatibleError as _IncompatibleError
from .._SireWrappers import Molecule as _Molecule

__all__ = ["decouple"]


def decouple(
    molecule, charge=(True, False), LJ=(True, False), property_map={}, intramol=True
):
    """
    Mark the molecule as being decoupled, where the interactions with the
    rest of the environment are removed, or annihilate, where the interactions
    within the molecule are removed as well (choose this mode with
    intramol=False).

    Parameters
    ----------

    molecule : BioSimSpace._SireWrappers.Molecule
        The molecule to be decoupled or annihilated.

    charge : Tuple[bool, bool]
        A Tuple of length two defines the charge of the molecule at the
        start and the end of the transformation. This allows user to
        selectively turn on or off the charge interactions
        e.g. (True, False)

    LJ : Tuple[bool, bool]
        A Tuple of length two defines the LJ of the molecule at the
        start and the end of the transformation. This allows user to
        selectively turn on or off the LJ interactions
        e.g. (True, False)

    property_map : dict
        A dictionary that maps "properties" to their user defined values.
        This allows the user to refer to properties with their own naming
        scheme, e.g. { "charge" : "my-charge" }

    intramol : bool
        Whether to couple the intra-molecule forces. Setting to ``False``
        means the intra-molecule forces will *not* be changed by the lambda
        thereby decouple the molecule instead of annihilate it.

    Returns
    -------

    decoupled : Sire.Mol.Molecule
        The molecule marked as being decoupled.
    """
    # Validate input.

    if not isinstance(molecule, _Molecule):
        raise TypeError(
            "'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'"
        )

    # Cannot decouple a perturbable molecule.
    if molecule.isDecoupled():
        raise _IncompatibleError("'molecule' has already been decoupled!")

    for field, name in zip((charge, LJ), ("charge", "LJ")):
        try:
            if len(field) != 2:
                raise ValueError(f"{name} must only have two values.")
        except TypeError:
            raise TypeError(f"{name} has to be a iterable.")

        for value in field:
            if not isinstance(value, bool):
                raise ValueError(f"{value} in {name} must be bool.")

    # Change names of charge and LJ tuples to avoid clashes with properties.
    charge_tuple = charge
    LJ_tuple = LJ

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    # Invert the user property mappings.
    inv_property_map = {v: k for k, v in property_map.items()}

    # Create a copy of this molecule and Sire object to check properties.
    mol = _Molecule(molecule)
    mol_sire = mol._sire_object

    # Get the user name for the required properties.
    ff = inv_property_map.get("forcefield", "forcefield")
    LJ = inv_property_map.get("LJ", "LJ")
    charge = inv_property_map.get("charge", "charge")
    element = inv_property_map.get("element", "element")
    ambertype = inv_property_map.get("ambertype", "ambertype")

    # Check for missing information.
    if not mol_sire.has_property(ff):
        raise _IncompatibleError("Cannot determine 'forcefield' of 'molecule'!")
    if not mol_sire.has_property(LJ):
        raise _IncompatibleError("Cannot determine LJ terms for molecule")
    if not mol_sire.has_property(charge):
        raise _IncompatibleError("Cannot determine charges for molecule")
    if not mol_sire.has_property(element):
        raise _IncompatibleError("Cannot determine elements in molecule")

    # Check for ambertype property (optional).
    has_ambertype = True
    if not mol_sire.has_property(ambertype):
        has_ambertype = False

    if not isinstance(intramol, bool):
        raise TypeError("'intramol' must be of type 'bool'")

    # Edit the molecule.
    mol_edit = mol_sire.edit()

    # Create dictionary to store charge and LJ tuples.
    mol_edit.set_property(
        "decouple", {"charge": charge_tuple, "LJ": LJ_tuple, "intramol": intramol}
    )

    # Set the "forcefield0" property.
    mol_edit.set_property("forcefield0", molecule._sire_object.property(ff))

    # Set starting properties based on fully-interacting molecule.
    mol_edit.set_property("charge0", molecule._sire_object.property(charge))
    mol_edit.set_property("LJ0", molecule._sire_object.property(LJ))
    mol_edit.set_property("element0", molecule._sire_object.property(element))
    if has_ambertype:
        mol_edit.set_property("ambertype0", molecule._sire_object.property(ambertype))

    mol_edit.set_property("annihilated", _SireBase.wrap(intramol))

    # Flag that this molecule is decoupled.
    mol_edit.set_property("is_decoupled", _SireBase.wrap(True))

    # Update the Sire molecule object of the new molecule.
    mol._sire_object = mol_edit.commit()

    return mol
