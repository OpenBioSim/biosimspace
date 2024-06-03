######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2024
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BioSimSpace is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""Functionality for merging molecules."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["merge"]

from sire.legacy import Base as _SireBase
from sire.legacy import IO as _SireIO
from sire.legacy import MM as _SireMM
from sire.legacy import Mol as _SireMol
from sire.legacy import Units as _SireUnits

from .._Exceptions import IncompatibleError as _IncompatibleError
from .._SireWrappers import Molecule as _Molecule


def merge(
    molecule0,
    molecule1,
    mapping,
    allow_ring_breaking=False,
    allow_ring_size_change=False,
    force=False,
    property_map0={},
    property_map1={},
):
    """
    Merge this molecule with 'other'.

    Parameters
    ----------

    molecule0 : BioSimSpace._SireWrappers.Molecule
        The reference molecule.

    molecule1 : BioSimSpace._SireWrappers.Molecule
        The molecule to merge with.

    mapping : dict
        The mapping between matching atom indices in the two molecules.

    allow_ring_breaking : bool
        Whether to allow the opening/closing of rings during a merge.

    allow_ring_size_change : bool
        Whether to allow changes in ring size.

    force : bool
        Whether to try to force the merge, even when the molecular
        connectivity changes not as the result of a ring transformation.
        This will likely lead to an unstable perturbation. This option
        takes precedence over 'allow_ring_breaking' and
        'allow_ring_size_change'.

    property_map0 : dict
        A dictionary that maps "properties" in this molecule to their
        user defined values. This allows the user to refer to properties
        with their own naming scheme, e.g. { "charge" : "my-charge" }

    property_map1 : dict
        A dictionary that maps "properties" in there other molecule to
        their user defined values.

    Returns
    -------

    merged : Sire.Mol.Molecule
        The merged molecule.
    """

    # Validate input.

    if not isinstance(molecule0, _Molecule):
        raise TypeError(
            "'molecule0' must be of type 'BioSimSpace._SireWrappers.Molecule'"
        )

    if not isinstance(molecule1, _Molecule):
        raise TypeError(
            "'molecule0' must be of type 'BioSimSpace._SireWrappers.Molecule'"
        )

    # Cannot merge a perturbable molecule.
    if molecule0._is_perturbable:
        raise _IncompatibleError("'molecule0' has already been merged!")
    if molecule1._is_perturbable:
        raise _IncompatibleError("'molecule1' has already been merged!")

    if not isinstance(property_map0, dict):
        raise TypeError("'property_map0' must be of type 'dict'")

    if not isinstance(property_map1, dict):
        raise TypeError("'property_map1' must be of type 'dict'")

    if not isinstance(allow_ring_breaking, bool):
        raise TypeError("'allow_ring_breaking' must be of type 'bool'")

    if not isinstance(allow_ring_size_change, bool):
        raise TypeError("'allow_ring_size_change' must be of type 'bool'")

    if not isinstance(force, bool):
        raise TypeError("'force' must be of type 'bool'")

    if not isinstance(mapping, dict):
        raise TypeError("'mapping' must be of type 'dict'.")
    else:
        # Make sure all key/value pairs are of type AtomIdx.
        for idx0, idx1 in mapping.items():
            if not isinstance(idx0, _SireMol.AtomIdx) or not isinstance(
                idx1, _SireMol.AtomIdx
            ):
                raise TypeError(
                    "key:value pairs in 'mapping' must be of type 'Sire.Mol.AtomIdx'"
                )

    # Set 'allow_ring_breaking' and 'allow_ring_size_change' to true if the
    # user has requested to 'force' the merge, i.e. the 'force' argument
    # takes precedence.
    if force:
        allow_ring_breaking = True
        allow_ring_size_change = True

    # Create a copy of this molecule.
    mol = _Molecule(molecule0)

    # Extract the two Sire molecule objects.
    molecule0 = molecule0._sire_object
    molecule1 = molecule1._sire_object

    # Get the atom indices from the mapping.
    idx0 = mapping.keys()
    idx1 = mapping.values()

    # Create the reverse mapping: molecule1 --> molecule0
    inv_mapping = {v: k for k, v in mapping.items()}

    # Invert the user property mappings.
    inv_property_map0 = {v: k for k, v in property_map0.items()}
    inv_property_map1 = {v: k for k, v in property_map1.items()}

    # Make sure that the molecules have a "forcefield" property and that
    # the two force fields are compatible.

    # Get the user name for the "forcefield" property.
    ff0 = inv_property_map0.get("forcefield", "forcefield")
    ff1 = inv_property_map1.get("forcefield", "forcefield")

    # Force field information is missing.
    if not molecule0.hasProperty(ff0):
        raise _IncompatibleError("Cannot determine 'forcefield' of 'molecule0'!")
    if not molecule1.hasProperty(ff1):
        raise _IncompatibleError("Cannot determine 'forcefield' of 'molecule1'!")

    # The force fields are incompatible.
    if not molecule0.property(ff0).isCompatibleWith(molecule1.property(ff1)):
        raise _IncompatibleError(
            "Cannot merge molecules with incompatible force fields!"
        )

    # Create lists to store the atoms that are unique to each molecule,
    # along with their indices.
    atoms0 = []
    atoms1 = []
    atoms0_idx = []
    atoms1_idx = []

    # Loop over each molecule to find the unique atom indices.

    # molecule0
    for atom in molecule0.atoms():
        if atom.index() not in idx0:
            atoms0.append(atom)
            atoms0_idx.append(atom.index())

    # molecule1
    for atom in molecule1.atoms():
        if atom.index() not in idx1:
            atoms1.append(atom)
            atoms1_idx.append(atom.index())

    # Create a new molecule to hold the merged molecule.
    molecule = _SireMol.Molecule("Merged_Molecule")

    # Add a single residue called LIG.
    res = molecule.edit().add(_SireMol.ResNum(1))
    res.rename(_SireMol.ResName("LIG"))

    # Create a single cut-group.
    cg = res.molecule().add(_SireMol.CGName("1"))

    # Counter for the number of atoms.
    num = 1

    # First add all of the atoms from molecule0.
    for atom in molecule0.atoms():
        # Add the atom.
        added = cg.add(atom.name())
        added.renumber(_SireMol.AtomNum(num))
        added.reparent(_SireMol.ResIdx(0))
        num += 1

    # Now add all of the atoms from molecule1 that aren't mapped from molecule0.
    for atom in atoms1:
        added = cg.add(atom.name())
        added.renumber(_SireMol.AtomNum(num))
        added.reparent(_SireMol.ResIdx(0))
        inv_mapping[atom.index()] = _SireMol.AtomIdx(num - 1)
        num += 1

    # Commit the changes to the molecule.
    molecule = cg.molecule().commit()

    # Make the molecule editable.
    edit_mol = molecule.edit()

    # Create lists of the actual property names in the molecules.
    props0 = []
    props1 = []

    # molecule0
    for prop in molecule0.propertyKeys():
        if prop in inv_property_map0:
            prop = inv_property_map0[prop]
        props0.append(prop)

    # molecule1
    for prop in molecule1.propertyKeys():
        if prop in inv_property_map1:
            prop = inv_property_map1[prop]
        props1.append(prop)

    # Determine the common properties between the two molecules.
    # These are the properties that can be perturbed.
    shared_props = list(set(props0).intersection(props1))

    # Now add default properties in each end state. This is required since Sire
    # removes units from atom properties that are zero valued. By pre-setting
    # a default property for the entire molecule, we ensure that atom properties
    # are added with the correct type.

    # Properties to ignore. (These will be handled later.)
    ignored_props = [
        "bond",
        "angle",
        "dihedral",
        "improper",
        "connectivity",
        "intrascale",
    ]

    # lambda = 0
    for prop in props0:
        if not prop in ignored_props:
            # This is a perturbable property.
            if prop in shared_props:
                name = f"{prop}0"
            # This property is unique to the lambda = 0 state.
            else:
                name = prop

            # Get the property from molecule0.
            property = molecule0.property(prop)

            # Try to set a default property at the lambda = 0 end state.
            try:
                default_prop = type(property)(molecule.info())
                edit_mol = edit_mol.setProperty(name, default_prop).molecule()
            except:
                pass

    # lambda = 1
    for prop in props1:
        if not prop in ignored_props:
            # This is a perturbable property.
            if prop in shared_props:
                name = f"{prop}1"
            # This property is unique to the lambda = 1 state.
            else:
                name = prop

            # Get the property from molecule0.
            property = molecule1.property(prop)

            # Try to set a default property at the lambda = 1 end state.
            try:
                default_prop = type(property)(molecule.info())
                edit_mol = edit_mol.setProperty(name, default_prop).molecule()
            except:
                pass

    del props0
    del props1

    # We now add properties to the merged molecule. The properties are used
    # to represent the molecule at two states along the alchemical pathway:
    #
    # lambda = 0:
    #   Here only molecule0 is active. The "charge" and "LJ" properties
    #   for atoms that are part of molecule1 are set to zero. Other force
    #   field properties, e.g. "bond", "angle", "dihedral", and "improper",
    #   are retained for the atoms in molecule1, although the indices of
    #   the atoms involved in the interactions must be re-mapped to their
    #   positions in the merged molecule. (Also, the interactions may now
    #   be between atoms of different type.) The properties are given the
    #   suffix "0", e.g. "charge0".
    #
    # lambda = 1:
    #   Here only molecule1 is active. We perform the same process as above,
    #   only we modify the properties of the atoms that are unique to
    #   molecule0. The properties are given the suffix "1", e.g. "charge1".
    #
    # Properties that aren't shared between the molecules (and thus can't
    # be merged) are set using their original names.

    ##############################
    # SET PROPERTIES AT LAMBDA = 0
    ##############################

    # Add the atom properties from molecule0.
    for atom in molecule0.atoms():
        # Add an "name0" property.
        edit_mol = (
            edit_mol.atom(atom.index())
            .setProperty("name0", atom.name().value())
            .molecule()
        )

        # Loop over all atom properties.
        for prop in atom.propertyKeys():
            # Get the actual property name.
            name = inv_property_map0.get(prop, prop)

            # This is a perturbable property. Rename to "property0", e.g. "charge0".
            if name in shared_props:
                name = name + "0"

            # Add the property to the atom in the merged molecule.
            edit_mol = (
                edit_mol.atom(atom.index())
                .setProperty(name, atom.property(prop))
                .molecule()
            )

    # Add the atom properties from molecule1.
    for atom in atoms1:
        # Get the atom index in the merged molecule.
        idx = inv_mapping[atom.index()]

        # Add an "name0" property.
        edit_mol = (
            edit_mol.atom(idx).setProperty("name0", atom.name().value()).molecule()
        )

        # Loop over all atom properties.
        for prop in atom.propertyKeys():
            # Get the actual property name.
            name = inv_property_map1.get(prop, prop)

            # Zero the "charge" and "LJ" property for atoms that are unique to molecule1.
            if name == "charge":
                edit_mol = (
                    edit_mol.atom(idx)
                    .setProperty("charge0", 0 * _SireUnits.e_charge)
                    .molecule()
                )
            elif name == "LJ":
                edit_mol = (
                    edit_mol.atom(idx)
                    .setProperty("LJ0", _SireMM.LJParameter())
                    .molecule()
                )
            elif name == "ambertype":
                edit_mol = edit_mol.atom(idx).setProperty("ambertype0", "du").molecule()
            elif name == "element":
                edit_mol = (
                    edit_mol.atom(idx)
                    .setProperty("element0", _SireMol.Element(0))
                    .molecule()
                )
            else:
                # This is a perturbable property. Rename to "property0", e.g. "charge0".
                if name in shared_props:
                    name = name + "0"

                # Add the property to the atom in the merged molecule.
                edit_mol = (
                    edit_mol.atom(idx).setProperty(name, atom.property(prop)).molecule()
                )

    # We now need to merge "bond", "angle", "dihedral", and "improper" parameters.
    # To do so, we extract the properties from molecule0, then add the additional
    # properties from molecule1, making sure to update the atom indices, and bond
    # atoms from molecule1 to the atoms to which they map in molecule0.

    # 1) bonds
    if "bond" in shared_props:
        # Get the user defined property names.
        prop0 = inv_property_map0.get("bond", "bond")
        prop1 = inv_property_map1.get("bond", "bond")

        # Get the "bond" property from the two molecules.
        bonds0 = molecule0.property(prop0)
        bonds1 = molecule1.property(prop1)

        # Get the molInfo object for each molecule.
        info0 = molecule0.info()
        info1 = molecule1.info()

        # Create the new set of bonds.
        bonds = _SireMM.TwoAtomFunctions(edit_mol.info())

        # Add all of the bonds from molecule0.
        for bond in bonds0.potentials():
            atom0 = info0.atomIdx(bond.atom0())
            atom1 = info0.atomIdx(bond.atom1())
            bonds.set(atom0, atom1, bond.function())

        # Loop over all bonds in molecule1.
        for bond in bonds1.potentials():
            # This bond contains an atom that is unique to molecule1.
            if (
                info1.atomIdx(bond.atom0()) in atoms1_idx
                or info1.atomIdx(bond.atom1()) in atoms1_idx
            ):
                # Extract the bond information.
                atom0 = info1.atomIdx(bond.atom0())
                atom1 = info1.atomIdx(bond.atom1())
                exprn = bond.function()

                # Map the atom indices to their position in the merged molecule.
                atom0 = inv_mapping[atom0]
                atom1 = inv_mapping[atom1]

                # Set the new bond.
                bonds.set(atom0, atom1, exprn)

        # Add the bonds to the merged molecule.
        edit_mol.setProperty("bond0", bonds)

    # 2) angles
    if "angle" in shared_props:
        # Get the user defined property names.
        prop0 = inv_property_map0.get("angle", "angle")
        prop1 = inv_property_map1.get("angle", "angle")

        # Get the "angle" property from the two molecules.
        angles0 = molecule0.property(prop0)
        angles1 = molecule1.property(prop1)

        # Get the molInfo object for each molecule.
        info0 = molecule0.info()
        info1 = molecule1.info()

        # Create the new set of angles.
        angles = _SireMM.ThreeAtomFunctions(edit_mol.info())

        # Add all of the angles from molecule0.
        for angle in angles0.potentials():
            atom0 = info0.atomIdx(angle.atom0())
            atom1 = info0.atomIdx(angle.atom1())
            atom2 = info0.atomIdx(angle.atom2())
            angles.set(atom0, atom1, atom2, angle.function())

        # Loop over all angles in molecule1.
        for angle in angles1.potentials():
            # This angle contains an atom that is unique to molecule1.
            if (
                info1.atomIdx(angle.atom0()) in atoms1_idx
                or info1.atomIdx(angle.atom1()) in atoms1_idx
                or info1.atomIdx(angle.atom2()) in atoms1_idx
            ):
                # Extract the angle information.
                atom0 = info1.atomIdx(angle.atom0())
                atom1 = info1.atomIdx(angle.atom1())
                atom2 = info1.atomIdx(angle.atom2())
                exprn = angle.function()

                # Map the atom indices to their position in the merged molecule.
                atom0 = inv_mapping[atom0]
                atom1 = inv_mapping[atom1]
                atom2 = inv_mapping[atom2]

                # Set the new angle.
                angles.set(atom0, atom1, atom2, exprn)

        # Add the angles to the merged molecule.
        edit_mol.setProperty("angle0", angles)

    # 3) dihedrals
    if "dihedral" in shared_props:
        # Get the user defined property names.
        prop0 = inv_property_map0.get("dihedral", "dihedral")
        prop1 = inv_property_map1.get("dihedral", "dihedral")

        # Get the "dihedral" property from the two molecules.
        dihedrals0 = molecule0.property(prop0)
        dihedrals1 = molecule1.property(prop1)

        # Get the molInfo object for each molecule.
        info0 = molecule0.info()
        info1 = molecule1.info()

        # Create the new set of dihedrals.
        dihedrals = _SireMM.FourAtomFunctions(edit_mol.info())

        # Add all of the dihedrals from molecule0.
        for dihedral in dihedrals0.potentials():
            atom0 = info0.atomIdx(dihedral.atom0())
            atom1 = info0.atomIdx(dihedral.atom1())
            atom2 = info0.atomIdx(dihedral.atom2())
            atom3 = info0.atomIdx(dihedral.atom3())
            dihedrals.set(atom0, atom1, atom2, atom3, dihedral.function())

        # Loop over all dihedrals in molecule1.
        for dihedral in dihedrals1.potentials():
            # This dihedral contains an atom that is unique to molecule1.
            if (
                info1.atomIdx(dihedral.atom0()) in atoms1_idx
                or info1.atomIdx(dihedral.atom1()) in atoms1_idx
                or info1.atomIdx(dihedral.atom2()) in atoms1_idx
                or info1.atomIdx(dihedral.atom3()) in atoms1_idx
            ):
                # Extract the dihedral information.
                atom0 = info1.atomIdx(dihedral.atom0())
                atom1 = info1.atomIdx(dihedral.atom1())
                atom2 = info1.atomIdx(dihedral.atom2())
                atom3 = info1.atomIdx(dihedral.atom3())
                exprn = dihedral.function()

                # Map the atom indices to their position in the merged molecule.
                atom0 = inv_mapping[atom0]
                atom1 = inv_mapping[atom1]
                atom2 = inv_mapping[atom2]
                atom3 = inv_mapping[atom3]

                # Set the new dihedral.
                dihedrals.set(atom0, atom1, atom2, atom3, exprn)

        # Add the dihedrals to the merged molecule.
        edit_mol.setProperty("dihedral0", dihedrals)

    # 4) impropers
    if "improper" in shared_props:
        # Get the user defined property names.
        prop0 = inv_property_map0.get("improper", "improper")
        prop1 = inv_property_map1.get("improper", "improper")

        # Get the "improper" property from the two molecules.
        impropers0 = molecule0.property(prop0)
        impropers1 = molecule1.property(prop1)

        # Get the molInfo object for each molecule.
        info0 = molecule0.info()
        info1 = molecule1.info()

        # Create the new set of impropers.
        impropers = _SireMM.FourAtomFunctions(edit_mol.info())

        # Add all of the impropers from molecule0.
        for improper in impropers0.potentials():
            atom0 = info0.atomIdx(improper.atom0())
            atom1 = info0.atomIdx(improper.atom1())
            atom2 = info0.atomIdx(improper.atom2())
            atom3 = info0.atomIdx(improper.atom3())
            impropers.set(atom0, atom1, atom2, atom3, improper.function())

        # Loop over all impropers in molecule1.
        for improper in impropers1.potentials():
            # This improper contains an atom that is unique to molecule1.
            if (
                info1.atomIdx(improper.atom0()) in atoms1_idx
                or info1.atomIdx(improper.atom1()) in atoms1_idx
                or info1.atomIdx(improper.atom2()) in atoms1_idx
                or info1.atomIdx(improper.atom3()) in atoms1_idx
            ):
                # Extract the improper information.
                atom0 = info1.atomIdx(improper.atom0())
                atom1 = info1.atomIdx(improper.atom1())
                atom2 = info1.atomIdx(improper.atom2())
                atom3 = info1.atomIdx(improper.atom3())
                exprn = improper.function()

                # Map the atom indices to their position in the merged molecule.
                atom0 = inv_mapping[atom0]
                atom1 = inv_mapping[atom1]
                atom2 = inv_mapping[atom2]
                atom3 = inv_mapping[atom3]

                # Set the new improper.
                impropers.set(atom0, atom1, atom2, atom3, exprn)

        # Add the impropers to the merged molecule.
        edit_mol.setProperty("improper0", impropers)

    ##############################
    # SET PROPERTIES AT LAMBDA = 1
    ##############################

    # Add the atom properties from molecule1.
    for atom in molecule1.atoms():
        # Get the atom index in the merged molecule.
        idx = inv_mapping[atom.index()]

        # Add an "name1" property.
        edit_mol = (
            edit_mol.atom(idx).setProperty("name1", atom.name().value()).molecule()
        )

        # Loop over all atom properties.
        for prop in atom.propertyKeys():
            # Get the actual property name.
            name = inv_property_map1.get(prop, prop)

            # This is a perturbable property. Rename to "property1", e.g. "charge1".
            if name in shared_props:
                name = name + "1"

            # Add the property to the atom in the merged molecule.
            edit_mol = (
                edit_mol.atom(idx).setProperty(name, atom.property(prop)).molecule()
            )

    # Add the properties from atoms unique to molecule0.
    for atom in atoms0:
        # Add an "name1" property.
        edit_mol = (
            edit_mol.atom(atom.index())
            .setProperty("name1", atom.name().value())
            .molecule()
        )

        # Loop over all atom properties.
        for prop in atom.propertyKeys():
            # Get the actual property name.
            name = inv_property_map0.get(prop, prop)

            # Zero the "charge" and "LJ" property for atoms that are unique to molecule0.
            if name == "charge":
                edit_mol = (
                    edit_mol.atom(atom.index())
                    .setProperty("charge1", 0 * _SireUnits.e_charge)
                    .molecule()
                )
            elif name == "LJ":
                edit_mol = (
                    edit_mol.atom(atom.index())
                    .setProperty("LJ1", _SireMM.LJParameter())
                    .molecule()
                )
            elif name == "ambertype":
                edit_mol = (
                    edit_mol.atom(atom.index())
                    .setProperty("ambertype1", "du")
                    .molecule()
                )
            elif name == "element":
                edit_mol = (
                    edit_mol.atom(atom.index())
                    .setProperty("element1", _SireMol.Element(0))
                    .molecule()
                )
            else:
                # This is a perturbable property. Rename to "property1", e.g. "charge1".
                if name in shared_props:
                    name = name + "1"

                # Add the property to the atom in the merged molecule.
                edit_mol = (
                    edit_mol.atom(atom.index())
                    .setProperty(name, atom.property(prop))
                    .molecule()
                )

    # Create a null LJParameter.
    null_lj = _SireMM.LJParameter()

    # Atoms with zero LJ sigma values need to have their sigma values set to the
    # value from the other end state.
    for atom in edit_mol.atoms():
        # Get the end state LJ sigma values.
        lj0 = atom.property("LJ0")
        lj1 = atom.property("LJ1")

        # Lambda = 0 state has a zero sigma value.
        if lj0.sigma() == null_lj.sigma():
            # Use the sigma value from the lambda = 1 state.
            edit_mol = (
                edit_mol.atom(atom.index())
                .set_property("LJ0", _SireMM.LJParameter(lj1.sigma(), lj0.epsilon()))
                .molecule()
            )

        # Lambda = 1 state has a zero sigma value.
        if lj1.sigma() == null_lj.sigma():
            # Use the sigma value from the lambda = 0 state.
            edit_mol = (
                edit_mol.atom(atom.index())
                .set_property("LJ1", _SireMM.LJParameter(lj0.sigma(), lj1.epsilon()))
                .molecule()
            )

    # We now need to merge "bond", "angle", "dihedral", and "improper" parameters.
    # To do so, we extract the properties from molecule1, then add the additional
    # properties from molecule0, making sure to update the atom indices, and bond
    # atoms from molecule0 to the atoms to which they map in molecule1.

    # 1) bonds
    if "bond" in shared_props:
        # Get the info objects for the two molecules.
        info0 = molecule0.info()
        info1 = molecule1.info()

        # Get the user defined property names.
        prop0 = inv_property_map0.get("bond", "bond")
        prop1 = inv_property_map1.get("bond", "bond")

        # Get the "bond" property from the two molecules.
        bonds0 = molecule0.property(prop0)
        bonds1 = molecule1.property(prop1)

        # Create the new set of bonds.
        bonds = _SireMM.TwoAtomFunctions(edit_mol.info())

        # Add all of the bonds from molecule1.
        for bond in bonds1.potentials():
            # Extract the bond information.
            atom0 = info1.atomIdx(bond.atom0())
            atom1 = info1.atomIdx(bond.atom1())
            exprn = bond.function()

            # Map the atom indices to their position in the merged molecule.
            atom0 = inv_mapping[atom0]
            atom1 = inv_mapping[atom1]

            # Set the new bond.
            bonds.set(atom0, atom1, exprn)

        # Loop over all bonds in molecule0
        for bond in bonds0.potentials():
            # This bond contains an atom that is unique to molecule0.
            if (
                info0.atomIdx(bond.atom0()) in atoms0_idx
                or info0.atomIdx(bond.atom1()) in atoms0_idx
            ):
                # Extract the bond information.
                atom0 = info0.atomIdx(bond.atom0())
                atom1 = info0.atomIdx(bond.atom1())
                exprn = bond.function()

                # Set the new bond.
                bonds.set(atom0, atom1, exprn)

        # Add the bonds to the merged molecule.
        edit_mol.setProperty("bond1", bonds)

    # 2) angles
    if "angle" in shared_props:
        # Get the info objects for the two molecules.
        info0 = molecule0.info()
        info1 = molecule1.info()

        # Get the user defined property names.
        prop0 = inv_property_map0.get("angle", "angle")
        prop1 = inv_property_map1.get("angle", "angle")

        # Get the "angle" property from the two molecules.
        angles0 = molecule0.property(prop0)
        angles1 = molecule1.property(prop1)

        # Create the new set of angles.
        angles = _SireMM.ThreeAtomFunctions(edit_mol.info())

        # Add all of the angles from molecule1.
        for angle in angles1.potentials():
            # Extract the angle information.
            atom0 = info1.atomIdx(angle.atom0())
            atom1 = info1.atomIdx(angle.atom1())
            atom2 = info1.atomIdx(angle.atom2())
            exprn = angle.function()

            # Map the atom indices to their position in the merged molecule.
            atom0 = inv_mapping[atom0]
            atom1 = inv_mapping[atom1]
            atom2 = inv_mapping[atom2]

            # Set the new angle.
            angles.set(atom0, atom1, atom2, exprn)

        # Loop over all angles in molecule0.
        for angle in angles0.potentials():
            # This angle contains an atom that is unique to molecule0.
            if (
                info0.atomIdx(angle.atom0()) in atoms0_idx
                or info0.atomIdx(angle.atom1()) in atoms0_idx
                or info0.atomIdx(angle.atom2()) in atoms0_idx
            ):
                # Extract the angle information.
                atom0 = info0.atomIdx(angle.atom0())
                atom1 = info0.atomIdx(angle.atom1())
                atom2 = info0.atomIdx(angle.atom2())
                exprn = angle.function()

                # Set the new angle.
                angles.set(atom0, atom1, atom2, exprn)

        # Add the angles to the merged molecule.
        edit_mol.setProperty("angle1", angles)

    # 3) dihedrals
    if "dihedral" in shared_props:
        # Get the info objects for the two molecules.
        info0 = molecule0.info()
        info1 = molecule1.info()

        # Get the user defined property names.
        prop0 = inv_property_map0.get("dihedral", "dihedral")
        prop1 = inv_property_map1.get("dihedral", "dihedral")

        # Get the "dihedral" property from the two molecules.
        dihedrals0 = molecule0.property(prop0)
        dihedrals1 = molecule1.property(prop1)

        # Create the new set of dihedrals.
        dihedrals = _SireMM.FourAtomFunctions(edit_mol.info())

        # Add all of the dihedrals from molecule1.
        for dihedral in dihedrals1.potentials():
            # Extract the dihedral information.
            atom0 = info1.atomIdx(dihedral.atom0())
            atom1 = info1.atomIdx(dihedral.atom1())
            atom2 = info1.atomIdx(dihedral.atom2())
            atom3 = info1.atomIdx(dihedral.atom3())
            exprn = dihedral.function()

            # Map the atom indices to their position in the merged molecule.
            atom0 = inv_mapping[atom0]
            atom1 = inv_mapping[atom1]
            atom2 = inv_mapping[atom2]
            atom3 = inv_mapping[atom3]

            # Set the new dihedral.
            dihedrals.set(atom0, atom1, atom2, atom3, exprn)

        # Loop over all dihedrals in molecule0.
        for dihedral in dihedrals0.potentials():
            # This dihedral contains an atom that is unique to molecule0.
            if (
                info0.atomIdx(dihedral.atom0()) in atoms0_idx
                or info0.atomIdx(dihedral.atom1()) in atoms0_idx
                or info0.atomIdx(dihedral.atom2()) in atoms0_idx
                or info0.atomIdx(dihedral.atom3()) in atoms0_idx
            ):
                # Extract the dihedral information.
                atom0 = info0.atomIdx(dihedral.atom0())
                atom1 = info0.atomIdx(dihedral.atom1())
                atom2 = info0.atomIdx(dihedral.atom2())
                atom3 = info0.atomIdx(dihedral.atom3())
                exprn = dihedral.function()

                # Set the new dihedral.
                dihedrals.set(atom0, atom1, atom2, atom3, exprn)

        # Add the dihedrals to the merged molecule.
        edit_mol.setProperty("dihedral1", dihedrals)

    # 4) impropers
    if "improper" in shared_props:
        # Get the info objects for the two molecules.
        info0 = molecule0.info()
        info1 = molecule1.info()

        # Get the user defined property names.
        prop0 = inv_property_map0.get("improper", "improper")
        prop1 = inv_property_map1.get("improper", "improper")

        # Get the "improper" property from the two molecules.
        impropers0 = molecule0.property(prop0)
        impropers1 = molecule1.property(prop1)

        # Create the new set of impropers.
        impropers = _SireMM.FourAtomFunctions(edit_mol.info())

        # Add all of the impropers from molecule1.
        for improper in impropers1.potentials():
            # Extract the improper information.
            atom0 = info1.atomIdx(improper.atom0())
            atom1 = info1.atomIdx(improper.atom1())
            atom2 = info1.atomIdx(improper.atom2())
            atom3 = info1.atomIdx(improper.atom3())
            exprn = improper.function()

            # Map the atom indices to their position in the merged molecule.
            atom0 = inv_mapping[atom0]
            atom1 = inv_mapping[atom1]
            atom2 = inv_mapping[atom2]
            atom3 = inv_mapping[atom3]

            # Set the new improper.
            impropers.set(atom0, atom1, atom2, atom3, exprn)

        # Loop over all impropers in molecule0.
        for improper in impropers0.potentials():
            # This improper contains an atom that is unique to molecule0.
            if (
                info0.atomIdx(improper.atom0()) in atoms0_idx
                or info0.atomIdx(improper.atom1()) in atoms0_idx
                or info0.atomIdx(improper.atom2()) in atoms0_idx
                or info0.atomIdx(improper.atom3()) in atoms0_idx
            ):
                # Extract the improper information.
                atom0 = info0.atomIdx(improper.atom0())
                atom1 = info0.atomIdx(improper.atom1())
                atom2 = info0.atomIdx(improper.atom2())
                atom3 = info0.atomIdx(improper.atom3())
                exprn = improper.function()

                # Set the new improper.
                impropers.set(atom0, atom1, atom2, atom3, exprn)

        # Add the impropers to the merged molecule.
        edit_mol.setProperty("improper1", impropers)

    # The number of potentials should be consistent for the "bond0"
    # and "bond1" properties, unless a ring is broken or changes size.
    if not (allow_ring_breaking or allow_ring_size_change):
        if (
            edit_mol.property("bond0").nFunctions()
            != edit_mol.property("bond1").nFunctions()
        ):
            raise _IncompatibleError(
                "Inconsistent number of bonds in merged molecule! "
                "A ring may have broken, or changed size. If you want to "
                "allow this perturbation, try using the 'allow_ring_breaking' "
                "or 'allow_ring_size_change' options."
            )

    # Create the connectivity objects.
    conn0 = _SireMol.Connectivity(edit_mol.info()).edit()
    conn1 = _SireMol.Connectivity(edit_mol.info()).edit()

    # Connect the bonded atoms in both end states.
    for bond in edit_mol.property("bond0").potentials():
        conn0.connect(bond.atom0(), bond.atom1())
    conn0 = conn0.commit()
    for bond in edit_mol.property("bond1").potentials():
        conn1.connect(bond.atom0(), bond.atom1())
    conn1 = conn1.commit()

    # Get the original connectivity of the two molecules.
    c0 = molecule0.property("connectivity")
    c1 = molecule1.property("connectivity")

    # Check that the merge hasn't modified the connectivity.

    # molecule0
    for x in range(0, molecule0.nAtoms()):
        # Convert to an AtomIdx.
        idx = _SireMol.AtomIdx(x)

        for y in range(x + 1, molecule0.nAtoms()):
            # Convert to an AtomIdx.
            idy = _SireMol.AtomIdx(y)

            # Was a ring opened/closed?
            is_ring_broken = _is_ring_broken(c0, conn1, idx, idy, idx, idy)

            # A ring was broken and it is not allowed.
            if is_ring_broken and not allow_ring_breaking:
                raise _IncompatibleError(
                    "The merge has opened/closed a ring. To allow this "
                    "perturbation, set the 'allow_ring_breaking' option "
                    "to 'True'."
                )

            # Did a ring change size?
            is_ring_size_change = _is_ring_size_changed(c0, conn1, idx, idy, idx, idy)

            # A ring changed size and it is not allowed.
            if (
                not is_ring_broken
                and is_ring_size_change
                and not allow_ring_size_change
            ):
                raise _IncompatibleError(
                    "The merge has changed the size of a ring. To allow this "
                    "perturbation, set the 'allow_ring_size_change' option "
                    "to 'True'. Be aware that this perturbation may not work "
                    "and a transition through an intermediate state may be "
                    "preferable."
                )

            # The connectivity has changed.
            if c0.connectionType(idx, idy) != conn1.connectionType(idx, idy):
                # The connectivity changed for an unknown reason.
                if not (is_ring_broken or is_ring_size_change) and not force:
                    raise _IncompatibleError(
                        "The merge has changed the molecular connectivity "
                        "but a ring didn't open/close or change size. "
                        "If you want to proceed with this mapping pass "
                        "'force=True'. You are warned that the resulting "
                        "perturbation will likely be unstable."
                    )
    # molecule1
    for x in range(0, molecule1.nAtoms()):
        # Convert to an AtomIdx.
        idx = _SireMol.AtomIdx(x)

        # Map the index to its position in the merged molecule.
        idx_map = inv_mapping[idx]

        for y in range(x + 1, molecule1.nAtoms()):
            # Convert to an AtomIdx.
            idy = _SireMol.AtomIdx(y)

            # Map the index to its position in the merged molecule.
            idy_map = inv_mapping[idy]

            # Was a ring opened/closed?
            is_ring_broken = _is_ring_broken(c1, conn0, idx, idy, idx_map, idy_map)

            # A ring was broken and it is not allowed.
            if is_ring_broken and not allow_ring_breaking:
                raise _IncompatibleError(
                    "The merge has opened/closed a ring. To allow this "
                    "perturbation, set the 'allow_ring_breaking' option "
                    "to 'True'."
                )

            # Did a ring change size?
            is_ring_size_change = _is_ring_size_changed(
                c1, conn0, idx, idy, idx_map, idy_map
            )

            # A ring changed size and it is not allowed.
            if (
                not is_ring_broken
                and is_ring_size_change
                and not allow_ring_size_change
            ):
                raise _IncompatibleError(
                    "The merge has changed the size of a ring. To allow this "
                    "perturbation, set the 'allow_ring_size_change' option "
                    "to 'True'. Be aware that this perturbation may not work "
                    "and a transition through an intermediate state may be "
                    "preferable."
                )

            # The connectivity has changed.
            if c1.connectionType(idx, idy) != conn0.connectionType(idx_map, idy_map):
                # The connectivity changed for an unknown reason.
                if not (is_ring_broken or is_ring_size_change) and not force:
                    raise _IncompatibleError(
                        "The merge has changed the molecular connectivity "
                        "but a ring didn't open/close or change size. "
                        "If you want to proceed with this mapping pass "
                        "'force=True'. You are warned that the resulting "
                        "perturbation will likely be unstable."
                    )

    # Set the "connectivity" property. If the end state connectivity is the same,
    # then we can just set the "connectivity" property.
    if conn0 == conn1:
        edit_mol.setProperty("connectivity", conn0)
    else:
        edit_mol.setProperty("connectivity0", conn0)
        edit_mol.setProperty("connectivity1", conn1)

    # Create the CLJNBPairs matrices.
    ff = molecule0.property(ff0)

    # Initialise the intrascale matrices for both end states.
    clj_nb_pairs0 = _SireMM.CLJNBPairs(edit_mol.info(), _SireMM.CLJScaleFactor(0, 0))
    clj_nb_pairs1 = _SireMM.CLJNBPairs(edit_mol.info(), _SireMM.CLJScaleFactor(0, 0))

    # Loop over all atoms unique to molecule0.
    for idx0 in atoms0_idx:
        # Loop over all atoms unique to molecule1.
        for idx1 in atoms1_idx:
            # Map the index to its position in the merged molecule.
            idx1 = inv_mapping[idx1]

            # Work out the connection type between the atoms in both end states.
            conn_type0 = conn0.connectionType(idx0, idx1)
            conn_type1 = conn1.connectionType(idx0, idx1)

            # Lambda = 0

            # The atoms aren't bonded.
            if conn_type0 == 0:
                clj_scale_factor = _SireMM.CLJScaleFactor(1, 1)
                clj_nb_pairs0.set(idx0, idx1, clj_scale_factor)

            # The atoms are part of a dihedral.
            elif conn_type0 == 4:
                clj_scale_factor = _SireMM.CLJScaleFactor(
                    ff.electrostatic14ScaleFactor(), ff.vdw14ScaleFactor()
                )
                clj_nb_pairs0.set(idx0, idx1, clj_scale_factor)

            # Lambda = 1

            # The atoms aren't bonded.
            if conn_type1 == 0:
                clj_scale_factor = _SireMM.CLJScaleFactor(1, 1)
                clj_nb_pairs1.set(idx0, idx1, clj_scale_factor)

            # The atoms are part of a dihedral.
            elif conn_type1 == 4:
                clj_scale_factor = _SireMM.CLJScaleFactor(
                    ff.electrostatic14ScaleFactor(), ff.vdw14ScaleFactor()
                )
                clj_nb_pairs1.set(idx0, idx1, clj_scale_factor)

    # Get the user defined "intrascale" property names.
    prop0 = inv_property_map0.get("intrascale", "intrascale")
    prop1 = inv_property_map1.get("intrascale", "intrascale")

    # Get the "intrascale" property from the two molecules.
    intrascale0 = molecule0.property(prop0)
    intrascale1 = molecule1.property(prop1)

    # Copy the intrascale from molecule1 into clj_nb_pairs0.

    # Perform a triangular loop over atoms from molecule1.
    for x in range(0, molecule1.nAtoms()):
        # Convert to an AtomIdx.
        idx = _SireMol.AtomIdx(x)

        # Map the index to its position in the merged molecule.
        idx = inv_mapping[idx]

        for y in range(x + 1, molecule1.nAtoms()):
            # Convert to an AtomIdx.
            idy = _SireMol.AtomIdx(y)

            # Map the index to its position in the merged molecule.
            idy = inv_mapping[idy]

            # Get the intrascale value.
            intra = intrascale1.get(_SireMol.AtomIdx(x), _SireMol.AtomIdx(y))

            # Only set if there is a non-zero value.
            # Set using the re-mapped atom indices.
            if not intra.coulomb() == 0:
                clj_nb_pairs0.set(idx, idy, intra)

    # Now copy in all intrascale values from molecule0 into both
    # clj_nb_pairs matrices.

    # Perform a triangular loop over atoms from molecule0.
    for x in range(0, molecule0.nAtoms()):
        for y in range(x + 1, molecule0.nAtoms()):
            # Get the intrascale value.
            intra = intrascale0.get(_SireMol.AtomIdx(x), _SireMol.AtomIdx(y))

            # Set the value in the new matrix, overwriting existing value.
            clj_nb_pairs0.set(_SireMol.AtomIdx(x), _SireMol.AtomIdx(y), intra)

            # Only set if there is a non-zero value.
            if not intra.coulomb() == 0:
                clj_nb_pairs1.set(_SireMol.AtomIdx(x), _SireMol.AtomIdx(y), intra)

    # Finally, copy the intrascale from molecule1 into clj_nb_pairs1.

    # Perform a triangular loop over atoms from molecule1.
    for x in range(0, molecule1.nAtoms()):
        # Convert to an AtomIdx.
        idx = _SireMol.AtomIdx(x)

        # Map the index to its position in the merged molecule.
        idx = inv_mapping[idx]

        for y in range(x + 1, molecule1.nAtoms()):
            # Convert to an AtomIdx.
            idy = _SireMol.AtomIdx(y)

            # Map the index to its position in the merged molecule.
            idy = inv_mapping[idy]

            # Get the intrascale value.
            intra = intrascale1.get(_SireMol.AtomIdx(x), _SireMol.AtomIdx(y))

            # Set the value in the new matrix, overwriting existing value.
            clj_nb_pairs1.set(idx, idy, intra)

    # Store the two molecular components.
    edit_mol.setProperty("molecule0", molecule0)
    edit_mol.setProperty("molecule1", molecule1)

    # Set the "intrascale" properties.
    edit_mol.setProperty("intrascale0", clj_nb_pairs0)
    edit_mol.setProperty("intrascale1", clj_nb_pairs1)

    # Set the "forcefield" properties.
    edit_mol.setProperty("forcefield0", molecule0.property(ff0))
    edit_mol.setProperty("forcefield1", molecule1.property(ff1))

    # Flag that this molecule is perturbable.
    edit_mol.setProperty("is_perturbable", _SireBase.wrap(True))

    # Update the Sire molecule object of the new molecule.
    mol._sire_object = edit_mol.commit()

    # Flag that the molecule is perturbable.
    mol._is_perturbable = True

    # Store the components of the merged molecule.
    mol._molecule0 = _Molecule(molecule0)
    mol._molecule1 = _Molecule(molecule1)

    # Return the new molecule.
    return mol


def _is_ring_broken(conn0, conn1, idx0, idy0, idx1, idy1):
    """
    Internal function to test whether a perturbation changes the connectivity
    around two atoms such that a ring is broken.

    Parameters
    ----------

    conn0 : Sire.Mol.Connectivity
        The connectivity object for the first end state.

    conn1 : Sire.Mol.Connectivity
        The connectivity object for the second end state.

    idx0 : Sire.Mol.AtomIdx
        The index of the first atom in the first state.

    idy0 : Sire.Mol.AtomIdx
        The index of the second atom in the first state.

    idx1 : Sire.Mol.AtomIdx
        The index of the first atom in the second state.

    idy1 : Sire.Mol.AtomIdx
        The index of the second atom in the second state.
    """

    # Have we opened/closed a ring? This means that both atoms are part of a
    # ring in one end state (either in it, or on it), whereas at least one
    # isn't in the other state.

    # Whether each atom is in a ring in both end states.
    in_ring_idx0 = conn0.inRing(idx0)
    in_ring_idy0 = conn0.inRing(idy0)
    in_ring_idx1 = conn1.inRing(idx1)
    in_ring_idy1 = conn1.inRing(idy1)

    # Whether each atom is on a ring in both end states.
    on_ring_idx0 = _is_on_ring(idx0, conn0)
    on_ring_idy0 = _is_on_ring(idy0, conn0)
    on_ring_idx1 = _is_on_ring(idx1, conn1)
    on_ring_idy1 = _is_on_ring(idy1, conn1)

    # Both atoms are in a ring in one end state and at least one isn't in the other.
    if (in_ring_idx0 & in_ring_idy0) ^ (in_ring_idx1 & in_ring_idy1):
        return True

    # Both atoms are on a ring in one end state and at least one isn't in the other.
    if (on_ring_idx0 & on_ring_idy0 & (conn0.connectionType(idx0, idy0) == 4)) ^ (
        on_ring_idx1 & on_ring_idy1 & (conn1.connectionType(idx1, idy1) == 4)
    ):
        # Make sure that the change isn't a result of ring growth, i.e. one of
        # the atoms isn't in a ring in one end state, while its "on" ring status
        # has changed between states.
        if not (
            (in_ring_idx0 | in_ring_idx1) & (on_ring_idx0 ^ on_ring_idx1)
            or (in_ring_idy0 | in_ring_idy1) & (on_ring_idy0 ^ on_ring_idy1)
        ):
            return True

    # Both atoms are in or on a ring in one state and at least one isn't in the other.
    if (
        (in_ring_idx0 | on_ring_idx0)
        & (in_ring_idy0 | on_ring_idy0)
        & (conn0.connectionType(idx0, idy0) == 3)
    ) ^ (
        (in_ring_idx1 | on_ring_idx1)
        & (in_ring_idy1 | on_ring_idy1)
        & (conn1.connectionType(idx1, idy1) == 3)
    ):
        iscn0 = set(conn0.connectionsTo(idx0)).intersection(
            set(conn0.connectionsTo(idy0))
        )
        if len(iscn0) != 1:
            return True
        common_idx = iscn0.pop()
        in_ring_bond0 = conn0.inRing(idx0, common_idx) | conn0.inRing(idy0, common_idx)
        iscn1 = set(conn1.connectionsTo(idx1)).intersection(
            set(conn1.connectionsTo(idy1))
        )
        if len(iscn1) != 1:
            return True
        common_idx = iscn1.pop()
        in_ring_bond1 = conn1.inRing(idx1, common_idx) | conn1.inRing(idy1, common_idx)
        if in_ring_bond0 ^ in_ring_bond1:
            return True

    # If we get this far, then a ring wasn't broken.
    return False


def _is_ring_size_changed(conn0, conn1, idx0, idy0, idx1, idy1, max_ring_size=12):
    """
    Internal function to test whether a perturbation changes the connectivity
    around two atoms such that a ring changes size.

    Parameters
    ----------

    conn0 : Sire.Mol.Connectivity
        The connectivity object for the first end state.

    conn1 : Sire.Mol.Connectivity
        The connectivity object for the second end state.

    idx0 : Sire.Mol.AtomIdx
        The index of the first atom in the first state.

    idy0 : Sire.Mol.AtomIdx
        The index of the second atom in the first state.

    idx1 : Sire.Mol.AtomIdx
        The index of the first atom in the second state.

    idy1 : Sire.Mol.AtomIdx
        The index of the second atom in the second state.

    max_ring_size : int
        The maximum size of what is considered to be a ring.
    """

    # Have a ring changed size? If so, then the minimum path size between
    # two atoms will have changed.

    # Work out the paths connecting the atoms in the two end states.
    paths0 = conn0.findPaths(idx0, idy0, max_ring_size)
    paths1 = conn1.findPaths(idx1, idy1, max_ring_size)

    # Initialise the ring size in each end state.
    ring0 = None
    ring1 = None

    # Determine the minimum path in the lambda = 0 state.
    if len(paths0) > 1:
        path_lengths0 = []
        for path in paths0:
            path_lengths0.append(len(path))
        ring0 = min(path_lengths0)

    if ring0 is None:
        return False

    # Determine the minimum path in the lambda = 1 state.
    if len(paths1) > 1:
        path_lengths1 = []
        for path in paths1:
            path_lengths1.append(len(path))
        ring1 = min(path_lengths1)

    # Return whether the ring has changed size.
    if ring1:
        return ring0 != ring1
    else:
        return False


def _is_on_ring(idx, conn):
    """
    Internal function to test whether an atom is adjacent to a ring.

    Parameters
    ----------

    idx : Sire.Mol.AtomIdx
        The index of the atom

    conn : Sire.Mol.Connectivity
        The connectivity object.

    Returns
    -------

    is_on_ring : bool
        Whether the atom is adjacent to a ring.
    """

    # Loop over all atoms connected to this atom.
    for x in conn.connectionsTo(idx):
        # The neighbour is in a ring.
        if conn.inRing(x) and (not conn.inRing(x, idx)):
            return True

    # If we get this far, then the atom is not adjacent to a ring.
    return False


def _removeDummies(molecule, is_lambda1):
    """
    Internal function which removes the dummy atoms from one of the endstates
    of a merged molecule.

    Parameters
    ----------

    molecule : BioSimSpace._SireWrappers.Molecule
        The molecule.

    is_lambda1 : bool
       Whether to use the molecule at lambda = 1.
    """
    if not molecule._is_perturbable:
        raise _IncompatibleError("'molecule' is not a perturbable molecule")

    # Always use the coordinates at lambda = 0.
    coordinates = molecule._sire_object.property("coordinates0")

    # Generate a molecule with all dummies present.
    molecule = molecule.copy()._toRegularMolecule(
        is_lambda1=is_lambda1, generate_intrascale=True
    )

    # Remove the parameters property, if it exists.
    if "parameters" in molecule._sire_object.propertyKeys():
        molecule._sire_object = (
            molecule._sire_object.edit().removeProperty("parameters").commit()
        )

    # Set the coordinates to those at lambda = 0
    molecule._sire_object = (
        molecule._sire_object.edit().setProperty("coordinates", coordinates).commit()
    )

    # Extract all the nondummy indices
    nondummy_indices = [
        i
        for i, atom in enumerate(molecule.getAtoms())
        if "du" not in atom._sire_object.property("ambertype")
    ]

    # Create an AtomSelection.
    selection = molecule._sire_object.selection()

    # Unselect all of the atoms.
    selection.selectNone()

    # Now add all of the nondummy atoms.
    for idx in nondummy_indices:
        selection.select(_SireMol.AtomIdx(idx))

    # Create a partial molecule and extract the atoms.
    partial_molecule = (
        _SireMol.PartialMolecule(molecule._sire_object, selection).extract().molecule()
    )

    # Remove the incorrect intrascale property.
    partial_molecule = (
        partial_molecule.edit().removeProperty("intrascale").molecule().commit()
    )

    # Recreate a BioSimSpace molecule object.
    molecule = _Molecule(partial_molecule)

    # Parse the molecule as a GROMACS topology, which will recover the intrascale
    # matrix.
    gro_top = _SireIO.GroTop(molecule.toSystem()._sire_object)

    # Convert back to a Sire system.
    gro_sys = gro_top.toSystem()

    # Add the intrascale property back into the merged molecule.
    edit_mol = molecule._sire_object.edit()
    edit_mol = edit_mol.setProperty(
        "intrascale", gro_sys[_SireMol.MolIdx(0)].property("intrascale")
    )
    molecule = _Molecule(edit_mol.commit())

    return molecule
