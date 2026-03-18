######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2025
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


def merge(
    molecule0,
    molecule1,
    mapping,
    allow_ring_breaking=False,
    allow_ring_size_change=False,
    force=False,
    roi=None,
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

    roi : list
                        The region of interest to merge. Consist of two lists of atom indices.

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
    from sire.legacy import IO as _SireIO
    from sire.legacy import MM as _SireMM
    from sire.legacy import Base as _SireBase
    from sire.legacy import Mol as _SireMol
    from sire.legacy import Units as _SireUnits

    from .._Exceptions import IncompatibleError as _IncompatibleError
    from .._SireWrappers import Molecule as _Molecule

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

    # Generate the mappings from each molecule to the merged molecule
    mol0_merged_mapping = {}
    mol1_merged_mapping = {}

    # Invert the user property mappings.
    inv_property_map0 = {v: k for k, v in property_map0.items()}
    inv_property_map1 = {v: k for k, v in property_map1.items()}

    # Make sure that the molecules have a "forcefield" property and that
    # the two force fields are compatible.

    # Get the user name for the "forcefield" property.
    ff0 = inv_property_map0.get("forcefield", "forcefield")
    ff1 = inv_property_map1.get("forcefield", "forcefield")

    # Force field information is missing.
    if not molecule0.has_property(ff0):
        raise _IncompatibleError("Cannot determine 'forcefield' of 'molecule0'!")
    if not molecule1.has_property(ff1):
        raise _IncompatibleError("Cannot determine 'forcefield' of 'molecule1'!")

    # The force fields are incompatible.
    if not molecule0.property(ff0).is_compatible_with(molecule1.property(ff1)):
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
    # Only part of the ligand is to be merged
    if roi is not None:
        if molecule0.num_residues() != molecule1.num_residues():
            raise ValueError(
                "The two molecules need to have the same number of residues"
            )

        num = 1
        for idx, (mol0_res, mol1_res) in enumerate(
            zip(molecule0.residues(), molecule1.residues())
        ):
            res = molecule.edit().add(_SireMol.ResNum(idx + 1))
            if mol0_res.name() == mol1_res.name():
                res.rename(mol0_res.name())
            else:
                resname = "MUT" if len(molecule0.residues()) > 1 else "LIG"
                res.rename(_SireMol.ResName(resname))

            cg = res.molecule().add(_SireMol.CGName(f"{idx}"))
            for atom in mol0_res.atoms():
                mol0_merged_mapping[atom.index()] = _SireMol.AtomIdx(num - 1)
                added = cg.add(atom.name())
                added.renumber(_SireMol.AtomNum(num))
                added.reparent(_SireMol.ResIdx(idx))
                num += 1

            for atom in mol1_res.atoms():
                if atom.index() in atoms1_idx:
                    added = cg.add(atom.name())
                    added.renumber(_SireMol.AtomNum(num))
                    added.reparent(_SireMol.ResIdx(idx))
                    mol1_merged_mapping[atom.index()] = _SireMol.AtomIdx(num - 1)
                    num += 1
                else:
                    mol1_merged_mapping[atom.index()] = mol0_merged_mapping[
                        inv_mapping[atom.index()]
                    ]
            molecule = cg.molecule().commit()
    else:
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
            mol0_merged_mapping[atom.index()] = _SireMol.AtomIdx(num - 1)
            num += 1

        # Now add all of the atoms from molecule1 that aren't mapped from molecule0.
        for atom in molecule1.atoms():
            if atom in atoms1:
                added = cg.add(atom.name())
                added.renumber(_SireMol.AtomNum(num))
                added.reparent(_SireMol.ResIdx(0))
                mol1_merged_mapping[atom.index()] = _SireMol.AtomIdx(num - 1)
                num += 1
            else:
                mol1_merged_mapping[atom.index()] = mol0_merged_mapping[
                    inv_mapping[atom.index()]
                ]

        # Commit the changes to the molecule.
        molecule = cg.molecule().commit()

    # Make the molecule editable.
    edit_mol = molecule.edit()

    # Create lists of the actual property names in the molecules.
    props0 = []
    props1 = []

    # molecule0
    for prop in molecule0.property_keys():
        if prop in inv_property_map0:
            prop = inv_property_map0[prop]
        props0.append(prop)

    # molecule1
    for prop in molecule1.property_keys():
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
        "cmap",
        "connectivity",
        "intrascale",
    ]

    # lambda = 0
    for prop in props0:
        if prop not in ignored_props:
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
                edit_mol = edit_mol.set_property(name, default_prop).molecule()
            except:
                pass

    # lambda = 1
    for prop in props1:
        if prop not in ignored_props:
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
                edit_mol = edit_mol.set_property(name, default_prop).molecule()
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
        # Get the atom index in the merged molecule.
        idx = mol0_merged_mapping[atom.index()]

        # Add an "name0" property.
        edit_mol = (
            edit_mol.atom(idx).set_property("name0", atom.name().value()).molecule()
        )

        # Loop over all atom properties.
        for prop in atom.property_keys():
            # Get the actual property name.
            name = inv_property_map0.get(prop, prop)

            # This is a perturbable property. Rename to "property0", e.g. "charge0".
            if name in shared_props:
                name = name + "0"

            # Add the property to the atom in the merged molecule.
            edit_mol = (
                edit_mol.atom(idx).set_property(name, atom.property(prop)).molecule()
            )

    # Add the atom properties from molecule1.
    for atom in atoms1:
        # Get the atom index in the merged molecule.
        idx = mol1_merged_mapping[atom.index()]

        # Add an "name0" property.
        edit_mol = (
            edit_mol.atom(idx).set_property("name0", atom.name().value()).molecule()
        )

        # Loop over all atom properties.
        for prop in atom.property_keys():
            # Get the actual property name.
            name = inv_property_map1.get(prop, prop)

            # Zero the "charge" and "LJ" property for atoms that are unique to molecule1.
            if name == "charge":
                edit_mol = (
                    edit_mol.atom(idx)
                    .set_property("charge0", 0 * _SireUnits.e_charge)
                    .molecule()
                )
            elif name == "LJ":
                edit_mol = (
                    edit_mol.atom(idx)
                    .set_property("LJ0", _SireMM.LJParameter())
                    .molecule()
                )
            elif name == "ambertype":
                edit_mol = (
                    edit_mol.atom(idx).set_property("ambertype0", "du").molecule()
                )
            elif name == "element":
                edit_mol = (
                    edit_mol.atom(idx)
                    .set_property("element0", _SireMol.Element(0))
                    .molecule()
                )
            else:
                # This is a perturbable property. Rename to "property0", e.g. "charge0".
                if name in shared_props:
                    name = name + "0"

                # Add the property to the atom in the merged molecule.
                edit_mol = (
                    edit_mol.atom(idx)
                    .set_property(name, atom.property(prop))
                    .molecule()
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
            atom0 = mol0_merged_mapping[info0.atom_idx(bond.atom0())]
            atom1 = mol0_merged_mapping[info0.atom_idx(bond.atom1())]
            bonds.set(atom0, atom1, bond.function())

        # Loop over all bonds in molecule1.
        for bond in bonds1.potentials():
            # This bond contains an atom that is unique to molecule1.
            if (
                info1.atom_idx(bond.atom0()) in atoms1_idx
                or info1.atom_idx(bond.atom1()) in atoms1_idx
            ):
                # Extract the bond information.
                atom0 = info1.atom_idx(bond.atom0())
                atom1 = info1.atom_idx(bond.atom1())
                exprn = bond.function()

                # Map the atom indices to their position in the merged molecule.
                atom0 = mol1_merged_mapping[atom0]
                atom1 = mol1_merged_mapping[atom1]

                # Set the new bond.
                bonds.set(atom0, atom1, exprn)

        # Add the bonds to the merged molecule.
        edit_mol.set_property("bond0", bonds)

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
            atom0 = mol0_merged_mapping[info0.atom_idx(angle.atom0())]
            atom1 = mol0_merged_mapping[info0.atom_idx(angle.atom1())]
            atom2 = mol0_merged_mapping[info0.atom_idx(angle.atom2())]
            angles.set(atom0, atom1, atom2, angle.function())

        # Loop over all angles in molecule1.
        for angle in angles1.potentials():
            # This angle contains an atom that is unique to molecule1.
            if (
                info1.atom_idx(angle.atom0()) in atoms1_idx
                or info1.atom_idx(angle.atom1()) in atoms1_idx
                or info1.atom_idx(angle.atom2()) in atoms1_idx
            ):
                # Extract the angle information.
                atom0 = info1.atom_idx(angle.atom0())
                atom1 = info1.atom_idx(angle.atom1())
                atom2 = info1.atom_idx(angle.atom2())
                exprn = angle.function()

                # Map the atom indices to their position in the merged molecule.
                atom0 = mol1_merged_mapping[atom0]
                atom1 = mol1_merged_mapping[atom1]
                atom2 = mol1_merged_mapping[atom2]

                # Set the new angle.
                angles.set(atom0, atom1, atom2, exprn)

        # Add the angles to the merged molecule.
        edit_mol.set_property("angle0", angles)

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
            atom0 = mol0_merged_mapping[info0.atom_idx(dihedral.atom0())]
            atom1 = mol0_merged_mapping[info0.atom_idx(dihedral.atom1())]
            atom2 = mol0_merged_mapping[info0.atom_idx(dihedral.atom2())]
            atom3 = mol0_merged_mapping[info0.atom_idx(dihedral.atom3())]
            dihedrals.set(atom0, atom1, atom2, atom3, dihedral.function())

        # Loop over all dihedrals in molecule1.
        for dihedral in dihedrals1.potentials():
            # This dihedral contains an atom that is unique to molecule1.
            if (
                info1.atom_idx(dihedral.atom0()) in atoms1_idx
                or info1.atom_idx(dihedral.atom1()) in atoms1_idx
                or info1.atom_idx(dihedral.atom2()) in atoms1_idx
                or info1.atom_idx(dihedral.atom3()) in atoms1_idx
            ):
                # Extract the dihedral information.
                atom0 = info1.atom_idx(dihedral.atom0())
                atom1 = info1.atom_idx(dihedral.atom1())
                atom2 = info1.atom_idx(dihedral.atom2())
                atom3 = info1.atom_idx(dihedral.atom3())
                exprn = dihedral.function()

                # Map the atom indices to their position in the merged molecule.
                atom0 = mol1_merged_mapping[atom0]
                atom1 = mol1_merged_mapping[atom1]
                atom2 = mol1_merged_mapping[atom2]
                atom3 = mol1_merged_mapping[atom3]

                # Set the new dihedral.
                dihedrals.set(atom0, atom1, atom2, atom3, exprn)

        # Add the dihedrals to the merged molecule.
        edit_mol.set_property("dihedral0", dihedrals)

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
            atom0 = mol0_merged_mapping[info0.atom_idx(improper.atom0())]
            atom1 = mol0_merged_mapping[info0.atom_idx(improper.atom1())]
            atom2 = mol0_merged_mapping[info0.atom_idx(improper.atom2())]
            atom3 = mol0_merged_mapping[info0.atom_idx(improper.atom3())]
            impropers.set(atom0, atom1, atom2, atom3, improper.function())

        # Loop over all impropers in molecule1.
        for improper in impropers1.potentials():
            # This improper contains an atom that is unique to molecule1.
            if (
                info1.atom_idx(improper.atom0()) in atoms1_idx
                or info1.atom_idx(improper.atom1()) in atoms1_idx
                or info1.atom_idx(improper.atom2()) in atoms1_idx
                or info1.atom_idx(improper.atom3()) in atoms1_idx
            ):
                # Extract the improper information.
                atom0 = info1.atom_idx(improper.atom0())
                atom1 = info1.atom_idx(improper.atom1())
                atom2 = info1.atom_idx(improper.atom2())
                atom3 = info1.atom_idx(improper.atom3())
                exprn = improper.function()

                # Map the atom indices to their position in the merged molecule.
                atom0 = mol1_merged_mapping[atom0]
                atom1 = mol1_merged_mapping[atom1]
                atom2 = mol1_merged_mapping[atom2]
                atom3 = mol1_merged_mapping[atom3]

                # Set the new improper.
                impropers.set(atom0, atom1, atom2, atom3, exprn)

        # Add the impropers to the merged molecule.
        edit_mol.set_property("improper0", impropers)

    # 5) cmap
    if "cmap" in shared_props:
        # Get the user defined property names.
        prop0 = inv_property_map0.get("cmap", "cmap")
        prop1 = inv_property_map1.get("cmap", "cmap")

        # Get the "cmap" property from the two molecules.
        cmaps0 = molecule0.property(prop0)
        cmaps1 = molecule1.property(prop1)

        # Get the molInfo object for each molecule.
        info0 = molecule0.info()
        info1 = molecule1.info()

        # Create the new set of CMAP functions.
        cmaps = _SireMM.CMAPFunctions(edit_mol.info())

        # Add all of the CMAP functions from molecule0.
        for func in cmaps0.parameters():
            atom0 = mol0_merged_mapping[info0.atom_idx(func.atom0())]
            atom1 = mol0_merged_mapping[info0.atom_idx(func.atom1())]
            atom2 = mol0_merged_mapping[info0.atom_idx(func.atom2())]
            atom3 = mol0_merged_mapping[info0.atom_idx(func.atom3())]
            atom4 = mol0_merged_mapping[info0.atom_idx(func.atom4())]
            cmaps.set(atom0, atom1, atom2, atom3, atom4, func.parameter())

        # Loop over all CMAP functions in molecule1.
        for func in cmaps1.parameters():
            # This CMAP function contains an atom that is unique to molecule1.
            if (
                info1.atom_idx(func.atom0()) in atoms1_idx
                or info1.atom_idx(func.atom1()) in atoms1_idx
                or info1.atom_idx(func.atom2()) in atoms1_idx
                or info1.atom_idx(func.atom3()) in atoms1_idx
                or info1.atom_idx(func.atom4()) in atoms1_idx
            ):
                atom0 = mol1_merged_mapping[info1.atom_idx(func.atom0())]
                atom1 = mol1_merged_mapping[info1.atom_idx(func.atom1())]
                atom2 = mol1_merged_mapping[info1.atom_idx(func.atom2())]
                atom3 = mol1_merged_mapping[info1.atom_idx(func.atom3())]
                atom4 = mol1_merged_mapping[info1.atom_idx(func.atom4())]
                cmaps.set(atom0, atom1, atom2, atom3, atom4, func.parameter())

        # Add the CMAP functions to the merged molecule.
        edit_mol.set_property("cmap0", cmaps)

    ##############################
    # SET PROPERTIES AT LAMBDA = 1
    ##############################

    # Add the atom properties from molecule1.
    for atom in molecule1.atoms():
        # Get the atom index in the merged molecule.
        idx = mol1_merged_mapping[atom.index()]

        # Add an "name1" property.
        edit_mol = (
            edit_mol.atom(idx).set_property("name1", atom.name().value()).molecule()
        )

        # Loop over all atom properties.
        for prop in atom.property_keys():
            # Get the actual property name.
            name = inv_property_map1.get(prop, prop)

            # This is a perturbable property. Rename to "property1", e.g. "charge1".
            if name in shared_props:
                name = name + "1"

            # Add the property to the atom in the merged molecule.
            edit_mol = (
                edit_mol.atom(idx).set_property(name, atom.property(prop)).molecule()
            )

    # Add the properties from atoms unique to molecule0.
    for atom in atoms0:
        # Get the atom index in the merged molecule.
        idx = mol0_merged_mapping[atom.index()]

        # Add an "name1" property.
        edit_mol = (
            edit_mol.atom(idx).set_property("name1", atom.name().value()).molecule()
        )

        # Loop over all atom properties.
        for prop in atom.property_keys():
            # Get the actual property name.
            name = inv_property_map0.get(prop, prop)

            # Zero the "charge" and "LJ" property for atoms that are unique to molecule0.
            if name == "charge":
                edit_mol = (
                    edit_mol.atom(idx)
                    .set_property("charge1", 0 * _SireUnits.e_charge)
                    .molecule()
                )
            elif name == "LJ":
                edit_mol = (
                    edit_mol.atom(idx)
                    .set_property("LJ1", _SireMM.LJParameter())
                    .molecule()
                )
            elif name == "ambertype":
                edit_mol = (
                    edit_mol.atom(idx).set_property("ambertype1", "du").molecule()
                )
            elif name == "element":
                edit_mol = (
                    edit_mol.atom(idx)
                    .set_property("element1", _SireMol.Element(0))
                    .molecule()
                )
            else:
                # This is a perturbable property. Rename to "property1", e.g. "charge1".
                if name in shared_props:
                    name = name + "1"

                # Add the property to the atom in the merged molecule.
                edit_mol = (
                    edit_mol.atom(idx)
                    .set_property(name, atom.property(prop))
                    .molecule()
                )

    # Tolerance for zero sigma values.
    null_lj_sigma = 1e-9

    # Atoms with zero LJ sigma values need to have their sigma values set to the
    # value from the other end state.
    for atom in edit_mol.atoms():
        # Get the end state LJ sigma values.
        lj0 = atom.property("LJ0")
        lj1 = atom.property("LJ1")

        # Lambda = 0 state has a zero sigma value.
        if abs(lj0.sigma().value()) <= null_lj_sigma:
            # Use the sigma value from the lambda = 1 state.
            edit_mol = (
                edit_mol.atom(atom.index())
                .set_property("LJ0", _SireMM.LJParameter(lj1.sigma(), lj0.epsilon()))
                .molecule()
            )

        # Lambda = 1 state has a zero sigma value.
        if abs(lj1.sigma().value()) <= null_lj_sigma:
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
            atom0 = info1.atom_idx(bond.atom0())
            atom1 = info1.atom_idx(bond.atom1())
            exprn = bond.function()

            # Map the atom indices to their position in the merged molecule.
            atom0 = mol1_merged_mapping[atom0]
            atom1 = mol1_merged_mapping[atom1]

            # Set the new bond.
            bonds.set(atom0, atom1, exprn)

        # Loop over all bonds in molecule0
        for bond in bonds0.potentials():
            # This bond contains an atom that is unique to molecule0.
            if (
                info0.atom_idx(bond.atom0()) in atoms0_idx
                or info0.atom_idx(bond.atom1()) in atoms0_idx
            ):
                # Extract the bond information.
                atom0 = mol0_merged_mapping[info0.atom_idx(bond.atom0())]
                atom1 = mol0_merged_mapping[info0.atom_idx(bond.atom1())]
                exprn = bond.function()

                # Set the new bond.
                bonds.set(atom0, atom1, exprn)

        # Add the bonds to the merged molecule.
        edit_mol.set_property("bond1", bonds)

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
            atom0 = info1.atom_idx(angle.atom0())
            atom1 = info1.atom_idx(angle.atom1())
            atom2 = info1.atom_idx(angle.atom2())
            exprn = angle.function()

            # Map the atom indices to their position in the merged molecule.
            atom0 = mol1_merged_mapping[atom0]
            atom1 = mol1_merged_mapping[atom1]
            atom2 = mol1_merged_mapping[atom2]

            # Set the new angle.
            angles.set(atom0, atom1, atom2, exprn)

        # Loop over all angles in molecule0.
        for angle in angles0.potentials():
            # This angle contains an atom that is unique to molecule0.
            if (
                info0.atom_idx(angle.atom0()) in atoms0_idx
                or info0.atom_idx(angle.atom1()) in atoms0_idx
                or info0.atom_idx(angle.atom2()) in atoms0_idx
            ):
                # Extract the angle information.
                atom0 = mol0_merged_mapping[info0.atom_idx(angle.atom0())]
                atom1 = mol0_merged_mapping[info0.atom_idx(angle.atom1())]
                atom2 = mol0_merged_mapping[info0.atom_idx(angle.atom2())]
                exprn = angle.function()

                # Set the new angle.
                angles.set(atom0, atom1, atom2, exprn)

        # Add the angles to the merged molecule.
        edit_mol.set_property("angle1", angles)

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
            atom0 = info1.atom_idx(dihedral.atom0())
            atom1 = info1.atom_idx(dihedral.atom1())
            atom2 = info1.atom_idx(dihedral.atom2())
            atom3 = info1.atom_idx(dihedral.atom3())
            exprn = dihedral.function()

            # Map the atom indices to their position in the merged molecule.
            atom0 = mol1_merged_mapping[atom0]
            atom1 = mol1_merged_mapping[atom1]
            atom2 = mol1_merged_mapping[atom2]
            atom3 = mol1_merged_mapping[atom3]

            # Set the new dihedral.
            dihedrals.set(atom0, atom1, atom2, atom3, exprn)

        # Loop over all dihedrals in molecule0.
        for dihedral in dihedrals0.potentials():
            # This dihedral contains an atom that is unique to molecule0.
            if (
                info0.atom_idx(dihedral.atom0()) in atoms0_idx
                or info0.atom_idx(dihedral.atom1()) in atoms0_idx
                or info0.atom_idx(dihedral.atom2()) in atoms0_idx
                or info0.atom_idx(dihedral.atom3()) in atoms0_idx
            ):
                # Extract the dihedral information.
                atom0 = mol0_merged_mapping[info0.atom_idx(dihedral.atom0())]
                atom1 = mol0_merged_mapping[info0.atom_idx(dihedral.atom1())]
                atom2 = mol0_merged_mapping[info0.atom_idx(dihedral.atom2())]
                atom3 = mol0_merged_mapping[info0.atom_idx(dihedral.atom3())]
                exprn = dihedral.function()

                # Set the new dihedral.
                dihedrals.set(atom0, atom1, atom2, atom3, exprn)

        # Add the dihedrals to the merged molecule.
        edit_mol.set_property("dihedral1", dihedrals)

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
            atom0 = info1.atom_idx(improper.atom0())
            atom1 = info1.atom_idx(improper.atom1())
            atom2 = info1.atom_idx(improper.atom2())
            atom3 = info1.atom_idx(improper.atom3())
            exprn = improper.function()

            # Map the atom indices to their position in the merged molecule.
            atom0 = mol1_merged_mapping[atom0]
            atom1 = mol1_merged_mapping[atom1]
            atom2 = mol1_merged_mapping[atom2]
            atom3 = mol1_merged_mapping[atom3]

            # Set the new improper.
            impropers.set(atom0, atom1, atom2, atom3, exprn)

        # Loop over all impropers in molecule0.
        for improper in impropers0.potentials():
            # This improper contains an atom that is unique to molecule0.
            if (
                info0.atom_idx(improper.atom0()) in atoms0_idx
                or info0.atom_idx(improper.atom1()) in atoms0_idx
                or info0.atom_idx(improper.atom2()) in atoms0_idx
                or info0.atom_idx(improper.atom3()) in atoms0_idx
            ):
                # Extract the improper information.
                atom0 = mol0_merged_mapping[info0.atom_idx(improper.atom0())]
                atom1 = mol0_merged_mapping[info0.atom_idx(improper.atom1())]
                atom2 = mol0_merged_mapping[info0.atom_idx(improper.atom2())]
                atom3 = mol0_merged_mapping[info0.atom_idx(improper.atom3())]
                exprn = improper.function()

                # Set the new improper.
                impropers.set(atom0, atom1, atom2, atom3, exprn)

        # Add the impropers to the merged molecule.
        edit_mol.set_property("improper1", impropers)

    # 5) cmap
    if "cmap" in shared_props:
        # Get the user defined property names.
        prop0 = inv_property_map0.get("cmap", "cmap")
        prop1 = inv_property_map1.get("cmap", "cmap")

        # Get the "cmap" property from the two molecules.
        cmaps0 = molecule0.property(prop0)
        cmaps1 = molecule1.property(prop1)

        # Get the molInfo object for each molecule.
        info0 = molecule0.info()
        info1 = molecule1.info()

        # Create the new set of CMAP functions.
        cmaps = _SireMM.CMAPFunctions(edit_mol.info())

        # Add all of the CMAP functions from molecule1.
        for func in cmaps1.parameters():
            atom0 = mol1_merged_mapping[info1.atom_idx(func.atom0())]
            atom1 = mol1_merged_mapping[info1.atom_idx(func.atom1())]
            atom2 = mol1_merged_mapping[info1.atom_idx(func.atom2())]
            atom3 = mol1_merged_mapping[info1.atom_idx(func.atom3())]
            atom4 = mol1_merged_mapping[info1.atom_idx(func.atom4())]
            cmaps.set(atom0, atom1, atom2, atom3, atom4, func.parameter())

        # Loop over all CMAP functions in molecule0.
        for func in cmaps0.parameters():
            # This CMAP function contains an atom that is unique to molecule0.
            if (
                info0.atom_idx(func.atom0()) in atoms0_idx
                or info0.atom_idx(func.atom1()) in atoms0_idx
                or info0.atom_idx(func.atom2()) in atoms0_idx
                or info0.atom_idx(func.atom3()) in atoms0_idx
                or info0.atom_idx(func.atom4()) in atoms0_idx
            ):
                atom0 = mol0_merged_mapping[info0.atom_idx(func.atom0())]
                atom1 = mol0_merged_mapping[info0.atom_idx(func.atom1())]
                atom2 = mol0_merged_mapping[info0.atom_idx(func.atom2())]
                atom3 = mol0_merged_mapping[info0.atom_idx(func.atom3())]
                atom4 = mol0_merged_mapping[info0.atom_idx(func.atom4())]
                cmaps.set(atom0, atom1, atom2, atom3, atom4, func.parameter())

        # Add the CMAP functions to the merged molecule.
        edit_mol.set_property("cmap1", cmaps)

    # The number of potentials should be consistent for the "bond0"
    # and "bond1" properties, unless a ring is broken or changes size.
    if not (allow_ring_breaking or allow_ring_size_change):
        if (
            edit_mol.property("bond0").num_functions()
            != edit_mol.property("bond1").num_functions()
        ):
            raise _IncompatibleError(
                "Inconsistent number of bonds in merged molecule! "
                "A ring may have broken, or changed size. If you want to "
                "allow this perturbation, try using the 'allow_ring_breaking' "
                "or 'allow_ring_size_change' options."
            )

    # Create the connectivity object.
    conn = _SireMol.Connectivity(edit_mol.info()).edit()

    # Connectivity in the merged molecule.
    conn0 = _SireMol.Connectivity(edit_mol.info()).edit()
    conn1 = _SireMol.Connectivity(edit_mol.info()).edit()

    bonds_atoms0 = {
        frozenset([bond.atom0(), bond.atom1()])
        for bond in edit_mol.property("bond0").potentials()
    }
    bonds_atoms1 = {
        frozenset([bond.atom0(), bond.atom1()])
        for bond in edit_mol.property("bond1").potentials()
    }

    for atom0, atom1 in bonds_atoms0:
        conn0.connect(atom0, atom1)

    for atom0, atom1 in bonds_atoms1:
        conn1.connect(atom0, atom1)

    # We only add a bond to the total connectivity if it's defined in both states
    # This results in a "broken" topology if one writes it in GROMACS
    # But GROMACS can't handle bond breaks anyway
    for atom0, atom1 in bonds_atoms0 & bonds_atoms1:
        conn.connect(atom0, atom1)

    conn = conn.commit()
    conn0 = conn0.commit()
    conn1 = conn1.commit()

    # Get the connectivity of the two molecules.
    c0 = molecule0.property("connectivity")
    c1 = molecule1.property("connectivity")

    # Check that the merge hasn't modified the connectivity.

    # The checking was blocked when merging a protein
    if roi is None:
        # Pre-compute ring membership as frozensets of integer indices for O(1)
        # lookup in the inner loops, avoiding repeated Sire in_ring calls.
        n_merged = edit_mol.info().num_atoms()
        c0_ring = frozenset(
            x for x in range(molecule0.num_atoms()) if c0.in_ring(_SireMol.AtomIdx(x))
        )
        c1_ring = frozenset(
            x for x in range(molecule1.num_atoms()) if c1.in_ring(_SireMol.AtomIdx(x))
        )
        conn_ring = frozenset(
            i for i in range(n_merged) if conn.in_ring(_SireMol.AtomIdx(i))
        )

        # Pre-compute merged-space AtomIdx arrays so the inner loops avoid
        # repeated dict lookups and AtomIdx construction.
        mol0_map = [
            mol0_merged_mapping[_SireMol.AtomIdx(i)]
            for i in range(molecule0.num_atoms())
        ]
        mol1_map = [
            mol1_merged_mapping[_SireMol.AtomIdx(i)]
            for i in range(molecule1.num_atoms())
        ]

        # molecule0
        for x in range(molecule0.num_atoms()):
            idx = _SireMol.AtomIdx(x)
            idx_map = mol0_map[x]
            x_in_ring = x in c0_ring or idx_map.value() in conn_ring

            for y in range(x + 1, molecule0.num_atoms()):
                idy = _SireMol.AtomIdx(y)
                idy_map = mol0_map[y]

                # Check whether the connectivity type has changed.
                conn_type_changed = c0.connection_type(
                    idx, idy
                ) != conn.connection_type(idx_map, idy_map)

                # Fast path: if neither atom is in a ring in either end state
                # and the connectivity hasn't changed, nothing to check.
                if (
                    not x_in_ring
                    and y not in c0_ring
                    and idy_map.value() not in conn_ring
                    and not conn_type_changed
                ):
                    continue

                # Combined ring check — calls find_paths once per connectivity.
                is_ring_broken, is_ring_size_change = _check_ring(
                    c0, conn, idx, idy, idx_map, idy_map
                )

                # A ring was broken and it is not allowed.
                if is_ring_broken and not allow_ring_breaking:
                    raise _IncompatibleError(
                        "The merge has opened/closed a ring. To allow this "
                        "perturbation, set the 'allow_ring_breaking' option "
                        "to 'True'."
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

                # The connectivity changed for an unknown reason.
                if (
                    conn_type_changed
                    and not (is_ring_broken or is_ring_size_change)
                    and not force
                ):
                    raise _IncompatibleError(
                        "The merge has changed the molecular connectivity "
                        "but a ring didn't open/close or change size. "
                        "If you want to proceed with this mapping pass "
                        "'force=True'. You are warned that the resulting "
                        "perturbation will likely be unstable."
                    )

        # molecule1
        for x in range(molecule1.num_atoms()):
            idx = _SireMol.AtomIdx(x)
            idx_map = mol1_map[x]
            x_in_ring = x in c1_ring or idx_map.value() in conn_ring

            for y in range(x + 1, molecule1.num_atoms()):
                idy = _SireMol.AtomIdx(y)
                idy_map = mol1_map[y]

                # Check whether the connectivity type has changed.
                conn_type_changed = c1.connection_type(
                    idx, idy
                ) != conn.connection_type(idx_map, idy_map)

                # Fast path: if neither atom is in a ring in either end state
                # and the connectivity hasn't changed, nothing to check.
                if (
                    not x_in_ring
                    and y not in c1_ring
                    and idy_map.value() not in conn_ring
                    and not conn_type_changed
                ):
                    continue

                # Combined ring check — calls find_paths once per connectivity.
                is_ring_broken, is_ring_size_change = _check_ring(
                    c1, conn, idx, idy, idx_map, idy_map
                )

                # A ring was broken and it is not allowed.
                if is_ring_broken and not allow_ring_breaking:
                    raise _IncompatibleError(
                        "The merge has opened/closed a ring. To allow this "
                        "perturbation, set the 'allow_ring_breaking' option "
                        "to 'True'."
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

                # The connectivity changed for an unknown reason.
                if (
                    conn_type_changed
                    and not (is_ring_broken or is_ring_size_change)
                    and not force
                ):
                    raise _IncompatibleError(
                        "The merge has changed the molecular connectivity "
                        "but a ring didn't open/close or change size. "
                        "If you want to proceed with this mapping pass "
                        "'force=True'. You are warned that the resulting "
                        "perturbation will likely be unstable."
                    )

    # Set the "connectivity" property.
    edit_mol.set_property("connectivity", conn)

    # Merge the intrascale properties of the two molecules.
    merged_intrascale = _SireIO.mergeIntrascale(
        molecule0.property("intrascale"),
        molecule1.property("intrascale"),
        edit_mol.info(),
        mol0_merged_mapping,
        mol1_merged_mapping,
    )

    # Store the two molecular components.
    edit_mol.set_property("molecule0", molecule0)
    edit_mol.set_property("molecule1", molecule1)

    # Set the "intrascale" properties.
    edit_mol.set_property("intrascale0", merged_intrascale[0])
    edit_mol.set_property("intrascale1", merged_intrascale[1])

    # Set the "forcefield" properties.
    edit_mol.set_property("forcefield0", molecule0.property(ff0))
    edit_mol.set_property("forcefield1", molecule1.property(ff1))

    # Flag that this molecule is perturbable.
    edit_mol.set_property("is_perturbable", _SireBase.wrap(True))

    # Update the Sire molecule object of the new molecule.
    mol._sire_object = edit_mol.commit()

    # Flag that the molecule is perturbable.
    mol._is_perturbable = True

    # Store the components of the merged molecule.
    mol._molecule0 = _Molecule(molecule0)
    mol._molecule1 = _Molecule(molecule1)

    # Return the new molecule.
    return mol


def _check_ring(conn0, conn1, idx0, idy0, idx1, idy1, max_path=50, max_ring_size=12):
    """
    Internal function to test whether a perturbation opens/closes a ring or
    changes its size for a given pair of atoms.

    Combines the work of the former _is_ring_broken and _is_ring_size_changed
    into a single function, calling find_paths only once per connectivity
    object rather than twice.

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

    max_path : int
        Maximum path length used when searching for rings. The default of 50
        covers typical macrocycles. Increase if larger rings need to be
        detected.

    max_ring_size : int
        The maximum ring size considered when checking for ring size changes.

    Returns
    -------

    (is_ring_broken, is_ring_size_changed) : (bool, bool)
    """

    # Find all paths between the two atoms in each end state. A single
    # find_paths call covers both the ring-broken and ring-size-changed checks.
    paths0 = conn0.find_paths(idx0, idy0, max_path)
    paths1 = conn1.find_paths(idx1, idy1, max_path)
    n0 = len(paths0)
    n1 = len(paths1)

    # --- Ring broken check ---

    # Two atoms share a ring iff there are ≥2 distinct paths between them.
    # This also handles atoms adjacent to a ring (not members themselves):
    # substituents on neighbouring ring atoms have ≥2 paths between them
    # through the ring, which collapses to 1 when the ring is broken.
    if (n0 >= 2) != (n1 >= 2):
        return True, False

    # Supplementary check for rings larger than max_path: find_paths may only
    # find the direct-bond path and miss the long way around the ring, giving
    # n=1 instead of n≥2. Sire's in_ring has no path-length limit and
    # correctly identifies ring membership in macrocycles.
    if (conn0.in_ring(idx0) and conn0.in_ring(idy0)) != (
        conn1.in_ring(idx1) and conn1.in_ring(idy1)
    ):
        return True, False

    # A direct bond was replaced by a ring path (or vice versa), leaving the
    # atoms completely disconnected in one topology (0 paths). This occurs
    # when a new ring forms and the old direct bond between two atoms is not
    # present in the ring topology.
    if (n0 == 0) != (n1 == 0):
        return True, False

    # --- Ring size changed check ---

    # Filter to paths within max_ring_size to avoid flagging macrocycle size
    # changes for rings beyond the threshold.
    short0 = [p for p in paths0 if len(p) <= max_ring_size]

    # No ring of relevant size in state 0 — nothing to compare.
    if len(short0) < 2:
        return False, False

    short1 = [p for p in paths1 if len(p) <= max_ring_size]

    # No ring of relevant size in state 1 — the ring was broken rather than
    # resized; let the ring-broken check handle this case.
    if len(short1) < 2:
        return False, False

    return False, min(len(p) for p in short0) != min(len(p) for p in short1)


def _removeDummies(molecule, is_lambda1):
    """
    Internal function which removes the dummy atoms from one of the endstates of a merged molecule.

    Parameters
    ----------

    molecule : BioSimSpace._SireWrappers.Molecule
        The molecule.

    is_lambda1 : bool
       Whether to use the molecule at lambda = 1.
    """
    from sire.legacy import IO as _SireIO
    from sire.legacy import Mol as _SireMol

    from .._Exceptions import IncompatibleError as _IncompatibleError
    from .._SireWrappers import Molecule as _Molecule

    if not molecule._is_perturbable:
        raise _IncompatibleError("'molecule' is not a perturbable molecule")

    # Always use the coordinates at lambda = 0.
    coordinates = molecule._sire_object.property("coordinates0")

    # Generate a molecule with all dummies present.
    molecule = molecule.copy()._toRegularMolecule(
        is_lambda1=is_lambda1, generate_intrascale=True
    )

    # Remove the parameters property, if it exists.
    if "parameters" in molecule._sire_object.property_keys():
        molecule._sire_object = (
            molecule._sire_object.edit().remove_property("parameters").commit()
        )

    # Set the coordinates to those at lambda = 0
    molecule._sire_object = (
        molecule._sire_object.edit().set_property("coordinates", coordinates).commit()
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
    selection.select_none()

    # Now add all of the nondummy atoms.
    for idx in nondummy_indices:
        selection.select(_SireMol.AtomIdx(idx))

    # Create a partial molecule and extract the atoms.
    partial_molecule = (
        _SireMol.PartialMolecule(molecule._sire_object, selection).extract().molecule()
    )

    # Remove the incorrect intrascale property.
    partial_molecule = (
        partial_molecule.edit().remove_property("intrascale").molecule().commit()
    )

    # Recreate a BioSimSpace molecule object.
    molecule = _Molecule(partial_molecule)

    # Parse the molecule as a GROMACS topology, which will recover the intrascale
    # matrix.
    gro_top = _SireIO.GroTop(molecule.toSystem()._sire_object)

    # Convert back to a Sire system.
    gro_sys = gro_top.to_system()

    # Add the intrascale property back into the merged molecule.
    edit_mol = molecule._sire_object.edit()
    edit_mol = edit_mol.set_property(
        "intrascale", gro_sys[_SireMol.MolIdx(0)].property("intrascale")
    )
    molecule = _Molecule(edit_mol.commit())

    return molecule
