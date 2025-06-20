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

"""Functionality for parameterising molecules."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

# A list of functions that will be exposed to the user. This list will be
# dynamically updated with the names of additional functions for the
# supported AMBER and Open Force Field models determined at import time.
__all__ = [
    "parameterise",
    "forceFields",
    "amberForceFields",
    "amberProteinForceFields",
    "openForceFields",
]

# A dictionary mapping AMBER protein force field names to their pdb2gmx
# compatibility. Note that the names specified below will be used for the
# parameterisation functions, so they should be suitably formatted.
_amber_protein_forcefields = {
    "ff03": True,
    "ff99": True,
    "ff99SB": False,
    "ff99SBildn": False,
    "ff14SB": False,
    "ff19SB": False,
}

from .. import _amber_home, _gmx_exe, _gmx_path, _isVerbose

from .._Exceptions import IncompatibleError as _IncompatibleError
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import Atom as _Atom
from .._SireWrappers import Molecule as _Molecule
from ..Solvent import waterModels as _waterModels
from ..Types import Charge as _Charge
from ..Types import Length as _Length
from .._Utils import _try_import, _have_imported
from .. import _Utils

from ._process import Process as _Process
from . import _Protocol


def parameterise(
    molecule,
    forcefield,
    ensure_compatible=True,
    work_dir=None,
    property_map={},
    **kwargs,
):
    """
    Parameterise a molecule using a specified force field.

    Parameters
    ----------

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
        The molecule to parameterise, either as a Molecule object or SMILES
        string.

    forcefield : str
        The force field. Run BioSimSpace.Parameters.forceFields() to get a
        list of the supported force fields.

    ensure_compatible : bool
        Whether to ensure that the topology of the parameterised molecule is
        compatible with that of the original molecule. An exception will be
        raised if this isn't the case, e.g. if atoms have been added. Set
        this to False is parameterising lone oxygen atoms corresponding to
        structural (crystal) water molecules. When True, the parameterised
        molecule will preserve the topology of the original molecule, e.g.
        the original atom and residue names will be kept.

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    kwargs : dict
        A dictionary of additional keyword arguments required for specific
        parameterisation functions.

    Returns
    -------

    process : :class:`Process <BioSimSpace.Parameters._process.Process>`
        A process to parameterise the molecule in the background. Call the
        .getMolecule() method on the returned process to block until the
        parameterisation is complete and get the parameterised molecule.
    """

    if not isinstance(forcefield, str):
        raise TypeError("'forcefield' must be of type 'str'")
    else:
        # Strip whitespace and convert to lower case.
        forcefield = forcefield.replace(" ", "").lower()

        if forcefield not in _forcefields_lower:
            raise ValueError("Supported force fields are: %s" % forceFields())

    return _forcefield_dict[forcefield](
        molecule, work_dir=work_dir, property_map=property_map, **kwargs
    )


def _parameterise_amber_protein(
    forcefield,
    molecule,
    tolerance=1.2,
    max_distance=_Length(6, "A"),
    water_model=None,
    pre_mol_commands=None,
    post_mol_commands=None,
    bonds=None,
    ensure_compatible=True,
    work_dir=None,
    property_map={},
    **kwargs,
):
    """
    Parameterise using the named AMBER protein force field.

    Parameters
    ----------

    forcefield : str
        The name of the AMBER force field.

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
        The molecule to parameterise, either as a Molecule object or SMILES
        string.

    tolerance : float
        The tolerance used when searching for disulphide bonds. Atoms will
        be considered to be bonded if they are a distance of less than
        tolerance times the sum of the equilibrium bond radii apart.

    max_distance : :class:`Length <BioSimSpace.Types.Length>`
        The maximum distance between atoms when searching for disulphide
        bonds.

    water_model : str
        The water model used to parameterise any structural ions.
        Run 'BioSimSpace.Solvent.waterModels()' to see the supported
        water models. This can also be used to parameterise water molecules,
        or lone oxygen atoms corresponding to structural (crystal) water
        molecules.

    pre_mol_commands : [str]
        A list of custom LEaP commands to be executed before loading the molecule.
        This can be used for loading custom parameter files, or sourcing additional
        scripts. Make sure to use absolute paths when specifying any external files.
        When this option is set, we can no longer fall back on GROMACS's pdb2gmx.

    post_mol_commands : [str]
        A list of custom LEaP commands to be executed after loading the molecule.
        This allows the use of additional commands that operate on the molecule,
        which is named "mol". When this option is set, we can no longer fall
        back on GROMACS's pdb2gmx.

    bonds : ((class:`Atom <BioSimSpace._SireWrappers.Atom>`, class:`Atom <BioSimSpace._SireWrappers.Atom>`))
        An optional tuple of atom pairs to specify additional atoms that
        should be bonded.

    ensure_compatible : bool
        Whether to ensure that the topology of the parameterised molecule is
        compatible with that of the original molecule. An exception will be
        raised if this isn't the case, e.g. if atoms have been added. Set
        this to False is parameterising lone oxygen atoms corresponding to
        structural (crystal) water molecules. When True, the parameterised
        molecule will preserve the topology of the original molecule, e.g.
        the original atom and residue names will be kept.

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    process : :class:`Process <BioSimSpace.Parameters._process.Process>`
        A process to parameterise the molecule in the background. Call the
        .getMolecule() method on the returned process to block until the
        parameterisation is complete and get the parameterised molecule.
    """

    if not isinstance(forcefield, str):
        raise TypeError("'forcefield' must be of type 'str'.")

    if forcefield not in _amber_protein_forcefields.keys():
        raise ValueError(
            f"Unsupported AMBER forcefield '{forcefield}' "
            f"options are {', '.join(_amber_protein_forcefields.keys())}."
        )

    # Extract the pdb2gmx support flag.
    is_pdb2gmx = _amber_protein_forcefields[forcefield]

    if _amber_home is None:
        if is_pdb2gmx and (_gmx_exe is None or _gmx_path is None):
            raise _MissingSoftwareError(
                "'BioSimSpace.Parameters.ff99' is not supported. "
                "Please install AmberTools (http://ambermd.org) or "
                "GROMACS (http://www.gromacs.org)."
            )
        else:
            raise _MissingSoftwareError(
                "'BioSimSpace.Parameters.ff99' is not supported. "
                "Please install AmberTools (http://ambermd.org)."
            )

    # Validate arguments.
    _validate(
        molecule=molecule,
        tolerance=tolerance,
        max_distance=max_distance,
        water_model=water_model,
        check_ions=True,
        pre_mol_commands=pre_mol_commands,
        post_mol_commands=post_mol_commands,
        bonds=bonds,
        ensure_compatible=ensure_compatible,
        property_map=property_map,
    )

    # Create the protocol.
    protocol = _Protocol.AmberProtein(
        forcefield=forcefield,
        pdb2gmx=is_pdb2gmx,
        tolerance=tolerance,
        water_model=water_model,
        max_distance=max_distance,
        pre_mol_commands=pre_mol_commands,
        post_mol_commands=post_mol_commands,
        bonds=bonds,
        ensure_compatible=ensure_compatible,
        property_map=property_map,
    )

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)


def gaff(
    molecule,
    work_dir=None,
    net_charge=None,
    charge_method="BCC",
    ensure_compatible=True,
    property_map={},
    **kwargs,
):
    """
    Parameterise using the GAFF force field.

    Parameters
    ----------

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
        The molecule to parameterise, either as a Molecule object or SMILES
        string.

    net_charge : int, :class:`Charge <BioSimSpace.Types.Charge>`
        The net charge on the molecule.

    charge_method : str
            The method to use when calculating atomic charges:
            "RESP", "CM2", "MUL", "BCC", "ESP", "GAS"

    ensure_compatible : bool
        Whether to ensure that the topology of the parameterised molecule is
        compatible with that of the original molecule. An exception will be
        raised if this isn't the case, e.g. if atoms have been added. When True,
        the parameterised molecule will preserve the topology of the original
        molecule, e.g. the original atom and residue names will be kept.

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    process : :class:`Process <BioSimSpace.Parameters._process.Process>`
        A process to parameterise the molecule in the background. Call the
        .getMolecule() method on the returned process to block until the
        parameterisation is complete and get the parameterised molecule.
    """

    if _amber_home is None:
        raise _MissingSoftwareError(
            "'BioSimSpace.Parameters.gaff' is not supported. "
            "Please install AmberTools (http://ambermd.org)."
        )

    # Validate arguments.
    _validate(
        molecule=molecule,
        check_ions=False,
        ensure_compatible=ensure_compatible,
        property_map=property_map,
    )

    if net_charge is not None:
        # Get the value of the charge.
        if isinstance(net_charge, _Charge):
            net_charge = net_charge.value()

        if isinstance(net_charge, float):
            if net_charge % 1 != 0:
                raise ValueError("'net_charge' must be integer valued.")

        # Try to convert to int.
        try:
            net_charge = int(net_charge)
        except:
            raise TypeError(
                "'net_charge' must be of type 'int', or `BioSimSpace.Types.Charge'"
            )

    # Create a default protocol.
    protocol = _Protocol.GAFF(
        version=1,
        net_charge=net_charge,
        charge_method=charge_method,
        ensure_compatible=ensure_compatible,
        property_map=property_map,
        **kwargs,
    )

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)


def gaff2(
    molecule,
    work_dir=None,
    net_charge=None,
    charge_method="BCC",
    ensure_compatible=True,
    property_map={},
    **kwargs,
):
    """
    Parameterise using the GAFF2 force field.

    Parameters
    ----------

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
        The molecule to parameterise, either as a Molecule object or SMILES
        string.

    net_charge : int, :class:`Charge <BioSimSpace.Types.Charge>`
        The net charge on the molecule.

    charge_method : str
            The method to use when calculating atomic charges:
            "RESP", "CM2", "MUL", "BCC", "ESP", "GAS"

    ensure_compatible : bool
        Whether to ensure that the topology of the parameterised molecule is
        compatible with that of the original molecule. An exception will be
        raised if this isn't the case, e.g. if atoms have been added. When True,
        the parameterised molecule will preserve the topology of the original
        molecule, e.g. the original atom and residue names will be kept.

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    process : :class:`Process <BioSimSpace.Parameters._process.Process>`
        A process to parameterise the molecule in the background. Call the
        .getMolecule() method on the returned process to block until the
        parameterisation is complete and get the parameterised molecule.
    """

    if _amber_home is None:
        raise _MissingSoftwareError(
            "'BioSimSpace.Parameters.gaff2' is not supported. "
            "Please install AmberTools (http://ambermd.org)."
        )

    # Validate arguments.
    _validate(
        molecule=molecule,
        check_ions=False,
        ensure_compatible=ensure_compatible,
        property_map=property_map,
    )

    if net_charge is not None:
        # Get the value of the charge.
        if isinstance(net_charge, _Charge):
            net_charge = net_charge.value()

        if isinstance(net_charge, float):
            if net_charge % 1 != 0:
                raise ValueError("'net_charge' must be integer valued.")

        # Try to convert to int.
        try:
            net_charge = int(net_charge)
        except:
            raise TypeError(
                "'net_charge' must be of type 'int', or `BioSimSpace.Types.Charge'"
            )

        if net_charge % 1 != 0:
            raise ValueError("'net_charge' must be integer valued.")

    # Create a default protocol.
    protocol = _Protocol.GAFF(
        version=2,
        net_charge=net_charge,
        charge_method=charge_method,
        ensure_compatible=ensure_compatible,
        property_map=property_map,
        **kwargs,
    )

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)


def _parameterise_openff(
    forcefield,
    molecule,
    ensure_compatible=True,
    use_nagl=True,
    work_dir=None,
    property_map={},
    **kwargs,
):
    """
    Parameterise a molecule using a force field from the Open Force Field
    Initiative.

    Parameters
    ----------

    forcefield : str
        The force field. Run BioSimSpace.Parameters.openForceFields() to get a
        list of the supported force fields.

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
        The molecule to parameterise, either as a Molecule object or SMILES
        string.

    ensure_compatible : bool
        Whether to ensure that the topology of the parameterised molecule is
        compatible with that of the original molecule. An exception will be
        raised if this isn't the case, e.g. if atoms have been added. When True,
        the parameterised molecule will preserve the topology of the original
        molecule, e.g. the original atom and residue names will be kept.

    use_nagl : bool
        Whether to use NAGL to compute AM1-BCC charges. If False, the default
        is to use AmberTools via antechamber and sqm. (This option is only
        used if NAGL is available.)

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    process : :class:`Process <BioSimSpace.Parameters._process.Process>`
        A process to parameterise the molecule in the background. Call the
        .getMolecule() method on the returned process to block until the
        parameterisation is complete and get the parameterised molecule.
    """

    from sire.legacy.Base import findExe as _findExe

    try:
        _findExe("antechamber")
    except:
        raise _MissingSoftwareError(
            f"'{forcefield}' is not supported. AmberTools "
            "(http://ambermd.org) is needed for charge "
            "calculation and 'antechamber' executable "
            "must be in your PATH."
        ) from None

    # Validate arguments.

    if not isinstance(molecule, (_Molecule, str)):
        raise TypeError(
            "'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule' or 'str'"
        )

    if not isinstance(forcefield, str):
        raise TypeError("'forcefield' must be of type 'str'")
    else:
        # Strip whitespace and convert to lower case.
        forcefield = forcefield.replace(" ", "").lower()

        if forcefield not in _forcefields_lower:
            raise ValueError("Supported force fields are: %s" % openForceFields())

    if not isinstance(ensure_compatible, bool):
        raise TypeError("'ensure_compatible' must be of type 'bool'.")

    if not isinstance(use_nagl, bool):
        raise TypeError("'use_nagl' must be of type 'bool'.")

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    # Create a default protocol.
    protocol = _Protocol.OpenForceField(
        forcefield,
        ensure_compatible=ensure_compatible,
        use_nagl=use_nagl,
        property_map=property_map,
    )

    # Run the parameterisation protocol in the background and return
    # a handle to the thread.
    return _Process(molecule, protocol, work_dir=work_dir, auto_start=True)


def forceFields():
    """
    Return a list of the supported force fields.

    Returns
    -------

    force_fields : [str]
        A list of the supported force fields.
    """
    return _forcefields


def amberForceFields():
    """
    Return a list of the supported AMBER force fields.

    Returns
    -------

    force_fields : [str]
        A list of the supported AMBER force fields.
    """
    return list(_amber_protein_forcefields.keys()) + ["gaff", "gaff2"]


def amberProteinForceFields():
    """
    Return a list of the supported AMBER protein force fields.

    Returns
    -------

    force_fields : [str]
        A list of the supported AMBER protein force fields.
    """
    return list(_amber_protein_forcefields.keys())


def openForceFields():
    """
    Return a list of the supported force fields from the Open Force Field
    Initiative.

    Returns
    -------

    force_fields : [str]
        A list of the supported force fields from the Open Force Field Initiative.
    """
    return _open_forcefields


def _validate_water_model(water_model):
    """
    Internal helper function used to validate the water model chosen
    for parameterising structural ions.

    Parameters
    ----------

    water_model : str
        The chosen water model.

    Returns
    -------

    is_valid : bool
        Whether the water model is supported.
    """
    # Srip whitespace and convert to lower case.
    water_model = water_model.replace(" ", "").lower()

    # Check that we support this water model.
    if water_model in _waterModels():
        return True
    else:
        return False


def _has_ions(molecule, property_map={}):
    """
    Internal helper function to check whether a molecule contains ions.

    Parameters
    ----------

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
        A molecule object.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }

    Returns
    -------

    has_ions : bool
        Whether the molecule contains ions.

    ions : [str]
        A list of the ions that were found.
    """

    # Store a list of ionic elements.
    # (Taken from an AMBER leaprc.water.* file.)

    elements = [
        "F",
        "Cl",
        "Br",
        "I",
        "Li",
        "Na",
        "K",
        "Rb",
        "Cs",
        "Mg",
        "Tl",
        "Cu",
        "Ag",
        "Be",
        "Cu",
        "Ni",
        "Pt",
        "Zn",
        "Co",
        "Pd",
        "Ag",
        "Cr",
        "Fe",
        "Mg",
        "V",
        "Mn",
        "Hg",
        "Cd",
        "Yb",
        "Ca",
        "Sn",
        "Pb",
        "Eu",
        "Sr",
        "Sm",
        "Ba",
        "Ra",
        "Al",
        "Fe",
        "Cr",
        "In",
        "Tl",
        "Y",
        "La",
        "Ce",
        "Pr",
        "Nd",
        "Sm",
        "Eu",
        "Gd",
        "Tb",
        "Dy",
        "Er",
        "Tm",
        "Lu",
        "Hf",
        "Zr",
        "Ce",
        "U",
        "Pu",
        "Th",
    ]

    # We need to search for ions individually since Sire can't
    # handle or'ed search strings beyond a certain size.

    # A list of ions that we've found.
    ions = []

    # Whether the molecule is a string.
    if isinstance(molecule, str):
        is_string = True
        molecule = molecule.upper()
    else:
        is_string = False

    for element in elements:
        if is_string:
            if element.upper() in molecule:
                ions.append(element)
        else:
            try:
                molecule.search(f"element {element}", property_map=property_map)
                ions.append(element)
            except:
                pass

    # Check whether we found any ions.
    if len(ions) > 0:
        return True, ions
    else:
        return False, ions


def _validate(
    molecule=None,
    tolerance=None,
    max_distance=None,
    water_model=None,
    check_ions=False,
    pre_mol_commands=None,
    post_mol_commands=None,
    bonds=None,
    ensure_compatible=True,
    work_dir=None,
    property_map=None,
):
    """
    Internal function to validate arguments.

    Parameters
    ----------

    molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
        The molecule to parameterise, either as a Molecule object or SMILES
        string.

    tolerance : float
        The tolerance used when searching for disulphide bonds. Atoms will
        be considered to be bonded if they are a distance of less than
        tolerance times the sum of the equilibrium bond radii apart.

    max_distance : :class:`Length <BioSimSpace.Types.Length>`
        The maximum distance between atoms when searching for disulphide
        bonds.

    water_model : str
        The water model used to parameterise any structural ions.
        Run 'BioSimSpace.Solvent.waterModels()' to see the supported
        water models. This can also be used to parameterise water molecules,
        or lone oxygen atoms corresponding to structural (crystal) water
        molecules.

    check_ions : bool
        Whether to check for the presence of structural ions. This is only
        required when parameterising with protein force fields.

    pre_mol_commands : [str]
        A list of custom LEaP commands to be executed before loading the molecule.
        This can be used for loading custom parameter files, or sourcing additional
        scripts. Make sure to use absolute paths when specifying any external files.
        When this option is set, we can no longer fall back on GROMACS's pdb2gmx.

    post_mol_commands : [str]
        A list of custom LEaP commands to be executed after loading the molecule.
        This allows the use of additional commands that operate on the molecule,
        which is named "mol". When this option is set, we can no longer fall
        back on GROMACS's pdb2gmx.

    bonds : ((class:`Atom <BioSimSpace._SireWrappers.Atom>`, class:`Atom <BioSimSpace._SireWrappers.Atom>`))
        An optional tuple of atom pairs to specify additional atoms that
        should be bonded.

    ensure_compatible : bool
        Whether to ensure that the topology of the parameterised molecule is
        compatible with that of the original molecule. An exception will be
        raised if this isn't the case, e.g. if atoms have been added. Set this
        to False is parameterising lone oxygen atoms corresponding to
        structural (crystal) water molecules. When True, the parameterised
        molecule will preserve the topology of the original molecule, e.g.
        the original atom and residue names will be kept.

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. { "charge" : "my-charge" }
    """

    if molecule is not None:
        if not isinstance(molecule, (_Molecule, str)):
            raise TypeError(
                "'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule' or 'str'"
            )

    if tolerance is not None:
        if not isinstance(tolerance, (int, float)):
            raise TypeError("'tolerance' must be of type 'float'.")
        if tolerance < 1:
            raise ValueError("'tolerance' must be >= 1.0.")

    if max_distance is not None:
        if not isinstance(max_distance, _Length):
            raise TypeError(
                "'max_distance' must be of type 'BioSimSpace.Types.Length'."
            )

    if water_model is not None:
        if not isinstance(water_model, str):
            raise TypeError("'water_model' must be of type 'str'.")
        if not _validate_water_model(water_model):
            water_models = ", ".join(_waterModels())
            raise ValueError(
                f"'{water_model}' is unsupported. Supported models are: {water_models}"
            )
    elif check_ions:
        has_ions, ions = _has_ions(molecule, property_map=property_map)
        if has_ions:
            ion_string = ", ".join(ions)
            raise ValueError(
                f"The molecule contains the following ions: {ion_string}. "
                "Please choose a 'water_model' for the ion parameters."
            )

    if pre_mol_commands is not None:
        if not isinstance(pre_mol_commands, (list, tuple)):
            raise TypeError("'pre_mol_commands' must be a 'list' of 'str' types.")
        else:
            if not all(isinstance(x, str) for x in pre_mol_commands):
                raise TypeError("'pre_mol_commands' must be a 'list' of 'str' types.")

    if post_mol_commands is not None:
        if not isinstance(post_mol_commands, (list, tuple)):
            raise TypeError("'post_mol_commands' must be a 'list' of 'str' types.")
        else:
            if not all(isinstance(x, str) for x in post_mol_commands):
                raise TypeError("'post_mol_commands' must be a 'list' of 'str' types.")

    if bonds is not None:
        if molecule is None:
            raise ValueError("Cannot add bonds when no 'molecule' is specified!")
        else:
            if not isinstance(bonds, (tuple, list)):
                raise TypeError("'bonds' must be of type 'tuple' or 'list'.")
            for bond in bonds:
                if not isinstance(bond, (tuple, list)):
                    raise TypeError(
                        "Each bond entry must be a 'tuple' or 'list' of atom pairs."
                    )
                else:
                    if len(bond) != 2:
                        raise ValueError("Each 'bonds' entry must contain two items.")
                    else:
                        # Extract the atoms in the bond.
                        atom0, atom1 = bond

                        # Make sure these are atoms.
                        if not isinstance(atom0, _Atom) or not isinstance(atom1, _Atom):
                            raise TypeError(
                                "'bonds' must contain tuples of "
                                "'BioSimSpace._SireWrappers.Atom' types."
                            )

                        # Make sure that they belong to the molecule being parameterised.
                        if not (
                            atom0._sire_object.molecule() == molecule._sire_object
                        ) or not (
                            atom1._sire_object.molecule() == molecule._sire_object
                        ):
                            raise ValueError(
                                "Atoms in 'bonds' don't belong to the 'molecule'."
                            )

    if not isinstance(ensure_compatible, bool):
        raise TypeError("'ensure_compatible' must be of type 'bool'.")

    if property_map is not None:
        if not isinstance(property_map, dict):
            raise TypeError("'property_map' must be of type 'dict'")


# Wrapper function to dynamically generate functions with a given name.
# Here "name" refers to the name of a supported AMBER protein force field.
# We could do this with functools.partial and functions.update_wrapper, but
# this approach preserves correct documentation for the partial function,
# both within the Python console, or via Sphinx.
import sys as _sys

_namespace = _sys.modules[__name__]


def _make_amber_protein_function(name):
    def _function(
        molecule,
        tolerance=1.2,
        max_distance=_Length(6, "A"),
        water_model=None,
        pre_mol_commands=None,
        post_mol_commands=None,
        bonds=None,
        ensure_compatible=True,
        work_dir=None,
        property_map={},
    ):
        """
        Parameterise a molecule using the named AMBER force field.

        Parameters
        ----------

        molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
            The molecule to parameterise, either as a Molecule object or SMILES
            string.

        tolerance : float
            The tolerance used when searching for disulphide bonds. Atoms will
            be considered to be bonded if they are a distance of less than
            tolerance times the sum of the equilibrium bond radii apart.

        max_distance : :class:`Length <BioSimSpace.Types.Length>`
            The maximum distance between atoms when searching for disulphide
            bonds.

        water_model : str
            The water model used to parameterise any structural ions.
            Run 'BioSimSpace.Solvent.waterModels()' to see the supported
            water models. This can also be used to parameterise water molecules,
            or lone oxygen atoms corresponding to structural (crystal) water
            molecules.

        pre_mol_commands : [str]
            A list of custom LEaP commands to be executed before loading the molecule.
            This can be used for loading custom parameter files, or sourcing additional
            scripts. Make sure to use absolute paths when specifying any external files.
            When this option is set, we can no longer fall back on GROMACS's pdb2gmx.

        post_mol_commands : [str]
            A list of custom LEaP commands to be executed after loading the molecule.
            This allows the use of additional commands that operate on the molecule,
            which is named "mol". When this option is set, we can no longer fall
            back on GROMACS's pdb2gmx.

        bonds : ((class:`Atom <BioSimSpace._SireWrappers.Atom>`, class:`Atom <BioSimSpace._SireWrappers.Atom>`))
            An optional tuple of atom pairs to specify additional atoms that
            should be bonded.

        ensure_compatible : bool
            Whether to ensure that the topology of the parameterised molecule is
            compatible with that of the original molecule. An exception will be
            raised if this isn't the case, e.g. if atoms have been added. Set
            this to False is parameterising lone oxygen atoms corresponding to
            structural (crystal) water molecules. When True, the parameterised
            molecule will preserve the topology of the original molecule, e.g.
            the original atom and residue names will be kept.

        work_dir : str
             The working directory for the process.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        process : :class:`Process <BioSimSpace.Parameters._process.Process>`
            A process to parameterise the molecule in the background. Call the
            .getMolecule() method on the returned process to block until the
            parameterisation is complete and get the parameterised molecule.
        """
        return _parameterise_amber_protein(
            name,
            molecule,
            tolerance=tolerance,
            max_distance=max_distance,
            water_model=water_model,
            pre_mol_commands=pre_mol_commands,
            post_mol_commands=post_mol_commands,
            bonds=bonds,
            ensure_compatible=ensure_compatible,
            work_dir=work_dir,
            property_map=property_map,
        )

    return _function


# Wrapper function to dynamically generate functions with a given name.
# Here "name" refers to the name of a supported force field from the Open
# Force Field Initiative. The force field name has been "tidied" so that
# it conforms to sensible function naming standards, i.e. "-" and "."
# characters replaced by underscores.
def _make_openff_function(name):
    def _function(
        molecule, ensure_compatible=True, use_nagl=True, work_dir=None, property_map={}
    ):
        """
        Parameterise a molecule using the named force field from the
        Open Force Field initiative.

        Parameters
        ----------

        molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
            The molecule to parameterise, either as a Molecule object or SMILES
            string.

        ensure_compatible : bool
            Whether to ensure that the topology of the parameterised molecule is
            compatible with that of the original molecule. An exception will be
            raised if this isn't the case, e.g. if atoms have been added. Set
            this to False is parameterising lone oxygen atoms corresponding to
            structural (crystal) water molecules. When True, the parameterised
            molecule will preserve the topology of the original molecule, e.g.
            the original atom and residue names will be kept.

        use_nagl : bool
            Whether to use NAGL to compute AM1-BCC charges. If False, the default
            is to use AmberTools via antechamber and sqm. (This option is only
            used if NAGL is available.)

        work_dir : str
            The working directory for the process.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        Returns
        -------

        process : :class:`Process <BioSimSpace.Parameters._process.Process>`
            A process to parameterise the molecule in the background. Call the
            .getMolecule() method on the returned process to block until the
            parameterisation is complete and get the parameterised molecule.
        """
        return _parameterise_openff(
            name,
            molecule,
            ensure_compatible=ensure_compatible,
            use_nagl=use_nagl,
            work_dir=work_dir,
            property_map=property_map,
        )

    return _function


# Dynamically create functions for all of the supported AMBER protein force fields.
for _ff in _amber_protein_forcefields.keys():
    # Generate the function and bind it to the namespace.
    _function = _make_amber_protein_function(_ff)
    _function.__name__ = _ff
    _function.__qualname = _ff
    setattr(_namespace, _ff, _function)

    # Expose the function to the user.
    __all__.append(_ff)

# Expose the two GAFF functions.
__all__.append("gaff")
__all__.append("gaff2")

# Create a list of the force field names (to date.)
# This needs to come after all of the force field functions.
_forcefields = []  # List of force fields (actual names).
_forcefields_lower = []  # List of lower case names.
_forcefield_dict = {}  # Mapping between lower case names and functions.
for _ff in list(_amber_protein_forcefields.keys()) + ["gaff", "gaff2"]:
    _forcefields.append(_ff)
    _forcefields_lower.append(_ff.lower())
    _forcefield_dict[_ff.lower()] = getattr(_namespace, _ff)

# Dynamically create functions for all available force fields from the Open
# Force Field Initiative.
_openforcefields = _try_import("openforcefields")

if _have_imported(_openforcefields):
    from glob import glob as _glob
    import os as _os

    _openff_dirs = _openforcefields.get_forcefield_dirs_paths()
    _open_forcefields = []
    # Loop over all force field directories.
    for _dir in _openff_dirs:
        # Glob all offxml files in the directory.
        _ff_list = _glob(f"{_dir}" + "/*.offxml")
        for _ff in _ff_list:
            # Get the force field name (base name minus extension).
            _base = _os.path.basename(_ff)
            _ff = _os.path.splitext(_base)[0]

            # Only include unconstrained force-fields since we need to go via
            # an intermediate ParmEd conversion. This means ParmEd must receive
            # bond parameters. We can then choose to constrain later, if required.
            # See, e.g: https://github.com/openforcefield/openff-toolkit/issues/603

            if "unconstrained" in _ff:
                # Append to the list of available force fields.
                _forcefields.append(_ff)
                _open_forcefields.append(_ff)

                # Create a sane function name, i.e. replace "-" and "."
                # characters with "_".
                _func_name = _ff.replace("-", "_")
                _func_name = _func_name.replace(".", "_")

                # Generate the function and bind it to the namespace.
                _function = _make_openff_function(_ff)
                _function.__name__ = _func_name
                _function.__qualname = _func_name
                setattr(_namespace, _func_name, _function)

                # Expose the function to the user.
                __all__.append(_func_name)

                # Convert force field name to lower case and map to its function.
                _forcefields_lower.append(_ff.lower())
                _forcefield_dict[_ff.lower()] = getattr(_namespace, _func_name)

    # Clean up redundant attributes.
    del _base
    del _dir
    del _ff_list
    del _func_name
    del _glob
    del _make_openff_function
    del _openff_dirs
    del _os
elif _isVerbose():
    print("openforcefields not available as this module cannot be loaded.")

# Clean up redundant attributes.
del _ff
del _function
del _namespace
del _openforcefields
del _sys
