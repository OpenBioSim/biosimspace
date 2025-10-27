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

"""
Functionality for handling parameterisation protocols
for AMBER force field models.
Author: Lester Hedges <lester.hedges@gmail.com>.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["AmberProtein", "GAFF"]


import os as _os
import subprocess as _subprocess

from ..._Utils import _try_import, _have_imported

# Temporarily redirect stderr to suppress import warnings.
import sys as _sys

_orig_stderr = _sys.stderr
_sys.stderr = open(_os.devnull, "w")

_openff = _try_import("openff")

if _have_imported(_openff):
    from openff.toolkit.topology import Molecule as _OpenFFMolecule
else:
    _OpenFFMolecule = _openff

# Reset stderr.
_sys.stderr = _orig_stderr
del _sys, _orig_stderr


from ... import _amber_home
from ... import _Utils
from ...Types import Length as _Length

from . import _protocol

# Set the tLEaP cmd directory.
if _amber_home is not None:
    _cmd_dir = _amber_home + "/dat/leap/cmd"
    if not _os.path.isdir(_cmd_dir):
        raise IOError("Missing tLEaP dat directory: '%s'" % _cmd_dir)

# Search for all of the required executables.

# Search for the tLEaP exe.

_tleap_exe = None

if _amber_home is not None:
    _tleap_exe = "%s/bin/tleap" % _amber_home
    if not _os.path.isfile(_tleap_exe):
        raise IOError("Missing tLEaP executable: '%s'" % _tleap_exe)

# Search for the Antechamber exe.

_antechamber_exe = None

if _amber_home is not None:
    _antechamber_exe = "%s/bin/antechamber" % _amber_home
    if not _os.path.isfile(_antechamber_exe):
        raise IOError("Missing Antechamber executable: '%s'" % _antechamber_exe)

# Search for the parmchk exe.

_parmchk_exe = None

if _amber_home is not None:
    _parmchk_exe = "%s/bin/parmchk2" % _amber_home
    if not _os.path.isfile(_parmchk_exe):
        raise IOError("Missing parmchk executable: '%s'" % _parmchk_exe)

# Check for ABCG2 support.
_has_abcg2 = False

try:
    cmd = "antechamber -L"
    proc = _subprocess.run(
        _Utils.command_split(cmd),
        shell=False,
        stdout=_subprocess.PIPE,
        stderr=_subprocess.PIPE,
        text=True,
    )
    if proc.returncode == 0:
        # Check the output for ABCG2 support.
        if "abcg2" in proc.stdout:
            _has_abcg2 = True
except:
    pass


class AmberProtein(_protocol.Protocol):
    """A class for handling AMBER protein force field models."""

    def __init__(
        self,
        forcefield,
        pdb2gmx=False,
        tolerance=1.2,
        max_distance=_Length(6, "A"),
        water_model=None,
        pre_mol_commands=None,
        post_mol_commands=None,
        bonds=None,
        ensure_compatible=True,
        property_map={},
    ):
        """
        Constructor.

        Parameters
        ----------

        forcefield : str
            The name of the force field.

        pdb2gmx : bool
            Whether the force field is supported by pdb2gmx.

        tolerance : float
            The tolerance used when searching for disulphide bonds. Atoms will
            be considered to be bonded if they are a distance of less than
            tolerance times the sum of the equilibrium bond radii apart.

        max_distance : :class:`Length <BioSimSpace.Types.Length>`
            The maximum distance between atoms when searching for disulphide
            bonds.

        water_model : str
            The water model used to parameterise any structural ions. This
            will be ignored when it is not supported by the chosen force field.
            Run 'BioSimSpace.Solvent.waterModels()' to see the supported
            water models.

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
            should be bonded. This is useful when the PDB CONECT record is
            incomplete.

        ensure_compatible : bool
            Whether to ensure that the topology of the parameterised molecule is
            compatible with that of the original molecule. An exception will be
            raised if this isn't the case, e.g. if atoms have been added. Set
            this to False is parameterising lone oxygen atoms corresponding to
            structural (crystal) water molecules. When True, the parameterised
            molecule will preserve the topology of the original molecule, e.g.
            the original atom and residue names will be kept.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """
        from ..._SireWrappers import Atom as _Atom

        # Call the base class constructor.
        super().__init__(
            forcefield=forcefield,
            ensure_compatible=ensure_compatible,
            property_map=property_map,
        )

        # Validate the pdb2gmx compatibility flag.
        if not isinstance(pdb2gmx, bool):
            raise TypeError("'pdb2gmx' must be of type 'bool'.")
        else:
            self._pdb2gmx = pdb2gmx

        # Validate the tolerance.
        if not isinstance(tolerance, (int, float)):
            raise TypeError("'tolerance' must be of type 'float'")
        tolerance = float(tolerance)
        if tolerance < 1:
            raise ValueError("'tolerance' must be >= 1.0.")
        self._tolerance = tolerance

        # Validate the max distance.
        if not isinstance(max_distance, _Length):
            raise ValueError(
                "'max_distance' must be of type 'BioSimSpace.Types.Length'"
            )
        self._max_distance = max_distance

        # Validate the water model.
        if water_model is not None and not isinstance(water_model, str):
            raise TypeError("'water_model' must be of type 'str'")
        else:
            self._water_model = water_model

        # Validate the additional leap commands.
        if pre_mol_commands is not None:
            if not isinstance(pre_mol_commands, (list, tuple)):
                raise TypeError("'pre_mol_commands must be a 'list' of 'str' types.")
            else:
                if not all(isinstance(x, str) for x in pre_mol_commands):
                    raise TypeError(
                        "'pre_mol_commands' must be a 'list' of 'str' types."
                    )
        self._pre_mol_commands = pre_mol_commands
        if post_mol_commands is not None:
            if not isinstance(post_mol_commands, (list, tuple)):
                raise TypeError("'post_mol_commands must be a 'list' of 'str' types.")
            else:
                if not all(isinstance(x, str) for x in post_mol_commands):
                    raise TypeError(
                        "'post_mol_commands' must be a 'list' of 'str' types."
                    )
        self._post_mol_commands = post_mol_commands

        # Validate the bond records.
        if bonds is not None:
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
        self._bonds = bonds

        # Set the compatibility flags.
        self._tleap = True
        if self._pre_mol_commands is not None or self._post_mol_commands is not None:
            self._pdb2gmx = False

    def run(self, molecule, work_dir=None, queue=None):
        """
        Run the parameterisation protocol.

        Parameters
        ----------

        molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
            The molecule to parameterise, either as a Molecule object or SMILES
            string.

        work_dir : :class:`WorkDir <BioSimSpace._Utils.WorkDir>`
            The working directory.

        queue : queue.Queue
            The thread queue is which this method has been run.

        Returns
        -------

        molecule : BioSimSpace._SireWrappers.Molecule
            The parameterised molecule.
        """
        from ... import _gmx_exe, _isVerbose
        from ...Convert import smiles as _smiles
        from ..._Exceptions import MissingSoftwareError as _MissingSoftwareError
        from ... import IO as _IO
        from ..._SireWrappers import Molecule as _Molecule
        from ..._Exceptions import ThirdPartyError as _ThirdPartyError
        import queue as _queue

        if not isinstance(molecule, (_Molecule, str)):
            raise TypeError(
                "'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule' or 'str'"
            )

        if work_dir is not None and not isinstance(work_dir, _Utils.WorkDir):
            raise TypeError("'work_dir' must be of type 'BioSimSpace._Utils.WorkDir'")

        if queue is not None and not isinstance(queue, _queue.Queue):
            raise TypeError("'queue' must be of type 'queue.Queue'")

        # Set work_dir to the current directory.
        if work_dir is None:
            work_dir = _os.getcwd()

        # Try to create a molecule from the SMILES string.
        if isinstance(molecule, str):
            try:
                smiles = molecule
                molecule = _smiles(molecule)
                edit_mol = molecule._sire_object.edit()
                edit_mol = edit_mol.rename(f"smiles:{smiles}").molecule()
                molecule._sire_object = edit_mol.commit()
            except Exception as e:
                msg = "Unable to convert SMILES to Molecule using RDKit."
                if _isVerbose():
                    msg += ": " + getattr(e, "message", repr(e))
                    raise _ThirdPartyError(msg) from e
                else:
                    raise _ThirdPartyError(msg) from None

        # Create a copy of the molecule.
        new_mol = molecule.copy()

        # Choose the program to run with depending on the force field compatibility.
        # If tLEaP and pdb2gmx are supported, default to tLEaP, then use pdb2gmx if
        # tLEaP fails to produce output.

        # First, try parameterise using tLEaP.
        if self._tleap:
            if _tleap_exe is not None:
                output = self._run_tleap(molecule, str(work_dir))
                new_mol._ion_water_model = self._water_model
            # Otherwise, try using pdb2gmx.
            elif self._pdb2gmx:
                if _gmx_exe is not None:
                    output = self._run_pdb2gmx(molecule, str(work_dir))
                else:
                    raise _MissingSoftwareError(
                        "Cannot parameterise. Missing AmberTools and GROMACS."
                    )

        # Parameterise using pdb2gmx.
        elif self._pdb2gmx:
            if _gmx_exe is not None:
                output = self._run_pdb2gmx(molecule, str(work_dir))
            else:
                raise _MissingSoftwareError(
                    "Cannot use pdb2gmx since GROMACS is not installed!"
                )

        # Prepend the working directory to the output file names.
        output = [
            _os.path.join(str(work_dir), output[0]),
            _os.path.join(str(work_dir), output[1]),
        ]

        try:
            # Load the parameterised molecule. (This could be a system of molecules.)
            par_mol = _IO.readMolecules(output)
            # Extract single molecules.
            if par_mol.nMolecules() == 1:
                par_mol = par_mol.getMolecules()[0]
        except Exception as e:
            msg = "Failed to read molecule from: '%s', '%s'" % (output[0], output[1])
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Make the molecule 'mol' compatible with 'par_mol'. This will create
        # a mapping between atom indices in the two molecules and add all of
        # the new properties from 'par_mol' to 'mol'.
        if self._ensure_compatible:
            new_mol.makeCompatibleWith(
                par_mol,
                property_map=self._property_map,
                overwrite=True,
                verbose=False,
            )
        else:
            try:
                new_mol.makeCompatibleWith(
                    par_mol,
                    property_map=self._property_map,
                    overwrite=True,
                    verbose=False,
                )
            except:
                new_mol = par_mol

        # Record the forcefield used to parameterise the molecule.
        new_mol._forcefield = self._forcefield

        if queue is not None:
            queue.put(new_mol)
        return new_mol

    def _run_tleap(self, molecule, work_dir):
        """
        Run using tLEaP.

        Parameters
        ----------

        molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            The molecule to parameterise.

        work_dir : str
            The working directory.
        """
        from ..._Exceptions import ParameterisationError as _ParameterisationError
        from ... import IO as _IO
        from ..._SireWrappers import Molecule as _Molecule
        from sire.legacy import IO as _SireIO
        from ... import _isVerbose

        # Write the system to a PDB file.
        try:
            # LEaP expects residue numbering to be ascending and continuous.
            renumbered_molecule = _SireIO.renumberConstituents(
                molecule.toSystem()._sire_object
            )[0]
            renumbered_molecule = _Molecule(renumbered_molecule)
            _IO.saveMolecules(
                _os.path.join(str(work_dir), "leap"),
                renumbered_molecule,
                "pdb",
                property_map=self._property_map,
            )
        except Exception as e:
            raise
            msg = "Failed to write system to 'PDB' format."
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Try to find a force field file.
        ff = _find_force_field(self._forcefield)

        # Check to see if any disulphide bonds are present.
        disulphide_bonds = self._get_disulphide_bonds(
            renumbered_molecule._sire_object,
            self._tolerance,
            self._max_distance,
            self._property_map,
        )

        # Get any additional bond records.
        bond_records = self._generate_bond_records(renumbered_molecule, self._bonds)

        # Remove any duplicate bonds, e.g. if disulphides are specified twice.
        pruned_bond_records = []
        for bond in bond_records:
            is_duplicate = False
            a0, a1 = bond.split()[1:]
            for disulphide in disulphide_bonds:
                d0, d1 = disulphide.split()[1:]

                # Make sure the bonded atoms don't match.
                if (a0 == d0 and a1 == d1) or (a0 == d1 and a1 == d0):
                    is_duplicate = True
                    break

            # Add the bond if it is not a duplicate.
            if not is_duplicate:
                pruned_bond_records.append(bond)

        # Write the LEaP input file.
        with open(_os.path.join(str(work_dir), "leap.txt"), "w") as file:
            file.write("source %s\n" % ff)
            if self._water_model is not None:
                if self._water_model in ["tip4p", "tip5p"]:
                    file.write("source leaprc.water.tip4pew\n")
                else:
                    file.write("source leaprc.water.%s\n" % self._water_model)
            # Write pre-mol user leap commands.
            if self._pre_mol_commands is not None:
                for command in self._pre_mol_commands:
                    file.write("%s\n" % command)
            file.write("mol = loadPdb leap.pdb\n")
            # Add any disulphide bond records.
            for bond in disulphide_bonds:
                file.write("%s\n" % bond)
            # Add any additional bond records.
            for bond in pruned_bond_records:
                file.write("%s\n" % bond)
            # Write post-mol user leap commands.
            if self._post_mol_commands is not None:
                for command in self._post_mol_commands:
                    file.write("%s\n" % command)
            file.write("saveAmberParm mol leap.top leap.crd\n")
            file.write("quit")

        # Generate the tLEaP command.
        command = "%s -f leap.txt" % _tleap_exe

        with open(_os.path.join(str(work_dir), "README.txt"), "w") as file:
            # Write the command to file.
            file.write("# tLEaP was run with the following command:\n")
            file.write("%s\n" % command)

        # Create files for stdout/stderr.
        stdout = open(_os.path.join(str(work_dir), "leap.out"), "w")
        stderr = open(_os.path.join(str(work_dir), "leap.err"), "w")

        # Run tLEaP as a subprocess.
        proc = _subprocess.run(
            _Utils.command_split(command),
            cwd=work_dir,
            shell=False,
            stdout=stdout,
            stderr=stderr,
        )
        stdout.close()
        stderr.close()

        # tLEaP doesn't return sensible error codes, so we need to check that
        # the expected output was generated.
        if _os.path.isfile(
            _os.path.join(str(work_dir), "leap.top")
        ) and _os.path.isfile(_os.path.join(str(work_dir), "leap.crd")):
            # Check the output of tLEaP for missing atoms.
            if self._ensure_compatible:
                if _has_missing_atoms(_os.path.join(str(work_dir), "leap.out")):
                    raise _ParameterisationError(
                        "tLEaP added missing atoms. The topology is now "
                        "inconsistent with the original molecule. Please "
                        "make sure that your initial molecule has a "
                        "complete topology."
                    )
            return ["leap.top", "leap.crd"]
        else:
            raise _ParameterisationError("tLEaP failed!")

    def _run_pdb2gmx(self, molecule, work_dir):
        """
        Run using pdb2gmx.

        Parameters
        ----------

        molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
            The molecule to parameterise.

        work_dir : str
            The working directory.
        """
        from ... import _gmx_exe, _isVerbose
        from ..._Exceptions import ParameterisationError as _ParameterisationError
        from ..._Exceptions import IncompatibleError as _IncompatibleError
        from ... import IO as _IO

        # A list of supported force fields, mapping to their GROMACS ID string.
        # GROMACS supports a sub-set of the AMBER force fields.
        supported_ff = {"ff99": "amber99", "ff99SB": "amber99sb", "ff03": "amber03"}

        if self._forcefield not in supported_ff:
            raise _IncompatibleError(
                "'pdb2gmx' does not support the '%s' force field." % self._forcefield
            )

        # Write the system to a PDB file.
        try:
            _IO.saveMolecules(
                _os.path.join(str(work_dir), "input"),
                molecule,
                "pdb",
                property_map=self._property_map,
            )
        except Exception as e:
            msg = "Failed to write system to 'PDB' format."
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Generate the pdb2gmx command.
        command = (
            "%s pdb2gmx -f input.pdb -o output.gro -p output.top -ignh -ff %s -water none"
            % (_gmx_exe, supported_ff[self._forcefield])
        )

        with open(_os.path.join(str(work_dir), "README.txt"), "w") as file:
            # Write the command to file.
            file.write("# pdb2gmx was run with the following command:\n")
            file.write("%s\n" % command)

        # Create files for stdout/stderr.
        stdout = open(_os.path.join(str(work_dir), "pdb2gmx.out"), "w")
        stderr = open(_os.path.join(str(work_dir), "pdb2gmx.err"), "w")

        # Run pdb2gmx as a subprocess.
        proc = _subprocess.run(
            _Utils.command_split(command),
            cwd=work_dir,
            shell=False,
            stdout=stdout,
            stderr=stderr,
        )
        stdout.close()
        stderr.close()

        # Check for the expected output.
        if _os.path.isfile(
            _os.path.join(str(work_dir), "output.gro")
        ) and _os.path.isfile(_os.path.join(str(work_dir), "output.top")):
            return ["output.gro", "output.top"]
        else:
            raise _ParameterisationError("pdb2gmx failed!")

    @staticmethod
    def _get_disulphide_bonds(
        molecule, tolerance=1.2, max_distance=_Length(6, "A"), property_map={}
    ):
        """
        Internal function to generate LEaP records for disulphide bonds.

        Parameters
        ----------

        molecule : Sire.Mol.Molecule
            The molecule of interest.

        Returns
        -------

        bond_records : [str]
            A list of LEaP formatted bond records.

        tolerance : float
            The tolerance to use when searching for bonds.

        max_distance : :class:`Length <BioSimSpace.Types.Length>`
            The maximum distance between atoms when searching for disulphide
            bonds.

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """
        import warnings as _warnings
        from sire.legacy import Mol as _SireMol

        if not isinstance(molecule, _SireMol.Molecule):
            raise TypeError("'molecule' must be of type 'Sire.Mol.Molecule'")

        if not isinstance(tolerance, (int, float)):
            raise TypeError("'tolerance' must be of type 'float'")
        tolerance = float(tolerance)
        if tolerance < 1:
            raise ValueError("'tolerance' must be >= 1.0.")

        if not isinstance(max_distance, _Length):
            raise ValueError(
                "'max_distance' must be of type 'BioSimSpace.Types.Length'"
            )
        max_radius2 = max_distance.angstroms().value() ** 2

        if not isinstance(property_map, dict):
            raise ValueError("'property_map' must be of type 'dict'")

        # Create a copy of the molecule.
        mol = molecule.__deepcopy__()

        # Get the connectivity of the molecule.
        conn = _SireMol.Connectivity(
            mol, _SireMol.CovalentBondHunter(tolerance, max_radius2), property_map
        )

        # Add this as a molecule property.
        mol = mol.edit().set_property("connectivity", conn).molecule().commit()

        # Create the search query.
        query = _SireMol.Select("bonds from element S to element S")

        # Try searching for disulphide bonds.
        try:
            disulphides = query(mol, property_map).bonds()
        except:
            disulphides = []

        # Create a list to store the LEaP bond record strings.
        bond_records = []

        # Whether a warning has already been raised.
        have_warned = False

        # Loop over the disulphide bonds and generate the records.
        for bond in disulphides:
            # Get the atoms in the bond.
            atom0 = bond.atom0()
            atom1 = bond.atom1()

            # Extract the residue numbers associated with the atoms.
            num0 = atom0.residue().number().value()
            num1 = atom1.residue().number().value()

            # Extract the residue indices associated with the atoms.
            # (Make sure these are one-indexed.)
            idx0 = atom0.residue().index().value() + 1
            idx1 = atom1.residue().index().value() + 1

            # If the residue numbers don't match the indices, then warn the user
            # since the molecule will require renumbering before passing to LEaP.
            if not have_warned:
                if idx0 != num0 or idx1 != num1:
                    _warnings.warn(
                        "Residue numbers are out of order. Use "
                        "'sire.legacy.IO.renumberConstituents' correct numbering."
                    )
                    have_warned = True

            # Extract the atom names.
            name0 = atom0.name().value()
            name1 = atom1.name().value()

            # Create the record.
            record = f"bond mol.{idx0}.{name0} mol.{idx1}.{name1}"

            # Append to the list.
            bond_records.append(record)

        return bond_records

    @staticmethod
    def _generate_bond_records(molecule, bonds):
        """
        Internal function to generate additional LEaP bond records.

        Parameters
        ----------

        molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
            The molecule of interest.

        bonds : ((class:`Atom <BioSimSpace._SireWrappers.Atom>`, class:`Atom <BioSimSpace._SireWrappers.Atom>`))
            An optional tuple of atom pairs to specify additional atoms that
            should be bonded. This is useful when the PDB CONECT record is
            incomplete.

        Returns
        -------

        bond_records : [str]
            A list of LEaP formatted bond records.
        """
        import warnings as _warnings
        from ..._SireWrappers import Atom as _Atom
        from ..._SireWrappers import Molecule as _Molecule

        if bonds is None:
            return []

        if not isinstance(molecule, _Molecule):
            raise TypeError(
                "'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'"
            )

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
                    ) or not (atom1._sire_object.molecule() == molecule._sire_object):
                        raise ValueError(
                            "Atoms in 'bonds' don't belong to the 'molecule'."
                        )

        # Create a list to store the LEaP bond record strings.
        bond_records = []

        # Whether a warning has already been raised.
        have_warned = False

        # Loop over the bonds and generate the records.
        for atom0, atom1 in bonds:
            # Extract the Sire objects.
            atom0 = atom0._sire_object
            atom1 = atom1._sire_object

            # Extract the residue numbers associated with the atoms.
            num0 = atom0.residue().number().value()
            num1 = atom1.residue().number().value()

            # Extract the residue indices associated with the atoms.
            # (Make sure these are one-indexed.)
            idx0 = atom0.residue().index().value() + 1
            idx1 = atom1.residue().index().value() + 1

            # If the residue numbers don't match the indices, then warn the user
            # since the molecule will require renumbering before passing to LEaP.
            if not have_warned:
                if idx0 != num0 or idx1 != num1:
                    _warnings.warn(
                        "Residue numbers are out of order. Use "
                        "'sire.legacy.IO.renumberConstituents' correct numbering."
                    )
                    have_warned = True

            # Extract the atom names.
            name0 = atom0.name().value()
            name1 = atom1.name().value()

            # Create the record.
            record = f"bond mol.{idx0}.{name0} mol.{idx1}.{name1}"

            # Append to the list.
            bond_records.append(record)

        return bond_records


class GAFF(_protocol.Protocol):
    """A class for handling protocols for the GAFF force field model."""

    # A list of supported charge methods.
    _charge_methods = ["RESP", "CM2", "MUL", "BCC", "ESP", "GAS"]

    # Add ABCG2 to the list of charge methods if it is supported.
    # This requires Antechamber > 24.0.
    if _has_abcg2:
        _charge_methods.append("ABCG2")

    def __init__(
        self,
        version,
        charge_method="BCC",
        net_charge=None,
        ensure_compatible=True,
        property_map={},
        **kwargs,
    ):
        """
        Constructor.

        Parameters
        ----------

        version : int
            Whether version 1 or 2 of GAFF.

        charge_method : str
            The method to use when calculating atomic charges:
            "RESP", "CM2", "MUL", "BCC", "ESP", "GAS"

        net_charge : int
            The net charge on the molecule.

        ensure_compatible : bool
            Whether to ensure that the topology of the parameterised molecule is
            compatible with that of the original molecule. An exception will be
            raised if this isn't the case, e.g. if atoms have been added. When
            True, the parameterised molecule will preserve the topology of the
            original molecule, e.g. the original atom and residue names will be
            kept.

            property_map : dict
            A dictionary that maps system "properties" to their user defined

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }

        **kwargs: dict
            Additional keyword arguments. These can be used to pass custom
            parameters to the Antechamber program.
        """
        from ...Types import Charge as _Charge

        if type(version) is not int:
            raise TypeError("'version' must be of type 'int'.")

        # Check that the version is valid.
        if not version in [1, 2]:
            raise ValueError("Unsupported version: options are 1 or 2.")

        if not isinstance(charge_method, str):
            raise TypeError("'charge_method' must be of type 'str'")

        # Strip whitespace and convert to upper case.
        charge_method = charge_method.replace(" ", "").upper()

        # Check that the charge method is valid.
        if not charge_method in self._charge_methods:
            raise ValueError(
                "Unsupported charge method: '%s'. Supported methods are: %s"
                % (charge_method, self._charge_methods)
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

        # Check the kwargs to see whether acdoctor is enabled.
        self._acdoctor = kwargs.get("acdoctor", True)
        if not isinstance(self._acdoctor, bool):
            raise TypeError("'acdoctor' must be of type 'bool'")

        # Set the version.
        self._version = version

        # Set the net molecular charge.
        self._net_charge = net_charge

        # Set the charge method.
        self._charge_method = charge_method

        # Call the base class constructor.
        super().__init__(
            forcefield="gaff",
            ensure_compatible=ensure_compatible,
            property_map=property_map,
        )

    def run(self, molecule, work_dir=None, queue=None):
        """
        Run the parameterisation protocol.

        Parameters
        ----------

        molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`, str
            The molecule to parameterise, either as a Molecule object or SMILES
            string.

        work_dir : :class:`WorkDir <BioSimSpace._Utils.WorkDir>`
            The working directory.

        queue : queue.Queue
            The thread queue is which this method has been run.

        Returns
        -------

        molecule : BioSimSpace._SireWrappers.Molecule
            The parameterised molecule.
        """
        from ..._Exceptions import ParameterisationError as _ParameterisationError
        from ...Convert import smiles as _smiles
        import warnings as _warnings
        from ...Parameters._utils import formalCharge as _formalCharge
        from ... import IO as _IO
        from ..._SireWrappers import Molecule as _Molecule
        from ..._Exceptions import ThirdPartyError as _ThirdPartyError
        import queue as _queue
        from ... import _isVerbose

        if not isinstance(molecule, (_Molecule, str)):
            raise TypeError(
                "'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule' or 'str'"
            )

        if work_dir is not None and not isinstance(work_dir, _Utils.WorkDir):
            raise TypeError("'work_dir' must be of type 'BioSimSpace._Utils.WorkDir'")

        if queue is not None and not isinstance(queue, _queue.Queue):
            raise TypeError("'queue' must be of type 'queue.Queue'")

        # Set work_dir to the current directory.
        if work_dir is None:
            work_dir = _os.getcwd()

        # Try to create a molecule from the SMILES string.
        if isinstance(molecule, str):
            try:
                smiles = molecule
                molecule = _smiles(molecule)
                edit_mol = molecule._sire_object.edit()
                edit_mol = edit_mol.rename(f"smiles:{smiles}").molecule()
                molecule._sire_object = edit_mol.commit()
            except Exception as e:
                msg = "Unable to convert SMILES to Molecule using RDKit."
                if _isVerbose():
                    msg += ": " + getattr(e, "message", repr(e))
                    raise _ThirdPartyError(msg) from e
                else:
                    raise _ThirdPartyError(msg) from None

        # Create a copy of the molecule.
        new_mol = molecule.copy()

        # Use the net molecular charge passed as an option.
        if self._net_charge is not None:
            charge = self._net_charge
        else:
            # The user will likely have passed a bare PDB or Mol2 file.
            # Antechamber expects the molecule to be uncharged, or integer
            # charged (where the charge, or number of electrons, is passed with
            # the -nc flag).

            # Get the total charge on the molecule.
            if "charge" in self._property_map:
                _property_map = {"charge": self._property_map["charge"]}
                prop = self._property_map["charge"]
            else:
                _property_map = {"charge": "charge"}
                prop = "charge"

            # The molecule has a charge property.
            if new_mol._getSireObject().has_property(prop):
                charge = new_mol.charge(property_map=_property_map).value()

                # Charge is non-integer, try to fix it.
                if abs(round(charge) - charge) > 0:
                    new_mol._fixCharge(property_map=_property_map)
                    charge = round(charge)
            else:
                charge = None

            # Only try "formal_charge" when "charge" is missing. Unlikely to have
            # both if this is a bare molecule, but the user could be re-parameterising
            # an existing molecule.
            if charge is None:
                # Get the total formal charge on the molecule.
                if "formal_charge" in self._property_map:
                    _property_map = {"charge": self._property_map["formal_charge"]}
                    prop = self._property_map["charge"]
                else:
                    _property_map = {"charge": "formal_charge"}
                    prop = "formal_charge"

                if new_mol._getSireObject().has_property(prop):
                    charge = new_mol.charge(property_map=_property_map).value()

                    # Compute the formal charge ourselves to check that it is consistent.
                    formal_charge = _formalCharge(new_mol).value()

                    if charge != formal_charge:
                        _warnings.warn(
                            "The formal charge on the molecule is %d "
                            "but we estimate it to be %d" % (charge, formal_charge)
                        )
                else:
                    msg = (
                        "The molecule has no 'charge' or 'formal_charge' information, and "
                        "no 'net_charge' option has been passed. You can use the "
                        "'BioSimSpace.Parameters.formalCharge' function to compute the "
                        "formal charge"
                    )
                    raise _ParameterisationError(msg)

        # If the molecule was loaded from an SDF file, then use this as the
        # intermediate format.
        fileformat_prop = self._property_map.get("fileformat", "fileformat")
        if (
            new_mol._sire_object.has_property(fileformat_prop)
            and "SDF" in new_mol._sire_object.property("fileformat").value()
        ):
            format = "sdf"
        else:
            format = "pdb"

        # Write the system to file.
        try:
            _IO.saveMolecules(
                _os.path.join(str(work_dir), "antechamber"),
                new_mol,
                format,
                property_map=self._property_map,
            )
        except Exception as e:
            msg = "Failed to write system to 'PDB' format."
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Generate the Antechamber command.
        command = (
            "%s -at %d -i antechamber.%s -fi %s "
            + "-o antechamber.mol2 -fo mol2 -c %s -s 2 -nc %d"
        ) % (
            _antechamber_exe,
            self._version,
            format,
            format,
            self._charge_method.lower(),
            charge,
        )

        # Disable acdoctor if requested.
        if not self._acdoctor:
            command += " -dr no"

        with open(_os.path.join(str(work_dir), "README.txt"), "w") as file:
            # Write the command to file.
            file.write("# Antechamber was run with the following command:\n")
            file.write("%s\n" % command)

        # Create files for stdout/stderr.
        stdout = open(_os.path.join(str(work_dir), "antechamber.out"), "w")
        stderr = open(_os.path.join(str(work_dir), "antechamber.err"), "w")

        # Run Antechamber as a subprocess.
        proc = _subprocess.run(
            _Utils.command_split(command),
            cwd=str(work_dir),
            shell=False,
            stdout=stdout,
            stderr=stderr,
        )
        stdout.close()
        stderr.close()

        # Antechamber doesn't return sensible error codes, so we need to check that
        # the expected output was generated.
        if _os.path.isfile(_os.path.join(str(work_dir), "antechamber.mol2")):
            # Run parmchk to check for missing parameters.
            command = (
                "%s -s %d -i antechamber.mol2 -f mol2 " + "-o antechamber.frcmod"
            ) % (_parmchk_exe, self._version)

            with open(_os.path.join(str(work_dir), "README.txt"), "a") as file:
                # Write the command to file.
                file.write("\n# ParmChk was run with the following command:\n")
                file.write("%s\n" % command)

            # Create files for stdout/stderr.
            stdout = open(_os.path.join(str(work_dir), "parmchk.out"), "w")
            stderr = open(_os.path.join(str(work_dir), "parmchk.err"), "w")

            # Run parmchk as a subprocess.
            proc = _subprocess.run(
                _Utils.command_split(command),
                cwd=str(work_dir),
                shell=False,
                stdout=stdout,
                stderr=stderr,
            )
            stdout.close()
            stderr.close()

            # The frcmod file was created.
            if _os.path.isfile(_os.path.join(str(work_dir), "antechamber.frcmod")):
                # Now call tLEaP using the partially parameterised molecule and the frcmod file.
                # tLEap will run in the same working directory, using the Mol2 file generated by
                # Antechamber.

                # Try to find a force field file.
                if self._version == 1:
                    ff = _find_force_field("gaff")
                else:
                    ff = _find_force_field("gaff2")

                # Write the LEaP input file.
                with open(_os.path.join(str(work_dir), "leap.txt"), "w") as file:
                    file.write("source %s\n" % ff)
                    file.write("mol = loadMol2 antechamber.mol2\n")
                    file.write("loadAmberParams antechamber.frcmod\n")
                    file.write("saveAmberParm mol leap.top leap.crd\n")
                    file.write("quit")

                # Generate the tLEaP command.
                command = "%s -f leap.txt" % _tleap_exe

                with open(_os.path.join(str(work_dir), "README.txt"), "a") as file:
                    # Write the command to file.
                    file.write("\n# tLEaP was run with the following command:\n")
                    file.write("%s\n" % command)

                # Create files for stdout/stderr.
                stdout = open(_os.path.join(str(work_dir), "leap.out"), "w")
                stderr = open(_os.path.join(str(work_dir), "leap.err"), "w")

                # Run tLEaP as a subprocess.
                proc = _subprocess.run(
                    _Utils.command_split(command),
                    cwd=str(work_dir),
                    shell=False,
                    stdout=stdout,
                    stderr=stderr,
                )
                stdout.close()
                stderr.close()

                # tLEaP doesn't return sensible error codes, so we need to check that
                # the expected output was generated.
                if _os.path.isfile(
                    _os.path.join(str(work_dir), "leap.top")
                ) and _os.path.isfile(_os.path.join(str(work_dir), "leap.crd")):
                    # Check the output of tLEaP for missing atoms.
                    if self._ensure_compatible:
                        if _has_missing_atoms(_os.path.join(str(work_dir), "leap.out")):
                            raise _ParameterisationError(
                                "tLEaP added missing atoms. The topology is now "
                                "inconsistent with the original molecule. Please "
                                "make sure that your initial molecule has a "
                                "complete topology."
                            )

                    # Load the parameterised molecule. (This could be a system of molecules.)
                    try:
                        par_mol = _IO.readMolecules(
                            [
                                _os.path.join(str(work_dir), "leap.top"),
                                _os.path.join(str(work_dir), "leap.crd"),
                            ],
                        )
                        # Extract single molecules.
                        if par_mol.nMolecules() == 1:
                            par_mol = par_mol.getMolecules()[0]
                    except Exception as e:
                        msg = "Failed to read molecule from: 'leap.top', 'leap.crd'"
                        if _isVerbose():
                            msg += ": " + getattr(e, "message", repr(e))
                            raise IOError(msg) from e
                        else:
                            raise IOError(msg) from None

                    # Make the molecule 'mol' compatible with 'par_mol'. This will create
                    # a mapping between atom indices in the two molecules and add all of
                    # the new properties from 'par_mol' to 'mol'.
                    if self._ensure_compatible:
                        new_mol.makeCompatibleWith(
                            par_mol,
                            property_map=self._property_map,
                            overwrite=True,
                            verbose=False,
                        )
                    else:
                        try:
                            new_mol.makeCompatibleWith(
                                par_mol,
                                property_map=self._property_map,
                                overwrite=True,
                                verbose=False,
                            )
                        except:
                            new_mol = par_mol

                    # Record the forcefield used to parameterise the molecule.
                    new_mol._forcefield = ff

                else:
                    raise _ParameterisationError("tLEaP failed!")
            else:
                raise _ParameterisationError("Parmchk failed!")
        else:
            raise _ParameterisationError("Antechamber failed!")

        if queue is not None:
            queue.put(new_mol)
        return new_mol


def _find_force_field(forcefield):
    """
    Internal helper function to search LEaP compatible force field files.

    Parameters
    ----------

    forcefield : str
        The name of the force field.

    Returns
    -------

    file : str
        The full path of the matching force field file.
    """
    import glob as _glob

    # Whether the force field is old.
    is_old = False

    # Search for a compatible force field file.
    ff = _glob.glob("%s/*.%s" % (_cmd_dir, forcefield))

    # Search the old force fields. First try a specific match.
    if len(ff) == 0:
        ff = _glob.glob("%s/oldff/leaprc.%s" % (_cmd_dir, forcefield))
        is_old = True

        # No matches, try globbing all files with matching extension.
        if len(ff) == 0:
            ff = _glob.glob("%s/oldff/*.%s" % (_cmd_dir, forcefield))

    # No force field found!
    if len(ff) == 0:
        raise ValueError("No force field file found for '%s'" % forcefield)

    # Multiple force fields found.
    elif len(ff) > 1:
        raise ValueError("Multiple force fields found for '%s': %s" % (forcefield, ff))

    # Create the force field name.
    ff = _os.path.basename(ff[0])
    if is_old:
        ff = "oldff/" + ff

    # Return the force field.
    return ff


def _has_missing_atoms(tleap_file):
    """
    Check whether tLEaP has added missing atoms.

    Parameters
    ----------

    tleap_file : str
        The path to the output file generated by tLEaP.

    Returns
    -------

    has_missing_atoms : bool
        Whether missing atoms were added.
    """

    if not isinstance(tleap_file, str):
        raise TypeError("'tleap_file' must be of type 'str'.")

    if not (_os.path.exists(tleap_file)):
        raise IOError(f"'tleap_file' doesn't exist: '{tleap_file}'")

    with open(tleap_file, "r") as file:
        for line in file:
            if "missing atoms" in line:
                return True

    return False
