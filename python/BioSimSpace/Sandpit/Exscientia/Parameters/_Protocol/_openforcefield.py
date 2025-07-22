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
for force fields from the Open Force Field Initiative.
Author: Lester Hedges <lester.hedges@gmail.com>.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["OpenForceField"]

# To override any protocols, just implement a custom "run" method in any
# of the classes.

from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from ..._Utils import _try_import, _have_imported

import os as _os

import queue as _queue
import subprocess as _subprocess

import warnings as _warnings

# Suppress duplicate to-Python converted warnings.
# Both Sire and RDKit register the same converter.
with _warnings.catch_warnings():
    _warnings.simplefilter("ignore")
    _rdkit = _try_import("rdkit")

    if _have_imported(_rdkit):
        from rdkit import Chem as _Chem
        from rdkit import RDLogger as _RDLogger

        # Disable RDKit warnings.
        _RDLogger.DisableLog("rdApp.*")
    else:
        _Chem = _rdkit
        _RDLogger = _rdkit

import sys as _sys

# Temporarily redirect stderr to suppress import warnings.
_orig_stderr = _sys.stderr
_sys.stderr = open(_os.devnull, "w")

_openmm = _try_import("openmm")

if _have_imported(_openmm):
    from openmm.app import PDBFile as _PDBFile
else:
    _PDBFile = _openmm

_openff = _try_import("openff")

# Initialise the NAGL support flag.
_has_nagl = False

if _have_imported(_openff):
    from openff.interchange import Interchange as _Interchange
    from openff.toolkit.topology import Molecule as _OpenFFMolecule
    from openff.toolkit.typing.engines.smirnoff import ForceField as _Forcefield

    try:
        from openff.toolkit.utils.nagl_wrapper import (
            NAGLToolkitWrapper as _NAGLToolkitWrapper,
        )

        _has_nagl = _NAGLToolkitWrapper.is_available()
        from openff.nagl_models import get_models_by_type as _get_models_by_type

        _models = _get_models_by_type("am1bcc")
        try:
            # Find the most recent AM1-BCC release candidate.
            _nagl = _NAGLToolkitWrapper()
            _nagl_model = sorted(
                [str(model) for model in _models if "rc" in str(model)], reverse=True
            )[0]
        except:
            _has_nagl = False
        del _models
    except:
        _has_nagl = False
else:
    _Interchange = _openff
    _OpenFFMolecule = _openff
    _Forcefield = _openff


try:
    from sire.legacy.Base import findExe as _findExe

    _findExe("antechamber")
    _has_antechamber = True
except:
    _has_antechamber = False

# Reset stderr.
_sys.stderr = _orig_stderr
del _sys, _orig_stderr

from sire.legacy import IO as _SireIO
from sire.legacy import Mol as _SireMol
from sire.legacy import System as _SireSystem

from ... import _isVerbose
from ... import Convert as _Convert
from ... import IO as _IO
from ..._Exceptions import ConversionError as _ConversionError
from ..._Exceptions import IncompatibleError as _IncompatibleError
from ..._Exceptions import ThirdPartyError as _ThirdPartyError
from ..._SireWrappers import Molecule as _Molecule
from ... import _Utils

from . import _protocol


class OpenForceField(_protocol.Protocol):
    """A class for handling protocols for Open Force Field models."""

    def __init__(
        self, forcefield, ensure_compatible=True, use_nagl=True, property_map={}
    ):
        """
        Constructor.

        Parameters
        ----------

        forcefield : str
            The name of the force field.

        ensure_compatible : bool
            Whether to ensure that the topology of the parameterised molecule is
            compatible with that of the original molecule. An exception will be
            raised if this isn't the case, e.g. if atoms have been added. When
            True, the parameterised molecule will preserve the topology of the
            original molecule, e.g. the original atom and residue names will be
            kept.

        use_nagl : bool
            Whether to use NAGL to compute AM1-BCC charges. If False, the default
            is to use AmberTools via antechamber and sqm. (This option is only
            used if NAGL is available.)

        property_map : dict
            A dictionary that maps system "properties" to their user defined
            values. This allows the user to refer to properties with their
            own naming scheme, e.g. { "charge" : "my-charge" }
        """

        # Call the base class constructor.
        super().__init__(
            forcefield=forcefield,
            ensure_compatible=ensure_compatible,
            property_map=property_map,
        )

        if not isinstance(use_nagl, bool):
            raise TypeError("'use_nagl' must be of type 'bool'")

        # Set the NAGL flag.
        self._use_nagl = use_nagl

        # Set the compatibility flags.
        self._tleap = False
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

        # Try to create BioSimSpace molecule from the SMILES string.
        if isinstance(molecule, str):
            is_smiles = True
            try:
                smiles = molecule
                molecule = _Convert.smiles(molecule)
            except Exception as e:
                msg = "Unable to convert SMILES to Molecule using RDKit."
                if _isVerbose():
                    msg += ": " + getattr(e, "message", repr(e))
                    raise _ThirdPartyError(msg) from e
                else:
                    raise _ThirdPartyError(msg) from None
        else:
            is_smiles = False

        # Try converting to RDKit format.
        try:
            rdmol = _Convert.toRDKit(molecule, property_map=self._property_map)
        except Exception as e:
            msg = "Failed to convert molecule to RDKit format."
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise (msg) from e
            else:
                raise _ConversionError(msg) from None

        # Create the Open Forcefield Molecule from the RDKit molecule.
        try:
            off_molecule = _OpenFFMolecule.from_rdkit(rdmol)
        except Exception as e:
            msg = "Unable to create OpenFF Molecule!"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _ThirdPartyError(msg) from e
            else:
                raise _ThirdPartyError(msg) from None

        # Apply AM1-BCC charges using NAGL.
        if _has_nagl and self._use_nagl:
            try:
                _nagl.assign_partial_charges(
                    off_molecule, partial_charge_method=_nagl_model
                )
            except Exception as e:
                msg = "Failed to assign AM1-BCC charges using NAGL."
                if _isVerbose():
                    msg += ": " + getattr(e, "message", repr(e))
                    raise _ThirdPartyError(msg) from e
                else:
                    raise _ThirdPartyError(msg) from None
            charge_from_molecules = [off_molecule]
        else:
            if not _has_antechamber:
                raise _MissingSoftwareError(
                    f"'{forcefield}' is not supported. AmberTools "
                    "(http://ambermd.org) is needed for charge "
                    "calculation and 'antechamber' executable "
                    "must be in your PATH."
                ) from None
            charge_from_molecules = None

        # Extract the molecular topology.
        try:
            off_topology = off_molecule.to_topology()
        except Exception as e:
            msg = "Unable to create OpenFF Topology!"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _ThirdPartyError(msg) from e
            else:
                raise _ThirdPartyError(msg) from None

        # Load the force field.
        try:
            ff = self._forcefield + ".offxml"
            forcefield = _Forcefield(ff)
        except Exception as e:
            msg = f"Unable to load force field: {ff}"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _ThirdPartyError(msg) from e
            else:
                raise _ThirdPartyError(msg) from None

        # Create an Interchange object.
        try:
            interchange = _Interchange.from_smirnoff(
                force_field=forcefield,
                topology=off_topology,
                charge_from_molecules=charge_from_molecules,
            )
        except Exception as e:
            msg = "Unable to create OpenFF Interchange object!"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _ThirdPartyError(msg) from e
            else:
                raise _ThirdPartyError(msg) from None

        # Export AMBER format files.
        try:
            interchange.to_prmtop(_os.path.join(str(work_dir), "interchange.prm7"))
            interchange.to_inpcrd(_os.path.join(str(work_dir), "interchange.rst7"))
        except Exception as e:
            msg = "Unable to write Interchange object to AMBER format!"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _ThirdPartyError(msg) from e
            else:
                raise _ThirdPartyError(msg) from None

        # Load the parameterised molecule. (This could be a system of molecules.)
        try:
            par_mol = _IO.readMolecules(
                [
                    _os.path.join(str(work_dir), "interchange.prm7"),
                    _os.path.join(str(work_dir), "interchange.rst7"),
                ],
            )
            # Extract single molecules.
            if par_mol.nMolecules() == 1:
                par_mol = par_mol.getMolecules()[0]
        except Exception as e:
            msg = "Failed to read molecule from: 'interchange.prm7', 'interchange.rst7'"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Make sure we retain stereochemistry information from the SMILES string.
        if is_smiles:
            new_mol = _Convert.smiles(smiles)
            edit_mol = new_mol._sire_object.edit()
            edit_mol = edit_mol.rename(f"smiles:{smiles}").molecule()
            new_mol._sire_object = edit_mol.commit()

        # Make the parameterised molecule compatible with the original topology.
        if self._ensure_compatible:
            new_mol = molecule.copy()
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
