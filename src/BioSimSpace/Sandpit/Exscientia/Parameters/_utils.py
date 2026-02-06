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

"""Utility functions."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["formalCharge"]


def formalCharge(molecule, property_map={}):
    """
    Compute the formal charge on a molecule. This function requires that
    the molecule has explicit hydrogen atoms.

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

    formal_charge : :class:`Charge <BioSimSpace.Types.Charge>`
        The total formal charge on the molecule.
    """
    from .. import IO as _IO
    import tempfile as _tempfile
    from ..Units.Charge import electron_charge as _electron_charge
    from .._SireWrappers import Molecule as _Molecule
    from .. import _Utils
    from .. import _isVerbose

    if not isinstance(molecule, _Molecule):
        raise TypeError(
            "'molecule' must be of type 'BioSimSpace._SireWrappers.Molecule'"
        )

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    from rdkit import Chem as _Chem
    from rdkit import RDLogger as _RDLogger

    # Disable RDKit warnings.
    _RDLogger.DisableLog("rdApp.*")

    # Create a temporary working directory.
    tmp_dir = _tempfile.TemporaryDirectory()
    work_dir = tmp_dir.name

    # Zero the total formal charge.
    formal_charge = 0

    # Get the fileformat property name.
    property = property_map.get("fileformat", "fileformat")

    # Preferentially use the file format that the molecule was loaded from.
    try:
        # Get the raw list of formats.
        raw_formats = molecule._sire_object.property(property).value().split(",")

        # Remove all formats other than PDB and SDF.
        formats = [f for f in raw_formats if f in ["PDB", "SDF"]]

        if len(formats) == 0:
            formats = ["PDB", "SDF"]
        elif len(formats) == 1:
            if formats[0] == "PDB":
                formats.append("SDF")
            else:
                formats.append("PDB")
    except:
        formats = ["PDB", "SDF"]

    # List of exceptions.
    exceptions = []

    # Try going via each format in turn.
    for format in formats:
        try:
            with _Utils.cd(work_dir):
                # Save the molecule in the given format.
                _IO.saveMolecules("tmp", molecule, format)

                # Load with RDKit.
                if format == "SDF":
                    rdmol = _Chem.MolFromMolFile("tmp.sdf")
                else:
                    rdmol = _Chem.MolFromPDBFile("tmp.pdb")

                # Compute the formal charge.
                formal_charge = _Chem.rdmolops.GetFormalCharge(rdmol)

                return formal_charge * _electron_charge

        except Exception as e:
            exceptions.append(e)

    # If we got this far, then we failed to compute the formal charge.
    msg = "Failed to compute the formal charge on the molecule."
    if _isVerbose():
        for e in exceptions:
            msg += "\n\n" + str(e)
        raise RuntimeError(msg)
    else:
        raise RuntimeError(msg) from None
