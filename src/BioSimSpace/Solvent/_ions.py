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

"""Functionality for adding ions to solvated molecular systems."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["addIons", "ions"]

# Map from user-facing ion name (lowercase) to (GROMACS residue name, charge).
# Names and charges are taken from the AMBER03 force field (amber03.ff/ions.itp)
# used as the default topology in BioSimSpace's solvation pipeline.
_ion_map = {
    "li": ("LI", 1),
    "lithium": ("LI", 1),
    "na": ("NA", 1),
    "sodium": ("NA", 1),
    "k": ("K", 1),
    "potassium": ("K", 1),
    "rb": ("RB", 1),
    "rubidium": ("RB", 1),
    "cs": ("CS", 1),
    "cesium": ("CS", 1),
    "cl": ("CL", -1),
    "chloride": ("CL", -1),
    "mg": ("MG", 2),
    "magnesium": ("MG", 2),
    "ca": ("CA", 2),
    "calcium": ("CA", 2),
    "zn": ("ZN", 2),
    "zinc": ("ZN", 2),
}

# Canonical list of supported ion names (GROMACS residue names, lowercase).
_ions = sorted(set(v[0].lower() for v in _ion_map.values()))


def ions():
    """
    Return a list of the supported ion names.

    Returns
    -------

    ions : [str]
        A list of the supported ion names.
    """
    return _ions


def addIons(
    system,
    ion,
    num_ions=0,
    ion_conc=0,
    is_neutral=True,
    preserved_waters=None,
    counter_ion=None,
    work_dir=None,
    property_map={},
):
    """
    Add ions to a pre-solvated molecular system using 'gmx genion'.

    Ions are added by replacing randomly selected water molecules (residue
    name SOL). The system must therefore already contain water.

    Parameters
    ----------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        A pre-solvated molecular system.

    ion : str
        The name of the ion to add. Case-insensitive. Use
        :func:`ions` for a list of supported values, e.g. ``"mg"``,
        ``"ca"``, ``"cl"``.

    num_ions : int
        The number of ions to add. Mutually exclusive with ``ion_conc``.

    ion_conc : float
        The concentration of the requested ion in mol/litre. The number
        of ions to add is calculated from the box volume. This correctly
        accounts for the volume of all molecules already present in the
        system. Mutually exclusive with ``num_ions``. Note that any
        counter-ions added by ``is_neutral=True`` are not included in
        this concentration — their count is determined solely by the net
        system charge.

    is_neutral : bool
        Whether to neutralise the system charge. When ``True``, genion
        adds Na\\ :sup:`+` or Cl\\ :sup:`-` (as appropriate) in addition
        to the requested ion to bring the total system charge to zero.
        Note that existing ions are never removed; neutralisation is
        achieved by adding counter-ions only. The number of counter-ions
        depends on the net system charge and is independent of
        ``ion_conc``.

    preserved_waters : [int] or [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
        A list of water molecules to preserve from replacement by genion.
        Each entry is either an integer (system-level molecule index) or
        a :class:`Molecule <BioSimSpace._SireWrappers.Molecule>` object.
        These molecules are temporarily removed from the system before
        genion runs and re-added to the result afterwards, guaranteeing
        that genion cannot select them for replacement. Useful for
        protecting crystallographic or binding-site water molecules.

    counter_ion : str
        The name of the ion to use for the opposite charge slot in the
        genion call. By default genion uses Na\\ :sup:`+` as the counter
        cation and Cl\\ :sup:`-` as the counter anion. Use this parameter
        to override the counter-ion species, e.g. ``"br"`` when adding a
        cation. The counter-ion must have the opposite charge sign to
        ``ion``. Requires ``is_neutral=True`` (otherwise no counter-ions
        are added and the parameter has no effect).

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values. This allows the user to refer to properties with their
        own naming scheme, e.g. ``{ "charge" : "my-charge" }``.

    Returns
    -------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The molecular system with ions added.
    """
    import math as _math

    from .. import _gmx_exe
    from .._Exceptions import MissingSoftwareError as _MissingSoftwareError

    if _gmx_exe is None:
        raise _MissingSoftwareError(
            "'BioSimSpace.Solvent.addIons' is not supported. "
            "Please install GROMACS (http://www.gromacs.org)."
        )

    from .._SireWrappers import Molecule as _Molecule
    from .._SireWrappers import System as _System

    if not isinstance(system, _System):
        raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")

    if system.nWaterMolecules() == 0:
        raise ValueError(
            "'system' contains no water molecules. 'addIons' requires a "
            "pre-solvated system."
        )

    if not isinstance(ion, str):
        raise TypeError("'ion' must be of type 'str'")

    ion_key = ion.strip().lower()
    if ion_key not in _ion_map:
        raise ValueError(f"Unsupported ion '{ion}'. Supported ions are: {ions()}")
    ion_name, ion_charge = _ion_map[ion_key]

    if not isinstance(num_ions, int):
        raise TypeError("'num_ions' must be of type 'int'")
    if num_ions < 0:
        raise ValueError("'num_ions' cannot be negative!")

    if not isinstance(ion_conc, (int, float)):
        raise TypeError("'ion_conc' must be of type 'float'")
    if ion_conc < 0:
        raise ValueError("'ion_conc' cannot be negative!")

    if num_ions > 0 and ion_conc > 0:
        raise ValueError("'num_ions' and 'ion_conc' are mutually exclusive.")

    if not isinstance(is_neutral, bool):
        raise TypeError("'is_neutral' must be of type 'bool'")

    if num_ions == 0 and ion_conc == 0 and not is_neutral:
        raise ValueError(
            "Nothing to do: 'num_ions' is 0, 'ion_conc' is 0, and 'is_neutral' is "
            "False. Set 'num_ions' > 0, 'ion_conc' > 0, or 'is_neutral' to True."
        )

    # Validate and resolve preserved_waters to a list of Molecule objects.
    preserved_mols = None
    if preserved_waters is not None:
        if not isinstance(preserved_waters, (list, tuple)):
            raise TypeError(
                "'preserved_waters' must be a list of int indices or Molecule objects."
            )
        # Build a set of mol_nums for fast water and membership tests.
        water_nums = {w._sire_object.number() for w in system.getWaterMolecules()}
        n_mols = system.nMolecules()
        preserved_mols = []
        for item in preserved_waters:
            if isinstance(item, int):
                if item < 0 or item >= n_mols:
                    raise ValueError(
                        f"Molecule index {item} is out of range for the system "
                        f"(nMolecules={n_mols})."
                    )
                mol = system[item]
                if mol._sire_object.number() not in water_nums:
                    raise ValueError(
                        f"Molecule at index {item} is not a water molecule."
                    )
                preserved_mols.append(mol)
            elif isinstance(item, _Molecule):
                if item._sire_object.number() not in system._mol_nums:
                    raise ValueError(
                        "A Molecule in 'preserved_waters' does not belong to "
                        "the system."
                    )
                if item._sire_object.number() not in water_nums:
                    raise ValueError(
                        "A Molecule in 'preserved_waters' is not a water molecule."
                    )
                preserved_mols.append(item)
            else:
                raise TypeError(
                    "'preserved_waters' items must be int indices or Molecule objects."
                )

    # Validate and resolve counter_ion.
    counter_ion_name = None
    counter_ion_charge = None
    if counter_ion is not None:
        if not isinstance(counter_ion, str):
            raise TypeError("'counter_ion' must be of type 'str'")
        counter_key = counter_ion.strip().lower()
        if counter_key not in _ion_map:
            raise ValueError(
                f"Unsupported counter_ion '{counter_ion}'. Supported ions are: {ions()}"
            )
        counter_ion_name, counter_ion_charge = _ion_map[counter_key]
        # The counter-ion must have the opposite charge sign to the requested ion.
        if ion_charge > 0 and counter_ion_charge >= 0:
            raise ValueError(
                f"'counter_ion' must be negatively charged when 'ion' is positively "
                f"charged. '{counter_ion}' has charge {counter_ion_charge:+d}."
            )
        if ion_charge < 0 and counter_ion_charge <= 0:
            raise ValueError(
                f"'counter_ion' must be positively charged when 'ion' is negatively "
                f"charged. '{counter_ion}' has charge {counter_ion_charge:+d}."
            )
        if not is_neutral:
            raise ValueError(
                "'counter_ion' has no effect when 'is_neutral=False'. Set "
                "'is_neutral=True' to enable counter-ion addition."
            )

    if work_dir is not None and not isinstance(work_dir, str):
        raise TypeError("'work_dir' must be of type 'str'")

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    # Compute num_ions from the box volume when ion_conc is requested.
    # This calculation uses the full triclinic box volume and does NOT rely on
    # gmx genion's -conc flag, which has a known issue where it ignores the
    # volume occupied by existing ions.
    if ion_conc > 0:
        box, angles = system.getBox(property_map=property_map)
        if box is None:
            raise ValueError(
                "'system' has no periodic box information. Cannot calculate "
                "ion count from 'ion_conc'."
            )
        a = box[0].angstroms().value()
        b = box[1].angstroms().value()
        c = box[2].angstroms().value()
        alpha = angles[0].radians().value()
        beta = angles[1].radians().value()
        gamma = angles[2].radians().value()
        # Triclinic box volume in litres (1 Å³ = 1e-27 L).
        V_L = (
            a
            * b
            * c
            * _math.sqrt(
                1
                - _math.cos(alpha) ** 2
                - _math.cos(beta) ** 2
                - _math.cos(gamma) ** 2
                + 2 * _math.cos(alpha) * _math.cos(beta) * _math.cos(gamma)
            )
        ) * 1e-27
        num_ions = max(1, round(ion_conc * V_L * 6.02214076e23))

    return _add_ions(
        system,
        ion_name,
        ion_charge,
        num_ions,
        is_neutral,
        preserved_mols,
        counter_ion_name,
        counter_ion_charge,
        work_dir,
        property_map,
    )


def _add_ions(
    system,
    ion_name,
    ion_charge,
    num_ions,
    is_neutral,
    preserved_mols=None,
    counter_ion_name=None,
    counter_ion_charge=None,
    work_dir=None,
    property_map={},
):
    """
    Internal function to add ions to a solvated system using 'gmx genion'.

    Parameters
    ----------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The pre-solvated molecular system.

    ion_name : str
        The GROMACS residue name of the ion (e.g. ``"MG"``, ``"CL"``).

    ion_charge : int
        The integer charge of the ion (e.g. ``2``, ``-1``).

    num_ions : int
        The number of ions to add.

    is_neutral : bool
        Whether to neutralise the system charge.

    preserved_mols : [:class:`Molecule <BioSimSpace._SireWrappers.Molecule>`]
        Molecules to exclude from genion's replacement pool. They are
        removed from the working copy before genion runs and reinserted
        between the non-water solute and the new water+ions in the result,
        mirroring the ordering used by ``_solvate`` for crystal waters.

    counter_ion_name : str
        The GROMACS residue name of the counter-ion (e.g. ``"BR"``), or
        ``None`` to use genion's default (NA/CL).

    counter_ion_charge : int
        The integer charge of the counter-ion, or ``None``.

    work_dir : str
        The working directory for the process.

    property_map : dict
        A dictionary that maps system "properties" to their user defined
        values.

    Returns
    -------

    system : :class:`System <BioSimSpace._SireWrappers.System>`
        The molecular system with ions added.
    """
    import os as _os
    import shutil as _shutil
    import subprocess as _subprocess

    from sire.legacy import Base as _SireBase

    from .. import IO as _IO
    from .. import _gmx_exe, _isVerbose, _Utils
    from .._SireWrappers import System as _System

    # Work on a copy so we don't modify the caller's system.
    _system = _System(system)

    # Convert all water molecules to GROMACS topology (SOL naming) so that
    # genion can identify and replace them. BioSimSpace switches water topology
    # depending on the backend engine; this makes it explicit for GROMACS.
    _system._set_water_topology("GROMACS", property_map=property_map)

    # The number of atoms per water molecule, needed to convert atom counts
    # to molecule counts when writing the topology file.
    num_point = _system.getWaterMolecules()[0].nAtoms()

    # Detect the water model name for the topology include.
    try:
        water_model = system._sire_object.property("water_model").to_string().lower()
    except Exception:
        if num_point == 3:
            water_model = "tip3p"
        elif num_point == 4:
            water_model = "tip4p"
        elif num_point == 5:
            water_model = "tip5p"
        else:
            water_model = "tip3p"

    # Remove preserved water molecules from the working copy before rearranging.
    # They are excluded from genion's replacement pool entirely and will be
    # reinserted into the result after the genion run.
    if preserved_mols:
        _system.removeMolecules(preserved_mols)

    # genion requires water molecules to be contiguous at the end of the
    # topology. Remove water, record how many non-water (solute) atoms there
    # are, then re-append the water. This also preserves the original ordering
    # of non-water molecules.
    waters = _system.getWaterMolecules()
    _system.removeWaterMolecules()
    num_solute_atoms = _system.nAtoms()
    _system = _system + waters

    # Create the working directory.
    work_dir = _Utils.WorkDir(work_dir)

    # Write to 6dp unless precision is already specified.
    _property_map = property_map.copy()
    if "precision" not in _property_map:
        _property_map["precision"] = _SireBase.wrap(6)

    with _Utils.cd(work_dir):
        # Write the rearranged system to GROMACS GRO and TOP files.
        try:
            _IO.saveMolecules(
                "system",
                _system,
                "gro87",
                match_water=False,
                property_map=_property_map,
            )
            _IO.saveMolecules(
                "system",
                _system,
                "grotop",
                match_water=False,
                property_map=_property_map,
            )
        except Exception as e:
            msg = (
                "Failed to write GROMACS topology file. Is your molecule parameterised?"
            )
            if _isVerbose():
                raise IOError(msg) from e
            else:
                raise IOError(msg) from None

        # Ensure amber03.ff/ions.itp is included in the TOP so that grompp
        # can resolve the ion atom types. Insert after the last #include if
        # missing.
        with open("system.top", "r") as file:
            top_lines = file.readlines()

        if not any("ions.itp" in line for line in top_lines):
            last_include = -1
            for i, line in enumerate(top_lines):
                if line.strip().startswith("#include"):
                    last_include = i
            if last_include >= 0:
                top_lines.insert(
                    last_include + 1,
                    '; Include ions\n#include "amber03.ff/ions.itp"\n\n',
                )
                with open("system.top", "w") as file:
                    file.writelines(top_lines)

        # Write a minimal mdp file required by grompp.
        with open("ions.mdp", "w") as file:
            file.write("; Neighbour searching\n")
            file.write("cutoff-scheme           = Verlet\n")
            file.write("rlist                   = 1.1\n")
            file.write("pbc                     = xyz\n")
            file.write("verlet-buffer-tolerance = -1\n")
            file.write("\n; Electrostatics\n")
            file.write("coulombtype             = cut-off\n")
            file.write("\n; VdW\n")
            file.write("rvdw                    = 1.0\n")

        # Create the grompp command.
        command = (
            "%s grompp -f ions.mdp -po ions.out.mdp -c system.gro -p system.top -o ions.tpr"
            % _gmx_exe
        )

        with open("README.txt", "w") as file:
            file.write("# gmx grompp was run with the following command:\n")
            file.write("%s\n" % command)

        # Create files for stdout/stderr.
        stdout = open("grompp.out", "w")
        stderr = open("grompp.err", "w")

        # Run grompp as a subprocess.
        _subprocess.run(
            _Utils.command_split(command), shell=False, stdout=stdout, stderr=stderr
        )
        stdout.close()
        stderr.close()

        # Check for the tpr output file.
        if not _os.path.isfile("ions.tpr"):
            raise RuntimeError(
                "'gmx grompp' failed to generate the required output for "
                "'gmx genion'. Check grompp.err for details."
            )

        # Build the genion command.
        # genion supports one positive ion type (-pname/-pq/-np) and one
        # negative ion type (-nname/-nq/-nn). We use the appropriate slot
        # for the requested ion. When is_neutral=True, genion adds counter-ions
        # to bring the total charge to zero; counter_ion_name overrides the
        # default counter-ion species (NA for cations, CL for anions).
        if ion_charge > 0:
            command = (
                "%s genion -s ions.tpr -o ions_out.gro -p system.top"
                " -pname %s -pq %d" % (_gmx_exe, ion_name, ion_charge)
            )
            if num_ions > 0:
                command += " -np %d" % num_ions
            # Override the default Cl- counter-ion if requested.
            if counter_ion_name is not None:
                command += " -nname %s -nq %d" % (counter_ion_name, counter_ion_charge)
        else:
            command = (
                "%s genion -s ions.tpr -o ions_out.gro -p system.top"
                " -nname %s -nq %d" % (_gmx_exe, ion_name, ion_charge)
            )
            if num_ions > 0:
                command += " -nn %d" % num_ions
            # Override the default Na+ counter-ion if requested.
            if counter_ion_name is not None:
                command += " -pname %s -pq %d" % (counter_ion_name, counter_ion_charge)

        if is_neutral:
            command += " -neutral"
        else:
            command += " -noneutral"

        with open("README.txt", "a") as file:
            file.write("\n# gmx genion was run with the following command:\n")
            file.write("%s\n" % command)

        # Create files for stdout/stderr.
        stdout = open("genion.out", "w")
        stderr = open("genion.err", "w")

        # Run genion as a subprocess. genion reads the group to replace
        # interactively; pipe "SOL" to select the solvent group.
        proc_echo = _subprocess.Popen(
            ["echo", "SOL"], shell=False, stdout=_subprocess.PIPE
        )
        proc = _subprocess.Popen(
            _Utils.command_split(command),
            shell=False,
            stdin=proc_echo.stdout,
            stdout=stdout,
            stderr=stderr,
        )
        proc.wait()
        proc_echo.stdout.close()
        stdout.close()
        stderr.close()

        # Check for the output GRO file.
        if not _os.path.isfile("ions_out.gro"):
            raise RuntimeError(
                "'gmx genion' failed to add ions! Check genion.err for details."
            )

        # Parse the output GRO file to extract the water and ion atoms.
        # We skip the first num_solute_atoms lines (which belong to non-water
        # molecules) and collect everything after that. The residue name is
        # read from the fixed GRO column layout (characters 5-10) to track
        # all ion types generically, including any counter-ions added by
        # -neutral.
        water_ion_lines = []
        ion_counts = {}  # residue name -> atom count

        with open("ions_out.gro", "r") as file:
            gro_lines = file.readlines()

        for line in gro_lines[num_solute_atoms + 2 : -1]:
            resname = line[5:10].strip()
            water_ion_lines.append(line)
            ion_counts[resname] = ion_counts.get(resname, 0) + 1

        # Add the box information (last line of the GRO file).
        water_ion_lines.append(gro_lines[-1])

        # Write a GRO file containing only the water and ion atoms.
        if len(water_ion_lines) - 1 > 0:
            with open("water_ions.gro", "w") as file:
                file.write("BioSimSpace %s water box\n" % water_model.upper())
                file.write("%d\n" % (len(water_ion_lines) - 1))
                for line in water_ion_lines:
                    file.write("%s" % line)

        # For OPC water, copy and strip the local topology file so it can
        # be included without its own [defaults] section.
        if water_model == "opc":
            template = _SireBase.getShareDir() + "/templates/water/opc"
            _shutil.copyfile(template + ".itp", "opc.top")
            with open("opc.top", "r") as file:
                opc_lines = file.readlines()
            with open("opc.top", "w") as file:
                for line in opc_lines[6:-1]:
                    file.write(line)

        # Write a TOP file for the water+ions subsystem with the updated
        # molecule counts. SOL atom counts are divided by num_point to give
        # the number of water molecules; all other residues are single-atom
        # ions and their counts are used directly.
        with open("water_ions.top", "w") as file:
            file.write("#define FLEXIBLE 1\n\n")
            file.write("; Include AmberO3 force field\n")
            file.write('#include "amber03.ff/forcefield.itp"\n\n')
            file.write("; Include %s water topology\n" % water_model.upper())
            if water_model == "opc":
                file.write('#include "opc.top"\n\n')
            else:
                file.write('#include "amber03.ff/%s.itp"\n\n' % water_model)
            file.write("; Include ions\n")
            file.write('#include "amber03.ff/ions.itp"\n\n')
            file.write("[ system ] \n")
            file.write("BioSimSpace %s water box\n\n" % water_model.upper())
            file.write("[ molecules ] \n")
            file.write(";molecule name    nr.\n")
            for resname, count in ion_counts.items():
                if resname == "SOL":
                    file.write("SOL               %d\n" % (count / num_point))
                else:
                    file.write("%-18s%d\n" % (resname, count))

        # Load the water+ions subsystem.
        water_ions = _IO.readMolecules(["water_ions.gro", "water_ions.top"])

        # Reconstruct the full system: take the non-water molecules from the
        # original system (preserving their original ordering and properties),
        # then insert any preserved waters, then append the new water+ions.
        # This mirrors _solvate's: molecule + crystal_waters + water_ions.
        solute = _System(system)
        solute.removeWaterMolecules()
        if preserved_mols:
            for mol in preserved_mols:
                solute = solute + mol
        result = solute + water_ions

        # Propagate space and other properties from the water+ions subsystem
        # to the result (this carries the box information from the genion run).
        for prop in water_ions._sire_object.property_keys():
            prop = _property_map.get(prop, prop)
            result._sire_object.set_property(
                prop, water_ions._sire_object.property(prop)
            )

        # Preserve the water model property.
        try:
            wm = system._sire_object.property("water_model")
            result._sire_object.set_property("water_model", wm)
        except Exception:
            pass

    return result
