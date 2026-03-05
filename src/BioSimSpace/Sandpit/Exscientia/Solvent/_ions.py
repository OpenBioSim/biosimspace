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
    "f": ("F", -1),
    "fluoride": ("F", -1),
    "cl": ("CL", -1),
    "chloride": ("CL", -1),
    "br": ("BR", -1),
    "bromide": ("BR", -1),
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


def addIons(system, ion, num_ions=0, is_neutral=True, work_dir=None, property_map={}):
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
        The number of ions to add.

    is_neutral : bool
        Whether to neutralise the system charge. When ``True``, genion
        adds Na\\ :sup:`+` or Cl\\ :sup:`-` (as appropriate) in addition
        to the requested ion to bring the total system charge to zero.
        Note that existing ions are never removed; neutralisation is
        achieved by adding counter-ions only.

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
    from .. import _gmx_exe
    from .._Exceptions import MissingSoftwareError as _MissingSoftwareError

    if _gmx_exe is None:
        raise _MissingSoftwareError(
            "'BioSimSpace.Solvent.addIons' is not supported. "
            "Please install GROMACS (http://www.gromacs.org)."
        )

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

    if not isinstance(is_neutral, bool):
        raise TypeError("'is_neutral' must be of type 'bool'")

    if num_ions == 0 and not is_neutral:
        raise ValueError(
            "Nothing to do: 'num_ions' is 0 and 'is_neutral' is False. "
            "Set 'num_ions' > 0 or 'is_neutral' to True."
        )

    if work_dir is not None and not isinstance(work_dir, str):
        raise TypeError("'work_dir' must be of type 'str'")

    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    return _add_ions(
        system, ion_name, ion_charge, num_ions, is_neutral, work_dir, property_map
    )


def _add_ions(
    system,
    ion_name,
    ion_charge,
    num_ions,
    is_neutral,
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
        # for the requested ion. When is_neutral=True, genion adds the
        # default NA/CL as counter-ions to bring the total charge to zero.
        if ion_charge > 0:
            command = (
                "%s genion -s ions.tpr -o ions_out.gro -p system.top"
                " -pname %s -pq %d" % (_gmx_exe, ion_name, ion_charge)
            )
            if num_ions > 0:
                command += " -np %d" % num_ions
        else:
            command = (
                "%s genion -s ions.tpr -o ions_out.gro -p system.top"
                " -nname %s -nq %d" % (_gmx_exe, ion_name, ion_charge)
            )
            if num_ions > 0:
                command += " -nn %d" % num_ions

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
        # original system (preserving their original ordering and properties)
        # and append the new water+ions.
        solute = _System(system)
        solute.removeWaterMolecules()
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
