from typing import Dict, Optional

import os
import pandas
import pytest

from BioSimSpace.Sandpit import Exscientia as BSS

from BioSimSpace.Sandpit.Exscientia.Process import Gromacs
from BioSimSpace.Sandpit.Exscientia.Protocol import Production
from BioSimSpace.Sandpit.Exscientia.Units.Temperature import kelvin
from BioSimSpace.Sandpit.Exscientia.Units.Time import femtosecond
from BioSimSpace.Sandpit.Exscientia._SireWrappers import Molecule

from BioSimSpace.Sandpit.Exscientia._Utils import _try_import, _have_imported

# Store the tutorial URL.
url = BSS.tutorialUrl()

# Make sure GROMACS is installed.
has_gromacs = BSS._gmx_exe is not None

# Make sure AMBER is installed.
if BSS._amber_home is not None:
    exe = "%s/bin/sander" % BSS._amber_home
    if os.path.isfile(exe):
        has_amber = True
    else:
        has_amber = False
else:
    has_amber = False

# Make sure NAMD is installed.
try:
    from sire.legacy.Base import findExe

    findExe("namd2")
    has_namd = True
except:
    has_namd = False

# Check whether AMBER parameterisation executables are installed.
has_tleap = BSS.Parameters._Protocol._amber._tleap_exe is not None
has_antechamber = BSS.Parameters._Protocol._amber._antechamber_exe is not None

# Check if openff is installed.
_openff = _try_import("openff")
has_openff = _have_imported(_openff)

# Make sure pyarrow is available as the pandas parquet engine. The parquet
# code does not work with fastparquet.
try:
    pandas.io.parquet.get_engine("pyarrow")
    has_pyarrow = True
except:
    has_pyarrow = False

# Check for alchemlyb.
try:
    import alchemlyb
    from alchemlyb.parsing.gmx import extract_u_nk

    has_alchemlyb = True

    # Check for parquet support.
    major, minor, *_ = alchemlyb.__version__.split(".")
    major = int(major)
    minor = int(minor)
    if major < 2 or (major < 3 and minor < 1):
        has_alchemlyb_parquet = False
    else:
        has_alchemlyb_parquet = True

except ModuleNotFoundError:
    has_alchemlyb = False
    has_alchemlyb_parquet = False

# Check for alchemtest.
try:
    from alchemtest.gmx import load_ABFE
    from alchemtest.amber import load_bace_example

    has_alchemtest = True
except ModuleNotFoundError:
    has_alchemtest = False

# Check for MDAnalysis
mda = _try_import("MDAnalysis")
has_mdanalysis = _have_imported(mda)

# Check for MDTraj.
mdtraj = _try_import("mdtraj")
has_mdtraj = _have_imported(mdtraj)


# Check for MDRestraintsGenerator.
_MDRestraintsGenerator = _try_import(
    "MDRestraintsGenerator",
    install_command="pip install MDRestraintsGenerator",
)

has_mdrestraints_generator = _have_imported(_MDRestraintsGenerator)


def get_energy(
    mol: Molecule,
    coord_mol: Optional[Molecule] = None,
    coord_mol_mapping: Optional[Dict[int, int]] = None,
) -> float:
    if coord_mol is not None:
        c = mol._sire_object.cursor()
        c_coord = coord_mol._sire_object.cursor()
        coordinates = c_coord["coordinates"].to_vector()
        if coord_mol_mapping is not None:
            sorted_keys = sorted(coord_mol_mapping)
            coordinates = [coordinates[coord_mol_mapping[k]] for k in sorted_keys]
        for atom, coord in zip(c.atoms(), coordinates):
            atom["coordinates"] = coord
        mol = Molecule(c.commit())

    system = mol.toSystem()

    protocol = Production(
        timestep=0.1 * femtosecond,
        runtime=0.1 * femtosecond,
        temperature=0 * kelvin,
    )

    proc = Gromacs(system, protocol)
    proc.start()
    proc.wait()

    res = proc.getPotentialEnergy().value()

    return res


@pytest.fixture
def benzene():
    return BSS.Parameters.parameterise(
        "c1ccccc1", "openff_unconstrained-2.0.0"
    ).getMolecule()


@pytest.fixture
def pyrrole():
    return BSS.Parameters.parameterise(
        "C1=CNC=C1", "openff_unconstrained-2.0.0"
    ).getMolecule()


@pytest.fixture
def mapping_benzene_pyrrole(benzene, pyrrole):
    return BSS.Align.matchAtoms(benzene, pyrrole, complete_rings_only=False)


@pytest.fixture
def merged_benzene_pyrrole(benzene, pyrrole, mapping_benzene_pyrrole):
    return BSS.Align.merge(
        benzene,
        pyrrole,
        mapping=mapping_benzene_pyrrole,
        allow_ring_breaking=True,
        allow_ring_size_change=True,
    )
