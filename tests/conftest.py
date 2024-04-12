collect_ignore_glob = ["*/out_test*.py"]

import os
import pytest

from pathlib import Path

import BioSimSpace as BSS

from BioSimSpace._Utils import _try_import, _have_imported

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

# Check for MDAnalysis.
mda = _try_import("MDAnalysis")
has_mdanalysis = _have_imported(mda)

# Check for MDTraj.
mdtraj = _try_import("mdtraj")
has_mdtraj = _have_imported(mdtraj)

# Check for alchemlyb.
_alchemlyb = _try_import("alchemlyb")
has_alchemlyb = _have_imported(_alchemlyb)

# Allow tests to be run from any directory.
root_fp = Path(__file__).parent.resolve()

# Fixtures for tests.


@pytest.fixture(scope="session")
def system():
    """Solvated alanine dipeptide system."""
    return BSS.IO.readMolecules(["tests/input/ala.top", "tests/input/ala.crd"])


@pytest.fixture(scope="module")
def perturbable_system():
    """A vacuum perturbable system."""
    return BSS.IO.readPerturbableSystem(
        f"{url}/perturbable_system0.prm7",
        f"{url}/perturbable_system0.rst7",
        f"{url}/perturbable_system1.prm7",
        f"{url}/perturbable_system1.rst7",
    )


@pytest.fixture(scope="module")
def solvated_perturbable_system():
    """A solvated perturbable system."""
    return BSS.IO.readPerturbableSystem(
        f"{url}/solvated_perturbable_system0.prm7",
        f"{url}/solvated_perturbable_system0.rst7",
        f"{url}/solvated_perturbable_system1.prm7",
        f"{url}/solvated_perturbable_system1.rst7",
    )
