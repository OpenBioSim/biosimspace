import math
import pytest
import requests
import tarfile
import tempfile

import BioSimSpace as BSS

from tests.conftest import url, has_alchemlyb, has_gromacs


@pytest.fixture(scope="module")
def fep_output():
    """Path to a temporary directory containing FEP output."""

    # Create URL to test test data.
    fep_url = f"{url}/fep_output.tar.gz"

    # Create a temporary directory.
    tmp_dir = tempfile.TemporaryDirectory()

    # Create the name of the local file.
    local_file = tmp_dir.name + "/fep_output"

    # Download the test data.
    req = requests.get(fep_url, stream=True)
    with open(local_file + ".tar.gz", "wb") as f:
        for chunk in req.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)
                f.flush()

    # Now extract the archive.
    tar = tarfile.open(local_file + ".tar.gz")
    tar.extractall(tmp_dir.name)
    tar.close()

    return tmp_dir


@pytest.fixture(scope="module")
def expected_results():
    """A dictionary of expected FEP results."""

    return {
        "amber": {"mbar": -12.5939, "ti": -13.6850},
        "gromacs": {"mbar": -6.0238, "ti": -8.4158},
        "somd": {"mbar": -6.3519, "ti": -6.3209},
    }


def test_setup_somd(perturbable_system):
    """Test setup for a relative alchemical free energy leg using SOMD."""
    free_nrg = BSS.FreeEnergy.Relative(perturbable_system, engine="somd")


@pytest.mark.skipif(has_gromacs is False, reason="Requires GROMACS to be installed.")
def test_setup_gromacs(perturbable_system):
    """Test setup for a relative alchemical free energy leg using GROMACS."""
    free_nrg = BSS.FreeEnergy.Relative(perturbable_system, engine="gromacs")


@pytest.mark.skipif(
    has_alchemlyb is False, reason="Requires alchemlyb to be installed."
)
@pytest.mark.parametrize("engine", ["amber", "gromacs", "somd"])
@pytest.mark.parametrize("estimator", ["mbar", "ti"])
def test_analysis(fep_output, engine, estimator, expected_results):
    """Test that the free energy analysis works as expected."""

    # Create the paths to the data.
    path_free = f"{fep_output.name}/fep_output/{engine}/free"
    path_vac = f"{fep_output.name}/fep_output/{engine}/vacuum"

    # Perform the analysis.
    free_nrg_free, _ = BSS.FreeEnergy.Relative.analyse(
        path_free, estimator=estimator, preprocess=False
    )
    free_nrg_vac, _ = BSS.FreeEnergy.Relative.analyse(
        path_vac, estimator=estimator, preprocess=False
    )

    # Compute the free-energy difference.
    delta_g = BSS.FreeEnergy.Relative.difference(free_nrg_free, free_nrg_vac)

    # Make sure the result is close to the expected value.
    assert math.isclose(
        delta_g[0].value(), expected_results[engine][estimator], rel_tol=1e-4
    )
