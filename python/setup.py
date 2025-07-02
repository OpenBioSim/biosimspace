import os
import sys
import platform

from setuptools import setup, find_packages

import versioneer

# A list of authors and their email addresses.

# Allow setting on the command line because setting
# environment variables on Windows is painful!
if "BSS_SKIP_DEPENDENCIES=1" in sys.argv:
    os.environ["BSS_SKIP_DEPENDENCIES"] = "1"
    sys.argv.remove("BSS_SKIP_DEPENDENCIES=1")

if "BSS_CONDA_INSTALL=1" in sys.argv:
    os.environ["BSS_CONDA_INSTALL"] = "1"
    sys.argv.remove("BSS_CONDA_INSTALL=1")

authors = (
    "Lester Hedges <lester.hedges@gmail.com, "
    "Christopher Woods <chryswoods@gmail.com>, "
    "Antonia Mey <antonia.mey@gmail.com"
)

_installed_list = None


# Function to check if a conda dependency has been installed
def is_installed(dep: str, conda: str):
    global _installed_list

    if _installed_list is None:
        p = subprocess.Popen([conda, "list", dep], stdout=subprocess.PIPE)
        lines = str(p.stdout.read())
        _installed_list = lines

    return _installed_list.find(dep) != -1


# Function to clear the cache of installed packages
def clear_installed_list():
    global _installed_list
    _installed_list = None


# Run the setup.
try:
    setup(
        name="BioSimSpace",
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        description="BioSimSpace: Making biomolecular simulation a breeze.",
        author=authors,
        url="https://github.com/openbiosim/biosimspace",
        license="GPLv3",
        packages=find_packages(),
        include_package_data=True,
        zip_safe=False,
    )

# Post setup configuration.
finally:
    import os

    if "install" in sys.argv and not (
        os.getenv("BSS_CONDA_INSTALL") or os.getenv("BSS_SKIP_DEPENDENCIES")
    ):
        import shlex
        import subprocess

        # Install Python dependencies and enable Jupyter widget extensions.
        print("\nSetting up python environment...")

        # Open files for stdout/stderr.
        stdout = sys.stdout
        stderr = sys.stderr

        # Create a list of the conda dependencies.
        conda_deps = [
            "alchemlyb",  # known not available on aarch64
            "configargparse",
            "ipywidgets",
            "kcombu_bss",
            "lomap2",
            "mdtraj",  # known not available on aarch64
            "mdanalysis",  # known not available on aarch64
            "networkx",
            "nglview",
            "openff-interchange-base",
            "openff-toolkit-base",
            "parmed",
            "pyarrow",
            "py3dmol",
            "pydot",
            "pygtail",
            "pyyaml",
            "rdkit",
            "sire",
        ]

        # Don't try to install things that are already installed...
        to_install_deps = []

        print("Checking for dependencies that are already installed...")

        from shutil import which

        conda_exe = which("mamba")
        real_conda_exe = which("conda")

        if conda_exe is None or not os.path.exists(conda_exe):
            conda_exe = real_conda_exe

        if conda_exe is None or not os.path.exists(conda_exe):
            raise IOError("Unable to install as cannot find conda or mamba!")

        for dep in conda_deps:
            if not is_installed(dep, conda=conda_exe):
                to_install_deps.append(dep)
                print(f"Need to install {dep}")
            else:
                print("Already installed %s" % dep)

        clear_installed_list()

        conda_deps = to_install_deps

        # Need to not use posix rules on windows with shlex.split, or path separator is escaped
        posix = sys.platform != "win32"

        print("Adding openbiosim channel")
        command = (
            "%s config --system --prepend channels openbiosim/label/dev"
            % real_conda_exe
        )
        print(command)

        print("Adding conda-forge channel")
        command = "%s config --system --prepend channels conda-forge" % real_conda_exe
        print(command)
        try:
            subprocess.run(
                shlex.split(command, posix=posix),
                shell=False,
                stdout=stdout,
                stderr=stderr,
            )
        except Exception as e:
            print(f"Something went wrong ({e}). Continuing regardless...")

        print("Disabling conda auto update")
        command = "%s config --system --set auto_update_conda false" % real_conda_exe
        print(command)
        try:
            subprocess.run(
                shlex.split(command, posix=posix),
                shell=False,
                stdout=stdout,
                stderr=stderr,
            )
        except Exception as e:
            print(f"Something went wrong ({e}). Continuing regardless...")

        print("Installing conda dependencies: %s" % ", ".join(conda_deps))
        command = "%s install -y -q %s" % (conda_exe, " ".join(conda_deps))
        print(command)

        all_installed_ok = True

        try:
            subprocess.run(
                shlex.split(command, posix=posix),
                shell=False,
                stdout=stdout,
                stderr=stderr,
                check=True,
            )
        except Exception:
            all_installed_ok = False

        if not all_installed_ok:
            print("There were errors installing some of the dependencies.")
            print(
                "We will now try to install them one-by-one. This may take some time..."
            )

            failures = []

            for dep in conda_deps:
                if not is_installed(dep, conda=conda_exe):
                    print("Trying again to install '%s'" % dep)
                    command = "%s install -y -q %s" % (conda_exe, dep)

                    try:
                        subprocess.run(
                            shlex.split(command, posix=posix),
                            shell=False,
                            stdout=stdout,
                            stderr=stderr,
                            check=True,
                        )
                    except Exception:
                        failures.append(dep)

            if len(failures) == 0:
                print("All dependencies installed successfully!")
            else:
                print(
                    "\n** Failed to install these dependencies: %s"
                    % ", ".join(failures)
                )
                print(
                    "** BioSimSpace will still install and run, but some functionality may not be available.\n"
                )
        else:
            print("All dependencies install successfully first time!")

        print("Activating notebook extension: nglview")

        if sys.platform == "win32":
            bin_dir = os.path.join(bin_dir, "Scripts")

        print("Cleaning conda environment")
        command = "%s clean --all --yes --quiet" % conda_exe
        subprocess.run(
            shlex.split(command, posix=posix),
            shell=False,
            stdout=stdout,
            stderr=stderr,
        )

        print("\nDone!")

        print(
            "\nIf you have problems with Jupyter permissions, try removing '$HOME/.jupyter' or '$HOME/.local/share/jupyter'"
        )

        print("\nFor optional package support...")
        print("AMBER:   http://ambermd.org")
        print("GROMACS: http://www.gromacs.org")
        print("NAMD:    http://www.ks.uiuc.edu/Research/namd")
