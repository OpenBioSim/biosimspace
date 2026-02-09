"""Upload built packages to the openbiosim Anaconda Cloud channel."""

import glob
import os
import subprocess
import sys

script = os.path.abspath(sys.argv[0])

# Go up one directory to get the source directory.
srcdir = os.path.dirname(os.path.dirname(script))

print(f"BioSimSpace source is in {srcdir}\n")

# Get the anaconda token to authorise uploads.
if "ANACONDA_TOKEN" in os.environ:
    conda_token = os.environ["ANACONDA_TOKEN"]
else:
    conda_token = "TEST"

# Get the anaconda channel labels.
if "ANACONDA_LABEL" in os.environ:
    conda_label = os.environ["ANACONDA_LABEL"]
else:
    conda_label = "dev"

# Search for rattler-build output first.
packages = glob.glob(os.path.join("output", "**", "*.conda"), recursive=True)

# Fall back to conda-bld output.
if not packages:
    if "CONDA" in os.environ:
        conda = os.environ["CONDA"]
        conda_bld = os.path.join(conda, "envs", "bss_build", "conda-bld")
        packages = glob.glob(
            os.path.join(conda_bld, "**", "biosimspace-*.tar.bz2"), recursive=True
        )

if not packages:
    print("No BioSimSpace packages to upload?")
    sys.exit(-1)

print("Uploading packages:")
for pkg in packages:
    print(f"  * {pkg}")

packages_str = " ".join(packages)

# Upload the packages to the openbiosim channel on Anaconda Cloud.
cmd = f"anaconda --token {conda_token} upload --user openbiosim --label {conda_label} --force {packages_str}"

print(f"\nUpload command:\n\n{cmd}\n")

if conda_token == "TEST":
    print("Not uploading as the ANACONDA_TOKEN is not set!")
    sys.exit(-1)


def run_cmd(cmd):
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    return str(p.stdout.read().decode("utf-8")).lstrip().rstrip()


output = run_cmd(cmd)
print(output)
print("Package uploaded!")
