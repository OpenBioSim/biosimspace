# Script that collects as much as it can from a failed conda build so that it can
# be stored as a GitHub Actions artifact for download and further debugging

import glob
import os
import sys
import tarfile

if "BUILD_DIR" not in os.environ:
    try:
        os.environ["BUILD_DIR"] = sys.argv[1]
    except Exception:
        print("You need to supply BUILD_DIR")
        sys.exit(-1)

build_dir = os.environ["BUILD_DIR"]

# We want to bzip up the last 'sire-*' and all 'broken-*' directories in build_dir
work_dirs = glob.glob(os.path.join(build_dir, "biosimspace_*", "work", "build"))
broken_dirs = glob.glob(os.path.join(build_dir, "broke*"))

if len(work_dirs) > 0:
    work_dirs = [work_dirs[-1]]

zipdirs = work_dirs + broken_dirs

output_filename = os.path.join(build_dir, "failed.tar.bz2")

print(f"Zipping up {zipdirs} to {output_filename}")


def filter_function(tarinfo):
    filename = tarinfo.name
    # print(filename)
    if filename.find(".git") != -1:
        # print("excluded!")
        return None
    else:
        return tarinfo


with tarfile.open(output_filename, "w:bz2") as tar:
    for dir in zipdirs:
        tar.add(dir, arcname=os.path.basename(dir), filter=filter_function)

print("Complete :-)")
