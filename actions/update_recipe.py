import sys
import os
import subprocess

script = os.path.abspath(sys.argv[0])

# we want to import the 'get_requirements' package from this directory
sys.path.insert(0, os.path.dirname(script))

from parse_requirements import parse_requirements

# go up one directories to get the source directory
# (this script is in BioSimSpace/actions/)
srcdir = os.path.dirname(os.path.dirname(script))

condadir = os.path.join(srcdir, "recipes", "biosimspace")

print(f"conda recipe in {condadir}")

# Store the name of the recipe and template YAML files.
recipe = os.path.join(condadir, "meta.yaml")

# If the BSS_TEMPLATE environment variable is set, use that as the template.
if "BSS_TEMPLATE" in os.environ:
    template = os.environ["BSS_TEMPLATE"]
else:
    template = os.path.join(condadir, "template.yaml")

# Now parse all of the requirements
run_reqs = parse_requirements(os.path.join(srcdir, "requirements.txt"))
print(run_reqs)
build_reqs = parse_requirements(os.path.join(srcdir, "requirements_build.txt"))
print(build_reqs)


def run_cmd(cmd):
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    return str(p.stdout.read().decode("utf-8")).lstrip().rstrip()


gitdir = os.path.join(srcdir, ".git")

# Get the BSS remote.
remote = run_cmd(
    f"git --git-dir={gitdir} --work-tree={srcdir} config --get remote.origin.url"
)
remote += ".git"
print(remote)

# Get the BSS branch.
branch = run_cmd(
    f"git --git-dir={gitdir} --work-tree={srcdir} rev-parse --abbrev-ref HEAD"
)
print(branch)

lines = open(template, "r").readlines()


def dep_lines(deps):
    lines = []

    for dep in deps:
        lines.append(f"    - {dep}\n")

    return "".join(lines)


run_reqs = dep_lines(run_reqs)

if len(build_reqs) > 0:
    build_reqs = f"  build:\n{dep_lines(build_reqs)}"
else:
    build_reqs = ""

with open(recipe, "w") as FILE:
    for line in lines:
        if line.find("BSS_BUILD_REQUIREMENTS") != -1:
            line = build_reqs
        elif line.find("BSS_RUN_REQUIREMENTS") != -1:
            line = run_reqs
        else:
            line = line.replace("BSS_REMOTE", remote)
            line = line.replace("BSS_BRANCH", branch)

        FILE.write(line)
