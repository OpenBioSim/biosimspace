import os
import sys

from setuptools import setup, find_packages

import versioneer

# A list of authors and their email addresses.
authors = (
    "Lester Hedges <lester.hedges@gmail.com, "
    "Christopher Woods <chryswoods@gmail.com>, "
    "Antonia Mey <antonia.mey@gmail.com"
)

# Run the setup.
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
