`BioSimSpace <https://biosimspace.openbiosim.org>`__
====================================================

.. image:: https://github.com/openbiosim/biosimspace/actions/workflows/devel.yaml/badge.svg
   :target: https://github.com/openbiosim/biosimspace/actions?query=workflow%3ARelease-Devel
   :alt: Build status

.. image:: https://anaconda.org/openbiosim/biosimspace/badges/downloads.svg
   :target: https://anaconda.org/openbiosim/biosimspace
   :alt: Conda Downloads

.. image:: https://img.shields.io/badge/License-GPL%20v3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0.html
   :alt: License

.. image:: https://joss.theoj.org/papers/4ba84ad443693b5dded90e35bf5f8225/status.svg
   :target: https://joss.theoj.org/papers/4ba84ad443693b5dded90e35bf5f8225
   :alt: Paper

About
-----

`BioSimSpace <https://biosimspace.openbiosim.org>`__ is an interoperable Python framework
for biomolecular simulation. With it you can:

* Write robust and portable biomolecular workflow components that work on
  different hardware, with different software packages, and that can be
  run in different ways, e.g. command-line, `Jupyter <https://jupyter.org>`__.
* Start, stop, and monitor molecular simulation processes within interactive Python environments.

Citation |DOI for Citing BioSimSpace|
=====================================

If you use BioSimSpace in any scientific software, please cite the following paper: ::

    @article{Hedges2019,
      doi = {10.21105/joss.01831},
      url = {https://doi.org/10.21105/joss.01831},
      year = {2019},
      publisher = {The Open Journal},
      volume = {4},
      number = {43},
      pages = {1831},
      author = {Lester Hedges and Antonia Mey and Charles Laughton and Francesco Gervasio and Adrian Mulholland and Christopher Woods and Julien Michel},
      title = {BioSimSpace: An interoperable Python framework for biomolecular simulation},
      journal = {Journal of Open Source Software}
    }

.. |DOI for Citing BioSimSpace| image:: https://joss.theoj.org/papers/4ba84ad443693b5dded90e35bf5f8225/status.svg
   :target: https://joss.theoj.org/papers/4ba84ad443693b5dded90e35bf5f8225

Documentation
-------------

Full documentation can be found `here <https://biosimspace.openbiosim.org>`__.

Installation
------------

Conda package
^^^^^^^^^^^^^

The easiest way to install BioSimSpace is using our `conda channel <https://anaconda.org/openbiosim/repo>`__.
BioSimSpace is built using dependencies from `conda-forge <https://conda-forge.org/>`__,
so please ensure that the channel takes strict priority. We recommend using
`Miniforge <https://github.com/conda-forge/miniforge>`__.

To create a new environment:

.. code-block:: bash

    conda create -n openbiosim -c conda-forge -c openbiosim biosimspace
    conda activate openbiosim

To install the latest development version you can use:

.. code-block:: bash

    conda create -n openbiosim-dev -c conda-forge -c openbiosim/label/dev biosimspace
    conda activate openbiosim-dev

When updating the development version it is generally advised to update `Sire <https://github.com/openbiosim/sire>`_
at the same time:

.. code-block:: bash

    conda update -c conda-forge -c openbiosim/label/dev biosimspace sire

Unless you add the required channels to your Conda configuration, then you'll
need to add them when updating, e.g., for the development package:

.. code-block:: bash

    conda update -c conda-forge -c openbiosim/label/dev biosimspace

Installing from source
^^^^^^^^^^^^^^^^^^^^^^

Alternatively, to install BioSimSpace from source:

(Before starting, you'll need a working `Git <https://git-scm.com>`__ installation.)

BioSimSpace is built on top of the `Sire <https://github.com/openbiosim/sire>`__
molecular simulation framework. To download and install Sire, follow the
instructions `here <https://github.com/openbiosim/sire#installation>`__, making
sure that BioSimSpace's dependencies are installed into the Sire conda
environment at the point at which Sire is installed.

Next you will need to download BioSimSpace and install it into your Sire
Conda environment.

.. code-block:: bash

   git clone https://github.com/openbiosim/biosimspace
   cd biosimspace/python
   pip install .

If you plan to develop and want an editable install, use:

.. code-block:: bash

   pip install -e .

If you want to skip the installation of BioSimSpace dependencies, e.g. if they
are already installed, then you can use:

.. code-block:: bash

   BSS_SKIP_DEPENDENCIES=1 pip install -e .

Once finished, you can test the installation by running:

.. code-block:: bash

   python

Then try importing the BioSimSpace package:

.. code-block:: python

   import BioSimSpace as BSS

If you don't want to install Sire from source, an alternative is to create a conda
environment containing only the dependencies of BioSimSpace, then install the
latest development code into that.

.. code-block:: bash

   conda create -n openbiosim-dev -c conda-forge -c openbiosim/label/dev biosimspace --only-deps
   conda activate openbiosim-dev
   git clone https://github.com/openbiosim/biosimspace
   cd biosimspace/python
   BSS_SKIP_DEPENDENCIES=1 pip install -e .

(You may also want to install optional dependencies, such as ``ambertools`` and
``gromacs`` into your environment.)

Developers
----------

Please follow the `developer's guide <https://biosimspace.openbiosim.org/contributing>`__.

Issues
------

Please report bugs and other issues using the GitHub `issue tracker <https://github.com/openbiosim/biosimspace/issues>`__.
When reporting issues please try to include a minimal code snippet that reproduces
the problem. Additional files can be also be uploaded as an archive, e.g. a zip
file. Please also report the branch on which you are experiencing the issue,
along with the BioSimSpace version number. This can be found by running:

.. code-block:: python

   import BioSimSpace as BSS
   print(BSS.__version__)
