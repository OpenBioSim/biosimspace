.. _ref_install:

============
Installation
============

1. No-installation - Run in a web browser
=========================================

We run a completely free `JupyterHub <https://try.openbiosim.org>`__ on
which we have BioSimSpace installed.

This is at `try.openbiosim.org <https://try.openbiosim.org>`__.
You only need a `GitHub account <https://github.com>`__, which is
used to log into the server.

Simply go to `try.openbiosim.org <https://try.openbiosim.org>`__ in your
web browser, and log in using your `GitHub account <https://github.com>`__.
This will start a Jupyter Lab instance. In here you can start a terminal,
and then use :mod:`BioSimSpace` directly via ``ipython``. Or you can start a Jupyter
notebook and use :mod:`BioSimSpace` there.

To import :mod:`BioSimSpace`, at a Python prompt type

>>> import BioSimSpace as BSS

If this imports without errors, then everything is working.
We encourage you to take a look at :doc:`the tutorials <tutorials/index>`
to learn how to use :mod:`BioSimSpace` or the
:doc:`quickstart guide <quickstart/index>` if you want an overview.

.. note::

   This free JupyterHub server is limited. You only have up to 2 GB of
   memory and at most 1 processor core. Disk storage is temporary,
   and any data will be lost when you log out. Because it only
   supports a limited number of concurrent users, inactive sessions will be
   automatically stopped and logged out after 20 minutes. Please only
   use this service to explore and learn :mod:`BioSimSpace`.
   Do not use it for production work.

2. Easy installation - Run in a conda environment
=================================================

The easiest way to install :mod:`BioSimSpace` is in a new
`conda environment <https://anaconda.org>`__.

You can use any conda environment or installation. We recommend using
`Miniforge <https://github.com/conda-forge/miniforge>`__,
as this is pre-configured to use `conda-forge <https://conda-forge.org>`__.

.. _Install_Miniforge:
Either... Install a new copy of ``Miniforge``
----------------------------------------------

To install a new copy of
`Miniforge <https://github.com/conda-forge/miniforge>`__,
first download a ``Miniforge`` from
`this page <https://github.com/conda-forge/miniforge>`__ that
matches your operating system and processor.

Install ``Miniforge`` following the
`instructions here <https://github.com/conda-forge/miniforge#install>`__.

Once installed, you should be able to run the ``conda`` command to
install other packages (e.g. ``conda -h`` will print out help on
how to use the ``conda`` command).

Or... Use an existing anaconda/miniconda install
------------------------------------------------

If you want to use an existing anaconda or miniconda installation,
then first open a terminal with that distribution activated.
For example, open a terminal via anaconda navigator, or
open a terminal and run
``source /path/to/conda/bin/activate``, where ``/path/to/conda`` is
the full path to your anaconda or miniconda installation.

You should now be able to run the ``conda`` command to install other
packages (e.g. ``conda -h`` will print out help on how to use the
``conda`` command).

And then... Install BioSimSpace into a new environment
------------------------------------------------------

We recommend that :mod:`BioSimSpace` is installed into a new (clean) environment.
This minimises the risk of failures caused by incompatible dependencies.

BioSimSpace is currently packaged for Python 3.8 and Python 3.9. We will start
by creating a Python 3.9 environment that we will call ``openbiosim``.

.. code-block:: bash

   $ conda create -n openbiosim "python<3.10"

.. note::

   We use ``python<3.10`` as this will install the most recent 3.9
   release of python.

We can now install :mod:`BioSimSpace` into that environment by typing

.. code-block:: bash

   $ conda install -n openbiosim -c openbiosim biosimspace

.. note::

   The option ``-n openbiosim`` tells ``conda`` to install :mod:`BioSimSpace`
   into the ``openbiosim`` environment. The option ``-c openbiosim``
   tells ``conda`` to install :mod:`BioSimSpace` from the ``openbiosim``
   conda channel.

If you want the latest development release, then install by typing

.. code-block:: bash

   $ conda install -n openbiosim -c "openbiosim/label/dev" biosimspace

To install the latest development version you can use:

.. code-block:: bash

    conda create -n openbiosim-dev -c conda-forge -c openbiosim/label/dev biosimspace
    conda activate openbiosim-dev

To run :mod:`BioSimSpace`, you must now activate the ``openbiosim`` environment.
You can do this by typing

.. code-block:: bash

   $ conda activate openbiosim

You can now start a Python session (e.g. running ``python``, or
``ipython`` or ``jupyter lab`` if you installed those). At the
Python prompt you can import :mod:`BioSimSpace` by typing

>>> import BioSimSpace as BSS

If this imports without errors, then everything is working.
We encourage you to take a look at :doc:`the tutorials <tutorials/index>`
to learn how to use :mod:`BioSimSpace` or the
:doc:`quickstart guide <quickstart/index>` if you want an overview.

3. Also easy installation - Run in a container
==============================================

Another route to install :mod:`BioSimSpace` is to download and run our
pre-built containers. These can be run via
`docker <https://www.docker.com>`__ (on Linux, MacOS and Windows)
or via `podman <https://podman.io>`__ (on Linux) on Intel (X86-64)
or ARM64 processors.

To run via `docker <https://www.docker.com>`__, simply type;

.. code-block:: bash

   $ docker run -p 8888:8888 -it openbiosim/biosimspace:latest

or, via `podman <https://podman.io>`__, type;

.. code-block:: bash

   $ podman run -p 8888:8888 -it openbiosim/biosimspace:latest

This will download the container from
`hub.docker.com <https://hub.docker.com/r/openbiosim/biosimspace>`__ and
will start a command prompt in that container.

You can now type ``python``, ``ipython`` or ``jupyter lab``
to start a python, ipython or jupyter lab session.

.. note::

   The option ``-p 8888:8888`` tells docker/podman to redirect
   port ``8888`` on your computer to port ``8888`` in the
   container. This will let you open a browser and navigate to
   the URL printed by ``jupyter lab`` if you are using jupyter.
   You can drop this option if you don't want to use
   ``jupyter lab``.

.. note::

   You can map directories from your computer into the container
   by using the ``-v`` option. For example,
   ``-v $HOME/input:/home/openbiosim/input`` would map your
   ``input`` folder in your home directory to the ``input`` folder
   in the home directory of the container. This will let :mod:`BioSimSpace`
   read and write files on your computer.

You can now start a Python session (e.g. running ``python``, or
``ipython`` or ``jupyter lab`` if you installed those). At the
Python prompt you can import :mod:`biosimspace` by typing

>>> import BioSimSpace as BSS

If this imports without errors, then everything is working.
We encourage you to take a look at :doc:`the tutorials <tutorials/index>`
to learn how to use :mod:`BioSimSpace` or the
:doc:`quickstart guide <quickstart/index>` if you want an overview.


4. Harder installation - Compile from source
============================================

The following documents a full installation of BioSimSpace from source. Before
starting, you'll need a working `Git <https://git-scm.com>`__ installation.

BioSimSpace is built on top of the `Sire <https://github.com/openbiosim/sire>`__
molecular simulation framework. To download and install Sire, follow the
instructions `here <https://sire.openbiosim.org/install.html>`__, making
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

5. Common issues
================

* If you experience problems with `Matplotlib <https://matplotlib.org>`__ when
  importing BioSimSpace on macOS, e.g.

.. code-block:: bash

   RuntimeError**: Python is not installed as a framework.

simply add the following to ``~/.matplotlib/matplotlibrc``

.. code-block:: bash

   backend: TkAgg

Note that plotting functionality will be disabled if you are using
BioSimSpace on a remote server without X forwarding.

* If you experience problems with `Jupyter <https://jupyter.org>`__ permissions,
  try removing ``$HOME/.jupyter`` or ``$HOME/.local/share/jupyter``

6. External dependencies
========================

Several additional packages are required for full access to all of BioSimSpace's
functionality. Please download and install these packages according to their
recommended installation instructions.

* `Amber / AmberTools <http://ambermd.org>`__ -- *Dynamics / Parameterisation*
* `Gromacs <http://www.gromacs.org>`__ -- *Dynamics / Parameterisation / Solvation*
* `Namd <http://www.ks.uiuc.edu/Research/namd>`__ -- *Dynamics*

(Note that BioSimSpace is built to be compatible with the ``ambertools`` and
``gromacs`` packages from conda-forge, but they are not included as hard
requirement. This means that BioSimSpace can be used in conda environments
with and without them. We've taken this approach to enable the use of stripped
down environments, and to avoid clashes with external versions of the packages,
which may be better optimised for specific hardware and usage requirements.)

For `Amber / AmberTools <http://ambermd.org>`__, we also recommend adding
``${AMBERHOME}/bin`` to your ``PATH`` to ensure that its binaries are
visible to third-party libraries, such as
`openff-toolkit <https://github.com/openforcefield/openff-toolkit>`__.

Please visit our :ref:`compatibility <ref_compatibility>` page to see which
versions of the external dependencies BioSimSpace has currently been tested
against.
