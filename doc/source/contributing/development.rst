=================
Developer's guide
=================

.. toctree::
   :maxdepth: 2

The source code for BioSimSpace is available on `GitHub <https://github.com/openbiosim/biosimspace>`__.

Python
======

BioSimSpace uses the Python programming language. Our aim is to provide a simple
and robust API where unnecessary implementation details are hidden from the user.

"\ *Hold on a second, this code isn't very Pythonic!*\ "

Indeed it is not, but this is a design choice. BioSimSpace is intended primarily
to be used by novices, who may be unfamiliar with Python, or programming in
general. We want to make it as easy as possible for these users to get up and
running with molecular simulation. BioSimSpace also needs to be robust and
portable, hence we need to use encapsulation to shield the user from unintended
consequences.

With this in mind, we use the following coding conventions:

Naming
------

We follow a C++ style naming convention.

* Packages: CamelCase
* Classes: CamelCase
* Methods: camelCase
* Functions: camelCase
* Variables: snake_case

For example, to instantiate a minimisation protocol from the ``Protocol`` package:

.. code-block:: python

   import BioSimSpace as BSS
   protocol = BSS.Protocol.Minimisation()

Modules
-------

BioSimSpace is a collection of packages, e.g. ``BioSimSpace.Gateway`` and
``BioSimSpace.Protocol``. Within each package is a set of modules that
implement the required functionality. Rather than directly exposing all of
the modules we choose to hide implementation details from the user. Instead
we use the package ``__init__.py`` to selectively import the required
classes and functions.


* Module files containing implementation details are prefixed with an underscore,
  i.e. ``_process.py``

* Where possible, external packages are hidden inside each module,
  e.g. ``import mdtraj as _mdtraj``

* Each module file contains an ``__all__`` variable that lists the specific items
  that should be imported.

* The package ``__init__.py`` can be used to safely expose the required
  functionality to the user with:

.. code-block:: python

   from module import *

This results in a clean API and documentation, with all extraneous information,
e.g. external modules, hidden from the user. This is important when working
interactively, since `IPython <https://ipython.org>`__ and `Jupyter <https://jupyter.org>`__
do not respect the ``__all__`` variable when auto-completing, meaning that the
user will see a full list of the available names when hitting tab. When
following the conventions above, the user will only be able to access the
exposed names. This greatly improves the clarity of the package, allowing
a new user to quickly determine the available functionality. Any user wishing
expose further implementation detail can, of course, type an underscore to
show the hidden names when searching.

Encapsulation
-------------

BioSimSpace aims to provide a means of writing robust and portable workflow
components (nodes). To this end, we choose to use an object oriented approach
where data is encapsulated, with getters used to retrieve data from an object.

To avoid unintended consequences, getters that return mutable data types, e.g.
lists and dictionaries, should return a copy of the data. This prevents the
user unintentionally modifying the private data contained in the object. Setters
should be used to explicitly modify member data.

For example:

.. code-block:: python

   # A class that holds a list of numbers.

   class MyClass():
       # A private class member variable containing a list of numbers.
       _list = [1, 2, 3, 4, 5]

       def getList(self):
           return self._list

   # Create an instance of the class.
   c = MyClass()
   n = c.getList()
   print(n)
   [1, 2, 3, 4, 5]

   # Update n.
   n.append(6)

   # The private member data has been modified!
   print(c.getList())
   [1, 2, 3, 4, 5, 6]

Instead use:

.. code-block:: python

   class MyClass():
       # A private class member variable containing a list of numbers.
       _list = [1, 2, 3, 4, 5]

       def getList(self):
           return self._list.copy()

   # Create an instance of the class.
   c = MyClass()
   n = c.getList()
   print(n)
   [1, 2, 3, 4, 5]

   # Update n.
   n.append(6)

   # The private member data is untouched.
   print(c.getList())
   [1, 2, 3, 4, 5]

Sire
----

BioSimSpace is built on top of the `Sire <https://github.com/openbiosim/sire>`__
and makes exstensive use of its API. Sire has recently been updated to a new,
Python compliant API. Within BioSimSpace, Sire is currently imported in a
*mixed* mode so that we can use *both* the new and legacy APIs while we
work to port things over. A consequence of this is that Sire imports
*must* follow those for BioSimSpace, e.g.:

.. code-block:: python

     from sire.mol import AtomIdx
     import BioSimSpace as BSS

will fail, while the following will work:

.. code-block:: python

     import BioSimSpace as BSS
     from sire.mol import AtomIdx


Workflow
========

Feature branches
----------------

First make sure that you are on the development branch of BioSimSpace:

.. code-block:: bash

   git checkout devel

Now create and switch to a feature branch. This should be prefixed with
*feature*, e.g.

.. code-block:: bash

   git checkout -b feature-process

While working on your feature branch you won't want to continually re-install
in order to make the changes active. To avoid this, install the package in
editable mode:

.. code-block:: bash

   pip install -e .

Pre-commit hooks
----------------

Pre-commit hooks are used to ensure consistent code formatting and linting.
To set up pre-commit in your development environment:

.. code-block:: bash

   pixi shell -e dev
   pre-commit install

This will run `ruff <https://docs.astral.sh/ruff/>`__ formatting and linting
checks automatically on each commit. To run the checks manually against all
files:

.. code-block:: bash

   pre-commit run --all-files

Testing
-------

When working on your feature it is important to write tests to ensure that it
does what is expected and doesn't break any existing functionality. Tests
should be placed inside the ``test`` directory, creating an appropriately named
sub-directory for any new packages.

The test suite is intended to be run using `pytest <https://docs.pytest.org/en/latest/contents.html>`__.
When run, ``pytest`` searches for tests in all directories and files below the current
directory, collects the tests together, then runs them. Pytest uses name matching
to locate the tests. Valid names start or end with *test*\ , e.g.:

.. code-block:: python

   # Files:
   test_file.py       file_test.py

   # Functions:
   def test_func():   def func_test():

We use the convention of ``test_*`` when naming files and functions.

Running tests
^^^^^^^^^^^^^

To run the full test suite, simply type:

.. code-block:: bash

   pytest tests

To run tests for a specific sub-module, e.g.:

.. code-block:: bash

   pytest tests/Process

To only run the unit tests in a particular file, e.g.:

.. code-block:: bash

   pytest tests/Process/test_namd.py

To run a specific unit tests in a particular file, e.g.:

.. code-block:: bash

   pytest tests/Process/test_namd.py::test_minimise

To get more detailed information about each test, run pytests using the
*verbose* flag, e.g.:

.. code-block:: bash

   pytest -v

More details regarding how to invoke ``pytest`` can be found `here <https://docs.pytest.org/en/latest/usage.html>`__.

Writing tests
^^^^^^^^^^^^^

Basics
""""""

Try to keep individual unit tests short and clear. Aim to test one thing, and
test it well. Where possible, try to minimise the use of ``assert`` statements
within a unit test. Since the test will return on the first failed assertion,
additional contextual information may be lost.

Floating point comparisons
""""""""""""""""""""""""""

Make use of the `approx <https://docs.pytest.org/en/latest/builtin.html#comparing-floating-point-numbers>`__
function from the ``pytest`` package for performing floating point comparisons, e.g:

.. code-block:: python

   from pytest import approx

   assert 0.1 + 0.2 == approx(0.3)

By default, the ``approx`` function compares the result using a relative tolerance
of 1e-6. This can be changed by passing a keyword argument to the function, e.g:

.. code-block:: python

   assert 2 + 3 == approx(7, rel=2)

Skipping tests
""""""""""""""

If you are using `test-driven development <https://en.wikipedia.org/wiki/Test-driven_development>`__
it might be desirable to write your tests before implementing the functionality,
i.e. you are asserting what the *output* of a function should be, not how it should
be *implemented*. In this case, you can make use of the ``pytest`` *skip* decorator
to flag that a unit test should be skipped, e.g.:

.. code-block:: python

   @pytest.mark.skip(reason="Not yet implemented.")
   def test_new_feature():
       # A unit test for an, as yet, unimplemented feature.
       ...

Parametrizing tests
"""""""""""""""""""

Often it is desirable to run a test for a range of different input parameters.
This can be achieved using the ``parametrize`` decorator, e.g.:

.. code-block:: python

   import pytest
   from operator import mul

   @pytest.mark.parametrize("x", [1, 2])
   @pytest.mark.parametrize("y", [3, 4])
   def test_mul(x, y):
       """ Test the mul function. """
       assert mul(x, y) == mul(y, x)

Here the function test_mul is parametrized with two parameters, ``x`` and ``y``.
By marking the test in this manner it will be executed using all possible
parameter pairs ``(x, y)``\ , i.e. ``(1, 3), (1, 4), (2, 3), (2, 4)``.

Alternatively:

.. code-block:: python

   import pytest
   from operator import sub
   @pytest.mark.parametrize("x, y, expected",
                           [(1, 2, -1),
                            (7, 3,  4),
                            (21, 58, -37)])
   def test_sub(x, y, expected):
       """ Test the sub function. """
       assert sub(x, y) == -sub(y, x) == expected

Here we are passing a list containing different parameter sets, with the names
of the parameters matched against the arguments of the test function.

Testing exceptions
""""""""""""""""""

Pytest provides a way of testing your code for known exceptions. For example,
suppose we had a function that raises an ``IndexError``\ :

.. code-block:: python

   def indexError():
       """ A function that raises an IndexError. """
       a = []
       a[3]

We could then write a test to validate that the error is thrown as expected:

.. code-block:: python

   def test_indexError():
       with pytest.raises(IndexError):
           indexError()

Custom attributes
"""""""""""""""""

It's possible to mark test functions with any attribute you like. For example:

.. code-block:: python

   @pytest.mark.slow
   def test_slow_function():
       """ A unit test that takes a really long time. """
       ...

Here we have marked the test function with the attribute ``slow`` in order to
indicate that it takes a while to run. From the command line it is possible
to run or skip tests with a particular mark.

.. code-block:: python

   pytest mypkg -m "slow"        # only run the slow tests
   pytest mypkg -m "not slow"    # skip the slow tests

The custom attribute can just be a label, as in this case, or could be your
own function decorator.

Documentation
-------------

BioSimSpace is fully documented using `NumPy <https://numpy.org>`__ style
docstrings. See `here <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`__
for details. The documentation is automatically built using
`Sphinx <http://sphinx-doc.org>`__ whenever a commit is pushed to devel, which
will then update this website.

To build the documentation locally you will first need to install some
additional packages.

.. code-block:: bash

   pip install sphinx sphinx_issues sphinx_rtd_theme

Then move to the ``doc`` directory and run:

.. code-block:: bash

   sphinx-build make html

When finished, point your browser to ``build/html/index.html``.

Committing
----------

If you create new tests, please make sure that they pass locally before
committing. When happy, commit your changes, e.g.

.. code-block:: bash

   git commit src/BioSimSpace/Feature/new_feature.py test/Feature/test_feature \
       -m "Implementation and test for new feature."

Remember that it is better to make small changes and commit frequently.

If your edits don't change the BioSimSpace source code, or documentation,
e.g. fixing typos, then please add ``ci skip`` to your commit message.
This will avoid unnecessarily triggering
`GitHub actions <https://github.com/openbiosim/biosimspace/actions?query=workflow%3ABuild>`__, e.g.
building a new BioSimSpace binary, updating the website, etc. To this end, we
have provided a git hook that will append ``ci skip`` if the commit only
modifies files in a blacklist that is specified in the file ``.ciignore``
(analogous to the ``.gitignore`` used to ignore untracked files). To enable
the hook, simply copy it into the ``.git/hooks`` directory:

.. code-block:: bash

    cp git_hooks/commit-msg .git/hooks

Any additional files or paths that shouldn't trigger a re-build can be added
to the ``.ciignore`` file.

Next, push your changes to the remote server, e.g.

.. code-block:: bash

   # Push to the feature branch on the main BioSimSpace repo, if you have access.
   git push origin feature

   # Push to the feature branch your own fork.
   git push fork feature

When the feature is complete, create a *pull request* on GitHub so that the
changes can be merged back into the development branch. For information, see
the documentation `here <https://help.github.com/articles/about-pull-requests>`__.
