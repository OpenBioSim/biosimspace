============
Coding Style
============

We really appreciate your help developing :mod:`BioSimSpace` and
welcome pull requests to the ``devel`` branch. To help us more
quickly review your pull request, and to keep a consistent coding
style throughout, we ask that you please follow the below coding
styles (and that you aren't offended if we modify your submission
so that it meets these styles).

Python
======

BioSimSpace uses the Python programming language. Our aim is to provide a simple
and robust API where unnecessary implementation details are hidden from the user.

"\ *Hold on a second, this code isn't very Pythonic!*\ "

Indeed it is not. BioSimSpace is built on top of an existing C++ framework and
the wrapped objects are designed to mimic the underlying C++ API. In addition,
this is partly a design choice, since BioSimSpace is intended primarily
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


(Note that `sire BioSimSpace repository <https://sire.openbiosim.org>`__, on top
of which BioSimSpace is built, has recently undergone a modernisation program
and now has a fully `PEP8-compliant <https://pep8.org>`__ API. We indent to
move BioSimSpace towards using Python compliant API too, while preserving backwards
compatibility.)

We use `ruff <https://docs.astral.sh/ruff>`__ to format and lint our Python
code. This is run automatically via `pre-commit <https://pre-commit.com>`__
hooks when you commit changes. See the :doc:`development` guide for setup
instructions.

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
