.. _ref-DevProcess:

===================
Development process
===================

:mod:`BioSimSpace` uses a ``main`` and ``devel`` development process,
using feature branches for all code development.

* ``main`` - this always contains the latest official release.
* ``devel`` - this always contains the latest development release, which will become the next official release.

Code should be developed on a fork or in a feature branch called ``feature_{feature}``.
When your feature is ready, please submit a pull request against ``devel``. This
will trigger our continuos integration (CI) process, which will build
:mod:`BioSimSpace` on a range of different platforms. All merge conflicts must
be fixed in the branch and all tests must pass before the pull request can be
merged into ``devel``.( NOTE THAT ONLY AUTHORISED DEVELOPERS CAN ACCEPT THE
PULL REQUEST.) Authorised developers will review the pull request as quickly
as they can.  They will be greatly helped if the feature is accompanied with
tests, examples and/or tutorial instructions.

.. note::

  We encourage that new functionality in a pull request is documented, e.g. by adding to
  the tutorials or writing a detailed description for the website.

Assuming the CI completes successfully, then one of the release team will
conduct a code review, with the code being merged into ``devel`` if it is
approved.

Bug fixes or issue fixes are developed on fix branches, called ``fix_{number}`` (again in
either the main repository or forks). If no `issue thread <https://github.com/OpenBioSim/biosimspace/issues>`__
exists for the bug that you are fixing, please open one to notify other users of the problem.
This also allows us to correctly cross-reference and categorise pull requests. Once a fix
is ready, please submit a pull request against ``devel``. Assuming that the CI passes and
the reviewers are happy, this will then be merged. The same approach should also be used
for updates to documentation.

Once a fix has been approved and merged into ``devel`` one of the core developers will
backport it to the previous major release, i.e. by applying the fix to the ``main`` branch.
To do so, create a new fix branch based on ``main``. If the files that were updated on
``devel`` as part of the fix don't contain any other differences to those in ``main``,
i.e. they haven't been updated during the merging of a new feature, then you can
directly apply those to main. This can be done via:

.. code-block:: bash

   git checkout devel file

where ``file`` is the path to the file that was changed. If multiple files in a directory
were changed, then you can use:

.. code-block:: bash

   git checkout devel directory

Note that the above approach will automatically stage the modified files. If you would
prefer to be able to check the diff and make additional edits, then use ``git restore``
instead, e.g.:

.. code-block:: bash

   git restore --source devel file

If the fix applied on ``devel`` was applied on top of a new feature, then you might be
able to use ``git cherry-pick`` to apply the specific commits that relate only to the
fix, e.g.

.. code-block:: bash

   git cherry-pick commitSha

where ``commitSha`` is a commit reference. Use the ``--edit`` option if you want
to change the commit message, and use ``--no-commit`` if you just want to apply
the edits to the files on ``main``, but not automatically make a new commit, e.g.
if you want to check the results. (For more details, see `here <https://git-scm.com/docs/git-cherry-pick>`__
or `here <https://www.atlassian.com/git/tutorials/cherry-pick>`__.)

Once the fix has been successfully applied, please raise a pull request against
``main`` for one of the development team to review. In cases where the updated
files are identical to ``devel``, it may not be necessary to ask for a review or
run CI, since this would have been done when the fix was applied to ``devel``.
Once approved, the fix be merged into ``main``.

Within a release cycle we will periodically create point releases on ``main``,
e.g. ``2023.1.1``. (The frequency will depend on the urgency of the fixes.)

.. note::

  The group of developers authorised to have access to the
  `main BioSimSpace repository <https://github.com/openbiosim/biosimspace>`__
  and to accept pull requests is not fixed,
  and will evolve over time. If you wish to join this group then
  please complete the tutorial and then demostrate your commitment
  by submitting good issues and pull requests from
  a personal fork of the repository. Please get in touch if you find
  this difficult, or follow
  `this workshop <https://chryswoods.com/introducing_git>`__
  and `this workshop <https://chryswoods.com/git_collaboration>`__ if you need
  to learn how to use Git, GitHub, feature branching, merging, pull
  requests etc.

=================
Creating releases
=================

We use a calendar-based version numbering system, based on regular release cadence
of several releases per year. Our aim is to make a major release every quarter
(so four releases per year). Each release will be numbered sequentially, e.g.
2023.1.0 is the first release of 2023, 2023.2.0 is the second release. Our aim
is that new functionality only appears in these “major” releases.

Major releases
--------------

There are a number of stages to go through to create a major release:

1. Make sure that all changes (features and fixes) that are required for the
   release have been merged into the ``devel`` branch, and the GitHub Action
   has run fully, building :mod:`BioSimSpace` on all supported platforms,
   running all the unit tests correctly, and building and uploading the conda
   packages to the ``dev`` channel on `anaconda.org <https://anaconda.org/openbiosim/biosimspace/files>`__.

2. Create feature branches from ``devel`` to synchronise updates from external
   sandpits, e.g. `here <https://github.com/Exscientia/BioSimSpace>`__. (This
   may also be submitted as an external PR from the industrial partner.)
   Once ready, submit a pull request so that these updates can be tested and
   reviewed. (It might be the case the an external partner needs to run a larger
   set of internal tests against the changes.) When ready, this can be merged
   into ``devel`` by an authorised developer.

3. The next task is to create a pull request to merge ``devel`` into ``main``.
   At this point there might be conflicts to resolve due to previous backporting
   of fixes from ``devel`` into ``main``. Since ``devel`` should now be the
   *source of truth*, an approach could be to perform the merge in the reverse
   fashion, taking *everything* from ``devel``, after which the ``release`` branch
   has the exact same tree as ``devel``, with all of the history and tags from
   ``main``.

.. code-block:: bash

   git checkout -b release devel
   git merge -s ours main

4. On the ``release`` branch, you will now need to update ``doc/source/changelog.rst``
   with the changes since the last release. Start the entry with a title that (will)
   link to the changes on GitHub,
   e.g.  `https://github.com/openbiosim/biosimspace/compare/2022.2.1...2023.1.0 <https://github.com/openbiosim/biosimspace/compare/2022.2.1...2023.1.0>`__
   links to the changes between the ``2022.2.1`` and ``2023.1.0`` releases.
   Follow a similar format for changes as already exist in this file. Try to
   link to pull requests and tutorials that describe new functionality if available.

5. Typically `Sire <https://sire.openbiosim.org/index.html>`__ will be release at the
   same time, so you will need to update the version number in the
   `requirements.txt <https://github.com/OpenBioSim/biosimspace/blob/devel/requirements.txt>`__
   file. Note that there are different versions of Sire for ``main`` and
   ``devel``, so you will need to ensure that *only* the ``main`` version is
   uncommented and increment it to the latest version of Sire available on ``main``.

6. Now create a pull request for the ``release`` branch against ``main``. Once the
   CI passes and it has been approved by an authorised developer it can be merged.

7. An authorised developer will now create a tag for this release on ``main``, e.g.
   using ``git tag -as {VERSION} -m "{VERSION} release"``, e.g.
   ``git tag -as 2023.1.0 -m "2023.1.0 release"`` would be the tag for the
   ``2023.1.0`` release. The tag can then be pushed to GitHub with
   ``git push origin tag``, where ``tag`` is the new tag that you've created,
   e.g. ``2023.1.0``.

8. Now we are finally ready to build the packages. This can be done by an
   authorised developer by triggering the *workflow dispatch* event for the
   `Release Main <https://github.com/OpenBioSim/biosimspace/actions/workflows/main.yaml>`__
   workflow. At this point you will need to choose ``yes`` to upload packages
   to `anaconda.org <https://anaconda.org/openbiosim/biosimspace>`__.

9. (*) GitHub Actions don't currently build the ARM64/aarch64 packages.
   These have to be built and uploaded manually. On a MacOS/M1 or
   Linux/aarch64 computer you should create build environments for
   the Python versions that :mod:`BioSimSpace` should support. Activate
   this environment, and then checkout the ``main`` branch, run
   ``python actions/update_recipe.py`` and then run ``conda-build``
   via the command ``conda mambabuild -c conda-forge -c openbiosim/label/main recipes/biosimspace``.
   This will result in a conda package in the ``conda-bld`` directory
   in the root directory of your conda environment. You then need
   to upload these packages to `anaconda.org <https://anaconda.org/openbiosim/biosimspace>`__,
   e.g. via the command ``anaconda --token {PASSWORD} upload --user openbiosim --label main {/path/to/biosimspace-packages}``
   (modified as appropriate to include the anaconda password and the
   path to the built conda package).

10. On GitHub, you can now create a release by using the
    `Draft a New Release <https://github.com/OpenBioSim/biosimspace/releases/new>`__
    link. Choose the version number for your release from the tag you created
    earlier. The text should be simple, e.g. titled ``BioSimSpace {VERSION}``,
    with the body ``This is the {VERSION} release of BioSimSpace.``, along with
    a link to the changelog for the release on the website (if a major release).

11. Next you should build the docker images for this release.
    Do this by following the instructions in the
    `containers repository <https://github.com/OpenBioSim/containers/blob/main/biosimspace/README.md>`__.
    You should make sure to run the extra command listed there
    to tag the container with the version number you used earlier.

12. Finally(!) you can now update the website. To do this, follow the
    instructions in the `website repository <https://github.com/OpenBioSim/biosimspace_website/blob/main/README.md>`__.

13. Bonus! Follow the instructions in the
    `containers repository <https://github.com/OpenBioSim/containers>`__
    to build the notebook container image and instruct
    `try.openbiosim.org <https://try.openbiosim.org>`__ to update
    and use that image.

14. Super-bonus! If you have time, please write a short news item piece
    that can be added to the `openbiosim website <https://openbiosim.org>`__
    to announce this new release.

(We are in the process of automating many of the above steps, so hope
that this process will become much easier in the future.)

Point releases
--------------

Between major releases it might be necessary to create *point* releases to fix
bugs. Fixes will have been applied following the :ref:`ref-DevProcess` and a
release can be made by following steps 7 through 11 above. There is no need
to create a changelog entry or update the website for a point release. (Unless
fixes apply to the documentation itself.)

Development releases
--------------------

Following a major release, the next commit to ``devel`` will be given a new
development tag based against the *next* major release number. For example,
following the ``2023.1.0`` release, ``devel`` will be tagged with ``2023.2.0.dev``.
Following this, future merges into ``devel`` will trigger builds of development
packages that will be pushed to our `anaconda.org <https://anaconda.org/openbiosim/biosimspace>`__
channel using the ``dev`` label. This allows users to pull and test features
that will be available in the next major release on ``main``.

Periodically an authorised user will update the development documentation on our
`website repository <https://github.com/OpenBioSim/biosimspace_website/blob/main/README.md>`__
by running the workflow dispatch event against the ``devel`` branch.
