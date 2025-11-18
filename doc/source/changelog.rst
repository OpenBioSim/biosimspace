Changelog
=========

OpenBioSim
----------

:mod:`BioSimSpace` is supported by `OpenBioSim <https://openbiosim.org>`__, a community interest
company supporting open-source development of fostering academic/industrial collaboration
within the biomolecular simulation community. Our software is hosted via the `OpenBioSim`
`GitHub <https://github.com/OpenBioSim/biosimspace>`__ organisation.

`2025.3.0 <https://github.com/openbiosim/biosimspace/compare/2025.2.0...2025.3.0>`_ - Nov 10 2025
-------------------------------------------------------------------------------------------------

* Fixed bug that caused trajectory frame reconstruction to be skipped for perturbable systems (`#460 <https://github.com/OpenBioSim/biosimspace/pull/460>`__).
* Fixed typo in RDKit ``rdMolDraw2D`` attribute name (`#463 <https://github.com/OpenBioSim/biosimspace/pull/463>`__).
* Fix incorrect context manager name in node runner (`#463 <https://github.com/OpenBioSim/biosimspace/pull/463>`__).
* Added functionality for quickly getting and setting the coordinates array of a :class:`System <BioSimSpace._SireWrappers.System>` (`#465 <https://github.com/OpenBioSim/biosimspace/pull/465>`__).
* Reduce depdency import overheads by switching from module level imports to function and method level (`#466 <https://github.com/OpenBioSim/biosimspace/pull/466>`__).
* Fully switch over to using the new Sire Python API naming convention, allowing BioSimSpace to be used within Sire (`#466 <https://github.com/OpenBioSim/biosimspace/pull/466>`__).
* Fixed logic used for setting the GPU device index for :class:`Process.OpenMM <BioSimSpace.Process.OpenMM>` (`#468 <https://github.com/OpenBioSim/biosimspace/pull/468>`__).

`2025.2.0 <https://github.com/openbiosim/biosimspace/compare/2025.1.0...2025.2.0>`_ - Oct 08 2025
-------------------------------------------------------------------------------------------------

* Add support for ``SOMD2`` FEP analysis using data frames with different numbers of samples (`#415 <https://github.com/OpenBioSim/biosimspace/pull/415>`__).
* Add support for the ABCG2 charge method for GAFF parameterisation (`#421 <https://github.com/OpenBioSim/biosimspace/pull/421>`__).
* Fixed f-string formatting error in FEP analysis exception message (`#423 <https://github.com/OpenBioSim/biosimspace/pull/423>`__).
* Fixed FEP energy trajectory slicing when intitial sample time is non-zero (`#424 <https://github.com/OpenBioSim/biosimspace/pull/424>`__).
* Add workaround for incompatibility between ``ParmEd`` and ``NumPy`` 2.3 (`#428 <https://github.com/OpenBioSim/biosimspace/pull/428>`__).
* Remove ``AnteChamber`` check when using ``OpenFF-NAGL`` (`#428 <https://github.com/OpenBioSim/biosimspace/pull/428>`__).
* Add missing OPC water solvation function docstring (`#431 <https://github.com/OpenBioSim/biosimspace/pull/431>`__).
* Use ``NumPy`` native ``array.tolist()`` rather than ``list(array)`` and add missing ATM protocol ``__str__`` methods (`#435 <https://github.com/OpenBioSim/biosimspace/pull/435>`__).
* Add missing OPC water solvation function docstring (`#431 <https://github.com/OpenBioSim/biosimspace/pull/431>`__).
* Allow user to skip writing of velocities to output files (`#437 <https://github.com/OpenBioSim/biosimspace/pull/437>`__).
* Add support for renumbering molecules when copying a system (`#440 <https://github.com/OpenBioSim/biosimspace/pull/440>`__).
* Added a molecular editing tutorial (`#447 <https://github.com/OpenBioSim/biosimspace/pull/447>`__).
* Fixed issue with a non-singleton being used in a set within a mapping pruning function (`#454 <https://github.com/OpenBioSim/biosimspace/pull/454>`__).
* Fix handling of trajectory frame reconstruction when using a perturbable reference system (`#456 <https://github.com/OpenBioSim/biosimspace/pull/456>`__).

`2025.1.0 <https://github.com/openbiosim/biosimspace/compare/2024.4.1...2025.1.0>`_ - Jul 01 2025
-------------------------------------------------------------------------------------------------

* Improved robustness of formal charge inference when reading molecules from PDB or SDF files (`#393 <https://github.com/OpenBioSim/biosimspace/pull/393>`__).
* Make sure the system extracted from AMBER trajectory frames during free-energy perturbation simulations is in the original, unsquashed format (`#403 <https://github.com/OpenBioSim/biosimspace/pull/403>`__).
* Add support for the ``ff19SB`` force field and OPC water (`#406 <https://github.com/OpenBioSim/biosimspace/pull/406>`__).
* Allow creation of ``SOMD`` perturbation files without modification to ghost atom bonded terms (`#407 <https://github.com/OpenBioSim/biosimspace/pull/407>`__).
* Support analysis of ``SOMD2`` energy trajectories with time varying lambda sampling (`#408 <https://github.com/OpenBioSim/biosimspace/pull/408>`__).

`2024.4.1 <https://github.com/openbiosim/biosimspace/compare/2024.4.0...2024.4.1>`_ - Feb 14 2025
-------------------------------------------------------------------------------------------------

* Add ``force_stereo_inference`` keyword argument to :func:`BioSimSpace.Convert.toRDKit <BioSimSpace.Convert.toRDKit>` to allow control over stereochemistry inference. This allows a user to match the old ``sire`` behaviour, where stereochemistry was always infered, which would override any information already present in a molecule (`#388 <https://github.com/OpenBioSim/biosimspace/pull/388>`__).
* Port optimised ``SOMD`` configuration settings to the Exscientia sandpit.

`2024.4.0 <https://github.com/openbiosim/biosimspace/compare/2024.3.0...2024.4.0>`_ - Feb 04 2025
-------------------------------------------------------------------------------------------------

* Improved robustness of restraint search and fixed torsion force constant (`@fjclark <https://github.com/fjclark>`_) (`#360 <https://github.com/OpenBioSim/biosimspace/pull/360>`__).
* Add support for ``acdoctor`` option for ``Antechamber`` to help with parameterisation of molecules with atypical valence (`#361 <https://github.com/OpenBioSim/biosimspace/pull/361>`__).
* Added missing AMBER atom mapping pruning options to the :mod:`BioSimSpace.Align` module (`#366 <https://github.com/OpenBioSim/biosimspace/pull/366>`__).
* Improved error messages when no lambda windows are found during alchemical free-energy analysis (`@fjclark <https://github.com/fjclark>`_) (`#368 <https://github.com/OpenBioSim/biosimspace/pull/368>`__).
* Decouple the selection used for aligment and RMSD for the :class:`RMSD <BioSimSpace.Metadynamics.CollectiveVariable.RMSD>` collective variable (`#379 <https://github.com/OpenBioSim/biosimspace/pull/379>`__).
* Added missing ``._to_default_unit()`` method to the :class:`GeneralUnit <BioSimSpace.Types._GeneralUnit>` class (`#385 <https://github.com/OpenBioSim/biosimspace/pull/385>`__).
* Clarified the meaning of the ``explicit_dummies`` keyword argument in the :class:`BioSimSpace.Process.Amber` class (`#385 <https://github.com/OpenBioSim/biosimspace/pull/385>`__).
* Added a new :class:`BioSimSpace.Types.Mass` class to handle atomic masses (`#385 <https://github.com/OpenBioSim/biosimspace/pull/385>`__).

`2024.3.0 <https://github.com/openbiosim/biosimspace/compare/2024.2.0...2024.3.0>`_ - Oct 24 2024
-------------------------------------------------------------------------------------------------

* Fixed link in protein free-energy perturbation tutorial (`#310 <https://github.com/OpenBioSim/biosimspace/pull/310>`__).
* Fixed alignment of monatomic molecules (`#313 <https://github.com/OpenBioSim/biosimspace/pull/313>`__ and (`#333 <https://github.com/OpenBioSim/biosimspace/pull/333>`__).
* Expose missing ``extra_args`` keyword argument for the :class:`Process.Somd <BioSimSpace.Process.Somd>` class (`#319 <https://github.com/OpenBioSim/biosimspace/pull/319>`__).
* Add support for the Alchemical Transfer Method (ATM) (`@mb2055 <https://github.com/mb2055>`_) (`#327 <https://github.com/OpenBioSim/biosimspace/pull/327>`__).
* Fixed :meth:`system.updateMolecules <BioSimSpace._SireWrappers.System.updateMolecules>` method when updating multiple molecules at once (`#336 <https://github.com/OpenBioSim/biosimspace/pull/336>`__).
* Added a new :meth:`system.removeBox <BioSimSpace._SireWrappers.System.removeBox>` method to remove the box from a system (`#338 <https://github.com/OpenBioSim/biosimspace/pull/338>`__).
* Fixed bug when using position restraints with OpenMM for perturbable systems (`#341 <https://github.com/OpenBioSim/biosimspace/pull/341>`__).
* Include both end states in the GROMACS restraint file when using position restraints with perturbable systems (`#342 <https://github.com/OpenBioSim/biosimspace/pull/342>`__).
* Set end state connectivity properties when a perturbation involves a change in connectivity (`#345 <https://github.com/OpenBioSim/biosimspace/pull/345>`__).
* Use direct Sire-to-RDKit conversion when parameterising with OpenFF to avoid issues with intermediate file formats (`#349 <https://github.com/OpenBioSim/biosimspace/pull/349>`__).
* Preserve SMILES based molecular properties during parameterisation (`#352 <https://github.com/OpenBioSim/biosimspace/pull/352>`__).
* Rotate perturbable and non-perturbable molecules independently when rotating box vectors (`#355 <https://github.com/OpenBioSim/biosimspace/pull/355>`__).

`2024.2.0 <https://github.com/openbiosim/biosimspace/compare/2024.1.0...2024.2.0>`_ - Jul 09 2024
-------------------------------------------------------------------------------------------------

* Fixed incorect use of ``self`` in :func:`Trajectory.getFrame <BioSimSpace.Trajectory.getFrame>` function (`#281 <https://github.com/OpenBioSim/biosimspace/pull/281>`__).
* Use SDF as an intermediate for ``antechamber`` if the original molecule was loaded from that format (`#287 <https://github.com/OpenBioSim/biosimspace/pull/287>`__).
* Detect dummy atoms by checking ``element`` *and* ``ambertype`` properties when creating ``SOMD`` pert files (`#289 <https://github.com/OpenBioSim/biosimspace/pull/289>`__).
* Add missing ``match_water`` kwarg to ``prepareFEP`` node (`#292 <https://github.com/OpenBioSim/biosimspace/pull/292>`__).
* Add protein free-energy perturbation functionality (`@akalpokas <https://github.com/akalpokas`_) (`#294 <https://github.com/OpenBioSim/biosimspace/pull/294>`__).
* Ensure that the LJ sigma parameter for perturbed atoms is non-zero (`#295 <https://github.com/OpenBioSim/biosimspace/pull/295>`__).
* Fixed return type docstrings for functions in the :mod:`BioSimSpace.Parameters` module (`#298 <https://github.com/OpenBioSim/biosimspace/pull/298>`__).
* Don't use ``sire.legacy.Base.wrap`` with the ``file_format`` property to avoid (incorrect) auto string to unit conversion of ``mol2`` to moles squared (`#300 <https://github.com/OpenBioSim/biosimspace/pull/300>`__).
* Expose ``SOMD`` torsion modification kwargs (`#302 <https://github.com/OpenBioSim/biosimspace/pull/302>`__).

`2024.1.0 <https://github.com/openbiosim/biosimspace/compare/2023.5.1...2024.1.0>`_ - Apr 15 2024
-------------------------------------------------------------------------------------------------

* Switch to using Langevin integrator for GROMACS free energy simulations (`#264 <https://github.com/OpenBioSim/biosimspace/pull/264>`__).
* Add support for clearing and disabling the IO file cache (`#266 <https://github.com/OpenBioSim/biosimspace/pull/266>`__).
* Add support for using ``openff-nagl`` to generate partial charges (`#267 <https://github.com/OpenBioSim/biosimspace/pull/267>`__).
* Fixed non-reproducible search for backbone restraint atom indices (`#270 <https://github.com/OpenBioSim/biosimspace/pull/270>`__).
* Add support for AMBER as an alchemical free-energy simulation engine (`#272 <https://github.com/OpenBioSim/biosimspace/pull/272>`__).
* Switch to using ``os.path.join`` to generate directory file names (`#276 <https://github.com/OpenBioSim/biosimspace/pull/276>`__).
* Make sure the ``fileformat`` property is preserved when creating single molecule systems (`#276 <https://github.com/OpenBioSim/biosimspace/pull/276>`__).
* Add a ``getRestraintType`` method to the base protocol that returns ``None`` (`#276 <https://github.com/OpenBioSim/biosimspace/pull/276>`__).
* Add missing ``themostat_time_constant`` kwarg to the :class:`Metadynamics <BioSimSpace.Protocol.Metdynamics>` protocol (`#276 <https://github.com/OpenBioSim/biosimspace/pull/276>`__).

`2023.5.1 <https://github.com/openbiosim/biosimspace/compare/2023.5.0...2023.5.1>`_ - Mar 20 2024
-------------------------------------------------------------------------------------------------

* Fixed path to user links file in the :func:`generateNetwork <BioSimSpace.Align.generateNetwork>` function (`#233 <https://github.com/OpenBioSim/biosimspace/pull/233>`__).
* Fixed redirection of stderr (`#233 <https://github.com/OpenBioSim/biosimspace/pull/233>`__).
* Switched to using ``AtomCoordMatcher`` to map parameterised molecules back to their original topology. This resolves issues where atoms moved between residues following parameterisation (`#235 <https://github.com/OpenBioSim/biosimspace/pull/235>`__).
* Make the GROMACS ``_generate_binary_run_file`` function static so that it can be used when initialising free energy simulations in setup-only mode (`#237 <https://github.com/OpenBioSim/biosimspace/pull/237>`__).
* Improve error handling and message when attempting to extract an all dummy atom selection (`#251 <https://github.com/OpenBioSim/biosimspace/pull/251>`__).
* Don't set SOMD specific end-state properties when decoupling a molecule (`#253 <https://github.com/OpenBioSim/biosimspace/pull/253>`__).
* Only convert to a end-state system when not running a free energy protocol with GROMACS so that hybrid topology isn't lost when using position restraints (`#257 <https://github.com/OpenBioSim/biosimspace/pull/257>`__).
* Exclude standard free ions from the AMBER position restraint mask (`#260 <https://github.com/OpenBioSim/biosimspace/pull/260>`__).
* Update the ``BioSimSpace.Types._GeneralUnit.__pow__`` operator to support fractional exponents (`#260 <https://github.com/OpenBioSim/biosimspace/pull/260>`__).

`2023.5.0 <https://github.com/openbiosim/biosimspace/compare/2023.4.1...2023.5.0>`_ - Dec 16 2023
-------------------------------------------------------------------------------------------------

* Add support for detecting nucleic acid backbones (`@fjclark <https://github.com/fjclark>`_) (`#189 <https://github.com/OpenBioSim/biosimspace/pull/189>`__).
* Added SOMD and GROMACS support for multiple distance restraints for ABFE calculations (`#178 <https://github.com/OpenBioSim/biosimspace/pull/178>`__).

`2023.4.1 <https://github.com/openbiosim/biosimspace/compare/2023.4.0...2023.4.1>`_ - Dec 14 2023
-------------------------------------------------------------------------------------------------

* Make sure ``match_water`` keyword argument is passed to specialised solvation functions (`#190 <https://github.com/OpenBioSim/biosimspace/pull/190>`__).
* Check perturbable molecules for velocities when combining molecules (`#192 <https://github.com/OpenBioSim/biosimspace/pull/192>`__).
* Make sure velocities are double counted when searching for velocity properties when combining molecules (`#197 <https://github.com/OpenBioSim/biosimspace/pull/197>`__).
* Remove redundant ``BioSimSpace.Types.Type.__ne__`` operator (`#201 <https://github.com/OpenBioSim/biosimspace/pull/201>`__).
* Minor internal updates due to Sire API fixes (`#203 <https://github.com/OpenBioSim/biosimspace/pull/203>`__).
* Fixed bug in the Boresch restraint search code (`@fjclark <https://github.com/fjclark>`_) (`#204 <https://github.com/OpenBioSim/biosimspace/pull/204>`__).
* Fixed ``renumber`` option in :meth:`extract <BioSimSpace._SireWrappers.Molecule.extract>` method (`#210 <https://github.com/OpenBioSim/biosimspace/pull/210>`__).
* Add workaround for fixing reconstruction of intrascale matrix in :func:`readPerturbableSystem <BioSimSpace.IO.readPerturbableSystem>` function (`#210 <https://github.com/OpenBioSim/biosimspace/pull/210>`__).
* Remove incorrect ``try_import`` statement in metadynamics driver script and make sure that global parameters in OpenMM script are unique (`#217 <https://github.com/OpenBioSim/biosimspace/pull/217>`__).
* Ensure the existing trajectory backend is used when getting the number of trajectory frames from a running process (`#219 <https://github.com/OpenBioSim/biosimspace/pull/219>`__).
* Fixed setting of ``igb`` config parameter	for PMEMD simulations (`@annamherz <https://github.com/annamherz>`_) (`#220 <https://github.com/OpenBioSim/biosimspace/pull/220>`__).
* Make sure AMBER restraint mask matches all hydrogen atoms (`#222 <https://github.com/OpenBioSim/biosimspace/pull/222>`__).
* Ensure all searches for disulphide bonds are convert to a ``SelectorBond`` object (`#224 <https://github.com/OpenBioSim/biosimspace/pull/224>`__).
* Fix injection of custom commands into ``LEaP`` script (`#226 <https://github.com/OpenBioSim/biosimspace/pull/226>`__).

`2023.4.0 <https://github.com/openbiosim/biosimspace/compare/2023.3.1...2023.4.0>`_ - Oct 13 2023
-------------------------------------------------------------------------------------------------

* Add support for computing trajectory RMSDs using Sire backend (`#152 <https://github.com/OpenBioSim/biosimspace/pull/152>`__).
* Add support for setting up systems containing crystal waters (`#154 <https://github.com/OpenBioSim/biosimspace/pull/154>`__).
* Add unified free-energy perturbation analysis using ``alchemlyb`` (`@annamherz <https://github.com/annamherz>`_) (`#155 <https://github.com/OpenBioSim/biosimspace/pull/155>`__).
* Fix handling of connectivity changes during molecular perturbations (`#157 <https://github.com/OpenBioSim/biosimspace/pull/157>`__).
* Fix issues related to new shared properties in Sire (`#160 <https://github.com/OpenBioSim/biosimspace/pull/160>`__).
* Fix issues in SOMD perturbation files for absolute binding free-energy simulations (`@fjclark <https://github.com/fjclark>`_) (`#164 <https://github.com/OpenBioSim/biosimspace/pull/164>`__).
* Don't generate velocities when performing a continuation with GROMACS (`#169 <https://github.com/OpenBioSim/biosimspace/pull/169>`__).
* Decouple custom parameters and additional commands in ``LEaP`` input (`#170 <https://github.com/OpenBioSim/biosimspace/pull/170>`__).
* Check for periodic space when updating box information from restart file or trajectory (`#173 <https://github.com/OpenBioSim/biosimspace/pull/173>`__).
* Add functionality to allow manual rotation and reduction of triclinic boxes, rather than performing automatically on read (`#175 <https://github.com/OpenBioSim/biosimspace/pull/175>`__).
* Allow unit-based protocol options to be passed as strings (`#179 <https://github.com/OpenBioSim/biosimspace/pull/179>`__).
* Fix assignment of ``gpu`` configuration option for SOMD (`#181 <https://github.com/OpenBioSim/biosimspace/pull/181>`__).

`2023.3.1 <https://github.com/openbiosim/biosimspace/compare/2023.3.0...2023.3.1>`_ - Aug 14 2023
-------------------------------------------------------------------------------------------------

* Check for non-periodic cartesian space when setting up vacuum simulations with all engines (`#125 <https://github.com/OpenBioSim/biosimspace/pull/125>`__).
* Fixed several issues caught by BioSimSpace tutorials suite (`#128 <https://github.com/OpenBioSim/biosimspace/pull/128>`__).
* Fixed import of incorrect ``alchemlyb`` extract function for GROMACS (`#132 <https://github.com/OpenBioSim/biosimspace/pull/132>`__).
* Handle issues with using certain triclinic box vectors with OpenMM by performing a pre lattice reduction using the internal OpenMM functionality (`#135 <https://github.com/OpenBioSim/biosimspace/pull/135>`__).
* Add support for OpenMM in example equilibration node (`@mb2055 <https://github.com/mb2055>`_) (`#138 <https://github.com/OpenBioSim/biosimspace/pull/138>`__).
* Fix use of ``totalSteps`` when using the OpenMM ``StateDataReporter`` (`#146 <https://github.com/OpenBioSim/biosimspace/pull/146>`__).
* Make sure ``alchemlyb`` is imported using ``try_import`` to avoid errors on platforms where it isn't available (`#151 <https://github.com/OpenBioSim/biosimspace/pull/151>`__).

`2023.3.0 <https://github.com/openbiosim/biosimspace/compare/2023.2.2...2023.3.0>`_ - Jun 30 2023
-------------------------------------------------------------------------------------------------

* Reinstate :data:`BioSimSpace.Stream <BioSimSpace.Stream>` sub-package (`#36 <https://github.com/OpenBioSim/biosimspace/pull/36>`__).
* Fixed ``setup.py`` file to work correctly on Windows (`#72 <https://github.com/OpenBioSim/biosimspace/pull/72>`__).
* Fixed bug with missing working directory when using ``rmsd_flex_align`` scoring function (`#75 <https://github.com/OpenBioSim/biosimspace/pull/75>`__).
* Use ``parmed`` to create ``openmm`` system to avoid issue parsing triclinic spaces with ``AmberPrmTopFile`` (`#77 <https://github.com/OpenBioSim/biosimspace/pull/77>`__).
* Fix parsing of AMBER free-energy perturbation standard output (`#79 <https://github.com/OpenBioSim/biosimspace/pull/79>`__).
* Fix bug in :data:`GeneralUnit <BioSimSpace.Types._GeneralUnit>` constructor (`#83 <https://github.com/OpenBioSim/biosimspace/pull/83>`__).
* Check molecule numbers in system when caching files to avoid issue when the UID and number of molecules are the same, but the actual molecules are different, e.g. after being edited (`#89 <https://github.com/OpenBioSim/biosimspace/pull/89>`__).
* Fix order of imports in ``prepareFEP`` node (`#90 <https://github.com/OpenBioSim/biosimspace/pull/90>`__).
* Recenter molecules following vacuum simulation with GROMACS to avoid precision overflow with molecular coordinates on write (`#95 <https://github.com/OpenBioSim/biosimspace/pull/95>`__).
* Fix expected angles used in unit test following updates to triclinic box code in Sire (`#99 <https://github.com/OpenBioSim/biosimspace/pull/99>`__).
* Add absolute binding free-energy support for SOMD (`@fjclark <https://github.com/fjclark>`_) (`#104 <https://github.com/OpenBioSim/biosimspace/pull/104>`__).
* Avoid streaming issues when reading binary AMBER restart files for a single frame (`#105 <https://github.com/OpenBioSim/biosimspace/pull/105>`__).
* Improve overlap matrix plotting functionality (`@fjclark <https://github.com/fjclark>`_) (`#107 <https://github.com/OpenBioSim/biosimspace/pull/107>`__).
* Handle updates to Sire parser format naming (`#108 <https://github.com/OpenBioSim/biosimspace/pull/108>`__).
* Wrap new Sire units grammar to improve parsing of units from strings (`#109 <https://github.com/OpenBioSim/biosimspace/pull/109>`__).
* Expose ``make_whole`` option in Sire to allow un-wrapping of molecular coordinates on read (`#110 <https://github.com/OpenBioSim/biosimspace/pull/110>`__).
* Make sure to call ``.value()`` on objects that now have units (`#110 <https://github.com/OpenBioSim/biosimspace/pull/110>`__).
* Handle missing values in AMBER standard output records (`#111 <https://github.com/OpenBioSim/biosimspace/pull/111>`__).
* Fix bug in ``plumed`` version requirement check (`#113 <https://github.com/OpenBioSim/biosimspace/pull/113>`__).
* Reinstate temperature control for all GROMACS simulation protocols (`#115 <https://github.com/OpenBioSim/biosimspace/pull/115>`__).
* Fix pre-processing selector in test section of ``conda`` recipe (`#117 <https://github.com/OpenBioSim/biosimspace/pull/117>`__).
* Fixed bug in SOMD free-energy perturbation analysis (`@fjclark <https://github.com/fjclark>`_) (`#119 <https://github.com/OpenBioSim/biosimspace/pull/119>`__).
* Catch exception when vacuum system has a cartesian space (`#120 <https://github.com/OpenBioSim/biosimspace/pull/120>`__).
* Add support for Sire as a trajectory backend (`#121 <https://github.com/OpenBioSim/biosimspace/pull/121>`__).

`2023.2.2 <https://github.com/openbiosim/biosimspace/compare/2023.2.1...2023.2.2>`_ - May 15 2023
-------------------------------------------------------------------------------------------------

* Rename tests directory to ``tests`` for compliance with ``pytest`` standard (`#51 <https://github.com/OpenBioSim/biosimspace/pull/51>`__).
* Fixed parsing of AMBER standard output records (`#56 <https://github.com/OpenBioSim/biosimspace/pull/56>`__).
* Re-add pre-minimisation stage to SOMD FEP configuration (`#59 <https://github.com/OpenBioSim/biosimspace/pull/59>`__).
* Fixed reference to ``plumed.dat`` file in AMBER configuration input for steered molecular dynamics (`#64 <https://github.com/OpenBioSim/biosimspace/pull/64>`__).
* Fixed :meth:`getDensity <BioSimSpace.Process.Amber.getDensity>` method (`#64 <https://github.com/OpenBioSim/biosimspace/pull/64>`__).

`2023.2.1 <https://github.com/openbiosim/biosimspace/compare/2023.2.0...2023.2.1>`_ - Apr 27 2023
-------------------------------------------------------------------------------------------------

* Update GitHub CI for our new release process (`#34 <https://github.com/OpenBioSim/biosimspace/pull/34>`__).
* Fixed :func:`readMolecules <BioSimSpace.IO.readMolecules>` so that can handle a tuple of input files again (`#38 <https://github.com/OpenBioSim/biosimspace/pull/38>`__).
* Fixed protocol mixin inheritance (`#41 <https://github.com/OpenBioSim/biosimspace/pull/41>`__).
* Update documentation for new development and release process (`#43 <https://github.com/OpenBioSim/biosimspace/pull/43>`__).
* Fixed SOMD inverse friction coefficient configuration parameter (`#49 <https://github.com/OpenBioSim/biosimspace/pull/49>`__).
* Fixes to the hydration free energy tutorial (`#49 <https://github.com/OpenBioSim/biosimspace/pull/49>`__).
* Fixed bug in SOMD test runner that caused it to return prior to assertions (`#49 <https://github.com/OpenBioSim/biosimspace/pull/49>`__).
* Expose ``extra_options`` and ``extra_lines`` parameters in :class:`BioSimSpace.FreeEnergy.Relative <BioSimSpace.FreeEnergy.Relative>` (`#49 <https://github.com/OpenBioSim/biosimspace/pull/49>`__).

`2023.2.0 <https://github.com/openbiosim/biosimspace/compare/2023.1.2...2023.2.0>`_ - Mar 30 2023
-------------------------------------------------------------------------------------------------

* Make sure that system properties are preserved when creating a new Sire system.
* Fixed an issue with the OpenMM minimisation protocol that meant that the number of steps was ignored (`#12 <https://github.com/OpenBioSim/biosimspace/pull/12>`__).
* Use native Sire PDB downloading functionality to remove ``pypdb`` dependency.
* Fixed an issue with SMILES characters in molecule names causing issues for ``gmx grompp`` (`#14 <https://github.com/OpenBioSim/biosimspace/pull/14>`__).
* Increase default SOMD cut-off since it uses reaction field (`#15 <https://github.com/OpenBioSim/biosimspace/pull/15>`__).
* No longer downcast molecules to single residues and atoms when searching (`#19 <https://github.com/OpenBioSim/biosimspace/pull/19>`__).
* Remove velocities when combining molecules if the property isn't present for all molecules (`#21 <https://github.com/OpenBioSim/biosimspace/pull/21>`__).
* Set default-valued properties when merging molecules to avoid issues with zero values when units are stripped (`#24 <https://github.com/OpenBioSim/biosimspace/pull/24>`__).
* Remove ``watchdog`` to avoid non-deterministic parsing of AMBER output (`#27 <https://github.com/OpenBioSim/biosimspace/pull/27>`__).
* Improved handling of disulphide bonds in multi-chain PDBs sharing the same residue numbers (`#28 <https://github.com/OpenBioSim/biosimspace/pull/28>`__).
* Allow keyword arguments to be passed through to ``lomap`` in :func:`generateNetwork <BioSimSpace.Align.generateNetwork>` (`#29 <https://github.com/OpenBioSim/biosimspace/pull/29>`__).
* Add mixin classes to allow position restraints to be used with a wider range of protocols (`@xiki-tempula <https://github.com/xiki-tempula>`_) and alchemical simulations for non-production protocols (`@msuruzhon <https://github.com/msuruzhon>`_). Switch to using ``gmx energy`` to parse GROMACS energy records (`@xiki-tempula <https://github.com/xiki-tempula>`_) (`#30 <https://github.com/OpenBioSim/biosimspace/pull/30>`__).
* Switch to using native RDKit conversion throughout to avoid conversion via an intermediate file format.
* Expose Sire to OpenMM conversion functionality in :mod:`BioSimSpace.Convert <BioSimSpace.Convert>`.
* Added Python 3.10 support and now build Python 3.10 packages. This is now the default version of Python for BioSimSpace, and the version we recommend for new workflows. Note that we will drop automatic building of Python 3.8 packages later this year (likely Q3 or Q4). This will be timed to co-incide with when we add Python 3.11 support, and when (we anticipate) conda-forge will drop Python 3.8. Our aim is to only build packages for a maximum of 3 Python versions at a time.

`2023.1.2 <https://github.com/openbiosim/biosimspace/compare/2023.1.1...2023.1.2>`_ - Feb 24 2023
-------------------------------------------------------------------------------------------------

* Refactor code to use a unified :class:`WorkDir <BioSimSpace._Utils.WorkDir>` class to simplify the creation of working directories (`#2 <https://github.com/OpenBioSim/biosimspace/pull/2>`__).
* Added :meth:`isSame <BioSimSpace._SireWrappers.System.isSame>` method to compare systems using a sub-set of system and molecular properties. This improves our file caching support, allowing a user to exclude properties when comparing cached systems prior to write, e.g. ignoring coordinates and velocities, if those are the only things that differ between the systems `(#3 <https://github.com/OpenBioSim/biosimspace/pull/3>`__).
* Added the initial version of :mod:`BioSimSpace.Convert <BioSimSpace.Convert>`, which provides support for converting between native `BioSimSpace`, `Sire <http://sire.openbiosim.org>`__, and `RDKit <https://www.rdkit.org>`__ objects (`#9 <https://github.com/OpenBioSim/biosimspace/pull/9>`__).
* Fixed several formatting issues with the website documentation.

`2023.1.1 <https://github.com/openbiosim/biosimspace/compare/2023.1.0...2023.1.1>`_ - Feb 07 2023
-------------------------------------------------------------------------------------------------

* Minor fixes to website documentation.
* Fixed issues with API documentation introduced by `pydocstringformatter <https://pypi.org/project/pydocstringformatter>`__.
* Fixed globbing of GROMACS trajectory files.

`2023.1.0 <https://github.com/openbiosim/biosimspace/compare/2022.3.0...2023.1.0>`_ - Feb 03 2023
-------------------------------------------------------------------------------------------------

* Wrapped the new `sire.load <https://sire.openbiosim.org/api/sire.html#sire.load>`__ function to allow loading of URLs.
* Add basic file caching support to avoid re-writing files for the same molecular system.
* Added :data:`BioSimSpace._Config` sub-package to simplify the generation of configuration files for molecular dynamics engines and improve flexiblity. (Adapted from code written by `@msuruzhon <https://github.com/msuruzhon>`_.)
* Deprecated ``BioSimSpace.IO.glob`` since globbing is now performed automatically.
* Autoformat entire codebase using `black <https://github.com/psf/black>`__.
* Fix issues following Sire 2023 API updates.
* Update documentation for new OpenBioSim website.

Michellab
---------

Prior to January 2023, :mod:`BioSimSpace` was hosted within the `michellab <https://github.com/michellab/BioSimSpace>`__
GitHub organisation. The following releases were made during that time.

`2022.3.0 <https://github.com/openbiosim/biosimspace/compare/2022.2.1...2022.3.0>`_ - Sep 28 2022 (Pre-release)
---------------------------------------------------------------------------------------------------------------

* Improved NAMD restraint implementation for consistency with other engines.
* Make sure we wait for ``trjconv`` to finish when calling as a sub-process.
* Added wrapper for ``Sire.Units.GeneralUnit``.
* Improved interoperability of ``BioSimSpace.Trajectory`` sub-package.
* Added ``BioSimSpace.Sandpit`` for experimental features from external collaborators.
* Added functionality to check for molecules in a ``BioSimSpace.System``.
* Added functionality to extract atoms and residues by absolute index.
* Allow continuation for GROMACS equilibration simulations. (`@kexul <https://github.com/kexul>`_)
* Update BioSimSpace to work with the new Sire 2023.0.0 Python API.

`2022.2.1 <https://github.com/openbiosim/biosimspace/compare/2022.2.0...2022.2.1>`_ - Mar 30 2022
-------------------------------------------------------------------------------------------------

* Fix performance issues when ensuring unique molecule numbering when adding molecules to ``BioSimSpace._SireWrappers.System`` and ``BioSimSpace._SireWrappers.Molecules`` objects.
* Fix extraction of box vector magnitudes for triclinic boxes.

`2022.2.0 <https://github.com/openbiosim/biosimspace/compare/2022.1.0...2022.2.0>`_ - Mar 24 2022
-------------------------------------------------------------------------------------------------

* Use fast C++ wrappers for updating coordinates and velocities during SOMD simulations.
* Fix import issues caused by change in module layout for conda-forge OpenMM package.
* Don't check for structural ions when parameterising with GAFF/GAFF2.
* Fix errors in funnel correction calculation.
* Switch to using conda-forge lomap2 package, removing need to vendor lomap code.
* Use py3Dmol to visualise maximum common substructure mappings.
* Rename ``.magnitude()`` method on ``BioSimSpace.Type`` objects to ``.value()`` to avoid confusion.
* Handle trjconv frame extraction failures within ``BioSimSpace.Process.Gromacs.getSystem()``.
* Catch and handle possible GSL error during singular valued decomposition routine used for molecular alignment.

`2022.1.0 <https://github.com/openbiosim/biosimspace/compare/2020.1.0...2022.1.0>`_ - Jan 26 2022
-------------------------------------------------------------------------------------------------

* Added basic support for cleaning PDB files with `pdb4amber <https://github.com/Amber-MD/pdb4amber>`_ prior to read.
* Added basic support for exporting BioSimSpace Nodes as Common Workflow Language wrappers.
* Added support for parameterising molecules using OpenForceField.
* Added support for using SMILES strings for input to parameterisation functions.
* Added support for funnel metadynamics simulations (`@dlukauskis <https://github.com/dlukauskis>`_).
* Added support for steered molecular dynamics simulations (`@AdeleLip <https://github.com/AdeleLip>`_).
* Added support for generating perturbation networks using LOMAP (`@JenkeScheen <https://github.com/JenkeScheen>`_).
* Fixed bug affecting certain improper/dihedral terms in SOMD perturbation file writer.
* Numerous performance improvements, particularly involving the manipulation and
  combination of molecular systems.
* Native Python pickling support for wrapped Sire types (`@chryswoods <https://github.com/chryswoods>`_).
* Numerous free-energy perturbation pipeline fixes and improvements. Thanks to `@kexul <https://github.com/kexul>`_ and `@msuruzhon <https://github.com/msuruzhon>`_ for their help testing and debugging.
* Switch continuous integration to GitHub actions using conda-forge compliant build and upload to Anaconda cloud.

`2020.1.0 <https://github.com/openbiosim/biosimspace/compare/2019.3.0...2020.1.0>`_ - July 28 2020
--------------------------------------------------------------------------------------------------

* Added logo to website and update theme (`@ppxasjsm <https://github.com/ppxasjsm>`_).
* Make sure potential terms are sorted when writing to SOMD perturbation files (`@ptosco <https://github.com/ptosco>`_).
* Switch to using ipywidgets.FileUpload to eliminate non-conda dependencies.
* Added support for single-leg free energy simulations.
* Created a KCOMBU mirror to avoid network issues during install.
* Allow AMBER simulations when system wasn't loaded from file.
* Handle GROMACS simulations with non-periodic boxes.
* Run vacuum simulations on a single thread when using GROMACS to avoid domain decomposition.
* Make sure BioSimSpace is always built against the latest version of Sire during conda build.

`2019.3.0 <https://github.com/openbiosim/biosimspace/compare/2019.2.0...2019.3.0>`_ - Nov 22 2019
-------------------------------------------------------------------------------------------------

* Make FKCOMBU download during conda build resilient to server downtime.
* Added support for xtc trajectory files and custom protocols with GROMACS.
* Fixed numerous typos in Sphinx documentation.
* Added Journal of Open Source Software paper.

`2019.2.0 <https://github.com/openbiosim/biosimspace/compare/2019.1.0...2019.2.0>`_ - Sep 11 2019
-------------------------------------------------------------------------------------------------

* Switched to using `RDKit <https://www.rdkit.org/>`_ for maximum common substructure (MCS) mappings.
* Handle perturbable molecules for non free-energy protocols with SOMD and GROMACS.
* Added basic metadynamics functionality with support for distance and torsion collective variables.
* Added support for inferring formal charge of molecules.
* Numerous MCS mapping fixes and improvements. Thanks to `@maxkuhn <https://github.com/maxkuhn>`_, `@dlukauskis <https://github.com/dlukauskis>`_, and `@ptosco <https://github.com/ptosco>`_ for help testing and debugging.
* Added Dockerfile to build thirdparty packages required by the BioSimSpace notebook server.
* Exposed Sire search functionality.
* Added thin-wrappers for several additional Sire objects, e.g. Residue, Atom, and Molecules container.
* Performance improvements for searching, indexing, and extracting objects from molecular containers, e.g. System, Molecule.

`2019.1.0 <https://github.com/openbiosim/biosimspace/compare/2018.1.1...2019.1.0>`_ - May 02 2019
-------------------------------------------------------------------------------------------------

* Added support for parameterising proteins and ligands.
* Added support for solvating molecular systems.
* Molecular dynamics drivers updated to support SOMD and GROMACS.
* Support free energy perturbation simulations with SOMD and GROMACS.
* Added Azure Pipeline to automatically build, test, document, and deploy BioSimSpace.
* Created automatic Conda package pipeline.

`2018.1.1 <https://github.com/openbiosim/biosimspace/compare/2018.1.0...2018.1.1>`_ - May 02 2018
-------------------------------------------------------------------------------------------------

* Fixed conda NetCDF issue on macOS. Yay for managing `python environments <https://xkcd.com/1987>`_\ !
* Install conda `ambertools <https://anaconda.org/AmberMD/ambertools>`_ during `setup <python/setup.py>`_.
* Search for bundled version of ``sander`` when running `AMBER <http://ambermd.org>`_ simulation processes.
* Pass executable found by ``BioSimSpace.MD`` to ``BioSimSpace.Process`` constructor.
* Fixed error in RMSD calculation within ``BioSimSpace.Trajectory`` class.
* Improved example scripts and notebooks.

2018.1.0 - May 01 2018
----------------------

* Initial public release of BioSimSpace.
