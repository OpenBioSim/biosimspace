.. _ref-FreeEnergy:

BioSimSpace.FreeEnergy
======================

The *FreeEnergy* package contains tools to configure, run, and analyse
*relative* free energy simulations.

.. automodule:: BioSimSpace.FreeEnergy

.. toctree::
   :maxdepth: 1

As well as the :class:`protocol <BioSimSpace.Protocol.FreeEnergy>` used for production

Free-energy perturbation simulations require a
:class:`System <BioSimSpace._SireWrappers.System>` containing a *merged*
molecule that can be *perturbed* between two molecular end states by use
of an *alchemical potential*. To create merged molecules, please use the
:ref:`ref-Align` package.

Relative free-energy calculations require the simulation of two perturbations,
typically referred to as *legs*. A potential of mean force (PMF) is computed
for each leg, which can then be used to computed the relative free-energy
difference. For generality and flexibility, BioSimSpace decouples the two legs,
allowing the use of difference molecular simulation engines,
:class:`protocols <BioSimSpace.Protocol.FreeEnergy>`, and for legs to be
re-used in different calculations.

Simulations are typically used to compute solvation (currently hydration only)
or binding free-energies. In the examples that follow, ``merged`` refers to a
perturbable molecule created by merging two ligands, ``ligA`` and ``ligB``,
``merged_sol`` refers to the same perturbable molecule in solvent, and
``complex_sol`` is a solvated protein-ligand complex containing the same
perturbable molecule. We assume that each molecule/system has been
appropriately minimised and equlibrated.

To setup, run, and analyse a binding free-energy calculation:

.. code-block:: python

   import BioSimSpace as BSS

   ...

   # Create two a protocol for the two legs of a binding free-energy simulation.
   # Use more lambda windows for the "bound" leg.
   protocol_bound = BSS.Protocol.FreeEnergy(num_lam=20)
   protocol_free  = BSS.Protocol.FreeEnergy(num_lam=12)

   # Setup the perturbations for each leg, using the SOMD engine. This will
   # create all of the input files and simulation processes that are required.
   fep_bound = BSS.FreeEnergy.Relative(
       complex_sol,
       protocol_bound,
       engine="somd",
       work_dir="ligA_ligB/bound"
   )
   fep_free  = BSS.FreeEnergy.Relative(
       merged_sol,
       protocol_bound,
       engine="somd",
       work_dir="ligA_ligB/free"
   )

   # Run all simulations for each leg. Note that the lambda windows are run
   # sequentially, so this is a sub-optimal way of executing the simulation
   # if you have access to HPC resources.

   # Bound leg.
   fep_bound.run()
   fep_bound.wait()

   # Free leg.
   fep_free.run()
   fep_free.wait()

   # Analyse the simulation data from each leg, returning the PMF and overlap
   # matrix.
   pmf_bound, overlap_bound = fep_bound.analyse()
   pmf_free,  overlap_free  = fep_free.analyse()

   # Compute the relative free-energy difference.
   free_nrg_binding = BSS.FreeEnergy.Relative.difference(pmf_bound, pmf_free)

Similarly, for a solvation free-energy calculation:

.. code-block:: python

   # Here we are assuming that we are using the same ligands, so will re-use
   # the free leg from the previous example.

   # Setup the perturbation for the vacuum leg using a default protocol.
   fep_vacuum = BSS.FreeEnergy.Relative(
       merged.toSystem(),
       engine="somd",
       work_dir="ligA_ligB/vacuum"
   )

   # Run the simulations for the perturbation.
   fep_vacuum.run()
   fep_vacuum.wait()

   # Analyse the simulation data.
   pmf_vacuum, overlap_vacuum = fep_vacuum.analyse()

   # Compute the relative free-energy difference.
   free_nrg_solvation = BSS.FreeEnergy.Relative.difference(pmf_free, pmf_vacuum)

Since it is usually preferable to run simulations intensive simulation such as
these on external HPC resources, the ``BioSimSpace.FreeEnergy`` package also
provides support for only creating the input files that are needed by passing
the ``setup_only=True`` argument. This saves the overhead of creating
:class:`Process <BioSimSpace.Process>` objects. The input files can then be
copied to a remote server, with the indivual simulations curated in a job
submission script. (We don't yet provide support for configuring and writing
submission scripts for you.)

To just setup the vacuum leg input files:

.. code-block:: python

   # Setup the input for the vacuum leg. No processes are created so the .run()
   # method won't do anything.
   fep_vacuum = BSS.FreeEnergy.Relative(
       merged.toSystem(),
       engine="somd",
       work_dir="ligA_ligB/vacuum",
       setup_only=True
   )


It is also possible to analyse existing simulation output directly by passing
the path to a working directory to :class:`FreeEnergy.Relative.analyse <BioSimSpace.FreeEnergy.Relative.analyse>`:

.. code-block:: python

   pmf_vacuum, overlap_vacuum = BSS.FreeEnergy.Relative.analyse("ligA_ligB/vacuum")

simulations, it is also possible to use
:class:`FreeEnergy.Relative <BioSimSpace.FreeEnergy.Relative>` to setup and run simulations
for minimising or equilibrating structures for each lambda window. See the
:class:`FreeEnergyMinimisation <BioSimSpace.Protocol.FreeEnergyMinimisation>` and
:class:`FreeEnergyEquilibration <BioSimSpace.Protocol.FreeEnergyEquilibration>`
protocols for details. At present, these protocols are only supported when not
using :class:`SOMD <BioSimSpace.Process.Somd>` as the simulation engine.


BioSimSpace.FreeEnergy.AToM
----------------------------

This package contains tools to configure, run, and analyse *relative* free
energy simulations using the *alchemical transfer method*  developed by the
`Gallicchio lab <https://www.compmolbiophysbc.org/atom-openmm>`.

Only available in the *OpenMM* engine, the *alchemical transfer method*
replaces the conventional notion of perturbing between two end states with
a single system containing both the free and bound ligand. The relative free
energy of binding is then associated with the swapping of the bound and free
ligands.

The *alchemical transfer method* has a few advantages over the conventional
approach, mainly arising from its relative simplicity and flexibility. The
method is particularly well-suited to the study of difficult ligand
transformations, such as scaffold-hopping and charge change perturbations.
The presence of both ligands in the same system also replaces the conventional
idea of _legs_, combining free, bound, forward and reverse legs into a
single simulation.

In order to perform a relative free energy calculation using the
*alchemical transfer method*, the user requires a protein and two ligands, as
well as knowledge of any common core shared between the two ligands.
AToM-compatible systems can be created from these elements using the
:class:`FreeEnergy.AToM <BioSimSpace.FreeEnergy.AToMSetup>` class.

.. code-block:: python

   from BioSimSpace.FreeEnergy import AToMSetup

   ...

   # Create an AToM setup object. 'protein', 'ligand1' and 'ligand2' must be
   # BioSimSpace Molecule objects.
   # 'ligand1' is bound in the lambda=0 state, 'ligand2' is bound in the lambda=1 state.
   atm_setup = AToMSetup(protein=protein, ligand1=ligand1, ligand2=ligand2)

   # Now create the BioSimSpace system. Here is where knowledge of the common core is required.
   # ligand1_rigid_core and ligand2_rigid_core are lists of integers, each of length three,
   # which define the indices of the common core atoms in the ligands.
   # Displacement is the desired distance between the centre of masses of the two ligands.
   system, data = atm_setup.prepare(
         ligand1_rigid_core=[1, 2, 3],
         ligand2_rigid_core=[1, 2, 3],
         displacement=22.0
   )

   # The prepare function returns two objects: a prepared BioSimSpace system that is ready
   # for AToM simulation, and a data dictionary containing information relevant to AToM calculations.
   # This dictionary does not need to be kept, as the information is also encoded in the system
   # object, but it may be useful for debugging.


Preparing the system for production runs is slightly more complex than in
the conventional approach, as the system will need to be annealed to an
intermediate lambda value, and then equilibrated at that value. The
:ref:`protocol <ref_protocols>` sub-module contains functionality for
equilibrating and annealing systems for AToM simulations.

Once the production simulations have been completed, the user can analyse
the data using the :func:`analyse <BioSimSpace.FreeEnergy.AToM.analyse>` function.

.. code-block:: python

   from BioSimSpace.FreeEnergy import AToM

   # Analyse the simulation data to get the free energy difference and associated error.
   ddg, error = AToM.analyse("path/to/working/directory")
