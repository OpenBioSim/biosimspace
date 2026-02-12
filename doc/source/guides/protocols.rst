.. _ref_protocols:

=========
Protocols
=========

.. toctree::
   :maxdepth: 2

BioSimSpace attempts to achieve interoperability by abstracting a range of
common biomolecular simulation tasks through generic *protocols*. These
present the user with a limited set of options which are handled by all of
the external packages which we support. For example, the various protocols
used for molecular dynamics simulations are defined :ref:`here <ref-Protocol>`.

Protocols are translated into a package specific implementation during the
construction of a :ref:`process <ref-Process>` object. It is here that all
of the required input files are created for the specific process, as well
as any command-line arguments that need to be passed to the executable.

Customising protocols
=====================

Often a user may wish to customise the protocols that are provided by
BioSimSpace. For a particular molecular dynamics engine, this can be
done during the instation of a process by using the ``extra_options`` and
``extra_lines`` parameters. These will override any defaults set by the
protocol. Alternatively, after creating a process you can use the
the helper methods, such as ``setConfig``. This is useful if you wish to
make small tweaks to the default protocols. For example, see the
:class:`BioSimSpace.Process.Amber <BioSimSpace.Process.Amber>`
documentation. (Command-line arguments can be configured in a similar way.)
Alternatively, a user can directly generate a completely custom procotol using
:class:`BioSimSpace.Protocol.Custom <BioSimSpace.Protocol.Custom>`, which
takes a list of configuration strings, or the path to a configuration file,
as a constructor argument. This will then need to be passed to the constructor
of the appropriate process object, (e.g. :class:`BioSimSpace.Process.Amber <BioSimSpace.Process.Amber>`
if it were a configuration file for AMBER).

Note that both of the methods above tie you in to a specific simulation
engine, i.e. a script written with a custom protocol is no longer
interoperable since it will only work on a different computer if the *same*
simulation engine happens to be installed there.

Another way to customise the existing protocols is to *override* the
``createConfig`` methods in the classes used for
:data:`configuration file generation <BioSimSpace._Config>`
for the molecular dynamics enegines.

Writing protocols
=================

Occasionally it might be desirable to write a protocol for a *new* simulation
that isn't currently supported by BioSimSpace. To do so requires the following
steps:

1. Create a new :ref:`Protocol <ref-Protocol>` class that accepts keyword
   arguments that are supported across *all* molecular dynamics engines that
   implement the protocol. (It is not necessary that the protocol is implemented
   by *all* of the engines, only that the arguments are supported by those
   that do.) The class should provide functionality for setting and getting
   the configuration options, along with checking the type and value of the
   user input.
2. Update the private ``_generate_config`` method in each of the :ref:`Process <ref-Process>`
   classes that support the protocol and add functionality to the
   ``createConfig`` method to the appropriate
   :data:`configuration file generator <BioSimSpace._Config>` to translate
   the high-level options into package specific configuration files. There
   is no need to make modifications to the process classes that don't
   support the protocol, since a
   :data:`BioSimSpace._Exceptions.IncompatibleError <BioSimSpace._Exceptions.IncompatibleError>`
   will be raised if the user attempts to create a process to implement
   the protocol.
3. Some protocols require specific command-line arguments to be passed to
   engine's executable. If this is the case, update the ``_generate_args``
   method of the process to modify the argument string when the protocol
   is used.

Protocols are also used in several other places within BioSimSpace, other
than for defining specific molecular dynamics processes described above.
In particular:

* Protocols are used to define the way in which a :ref:`molecular parameterisation <ref-Parameters>`
  is performed. The :class:`BioSimSpace.Parameters <BioSimSpace.Parameters>`
  package comes with its own set of `protocols <https://github.com/openbiosim/biosimspace/tree/devel/src/BioSimSpace/Parameters/Protocol>`__,
  which can be customised by the user. To do so, either override the ``run``
  method of an existing protocol, or define an entirely new protocol that
  inherits from the base class (you would need to do this if, for example,
  you were adding support for a new force field).
* Complex, multi-stage, molecular simulation protocols, such as free-energy
  perturbation, are typically executed by a separate package, e.g.
  :class:`BioSimSpace.FreeEnergy <BioSimSpace.FreeEnergy>`. Here a ``run``
  method orchestrates the separate processes that are required by overall
  simulation. This means that the free-energy objects themselves, e.g.
  :class:`BioSimSpace.FreeEnergy.Relative <BioSimSpace.FreeEnergy.Relative>`,
  can be thought of as a high-level protocol for the simulation method, with
  low-level protocols used to implement the specific stages of the simulation.
* :ref:`Metadynamics <ref-Metadynamics>` support is provided via a molecular
  simulation engine patched by `PLUMED <https://www.plumed.org>`__. Since it
  is not possible to use PLUMED in isolation, the ``BioSimSpace.Process.Plumed``
  class is not exposed to the user, instead being indirectly created as a
  member of another molecular simulation process, e.g.
  :class:`BioSimSpace.Process.Gromacs <BioSimSpace.Process.Gromacs>`, when
  a metadynamics protocol is chosen. Modifications for specific metadynamics
  protocols, e.g. adding support for new collective variables, needs to
  be implemented in the ``createConfig`` method of the `BioSimSpace.Process.Plumed <https://github.com/openbiosim/biosimspace/blob/devel/src/BioSimSpace/Process/_plumed.py>`__ class.
