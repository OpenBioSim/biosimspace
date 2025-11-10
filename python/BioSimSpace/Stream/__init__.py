######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2025
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BioSimSpace is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""
.. currentmodule:: BioSimSpace.Stream

Functions
=========

.. autosummary::
    :toctree: generated/

    save
    load

Examples
========

Stream a :class:`System <BioSimSpace._SireWrappers.System>` object to and from
file.

.. code-block:: python

   import BioSimSpace as BSS

   files = BSS.IO.expand(BSS.tutorialUrl(), ["ala.top", "ala.crd"], ".bz2")

   # Load a molecular system.
   system0 = BSS.IO.readMolecules(files)

   # Stream to file.
   BSS.Stream.save(system0, "system")

   # Alternatively, stream the system object directly.
   system0.save("system")

   # Stream from file.
   system1 = BSS.Stream.load("system.s3")

   # Check the metadata associated with an existing stream. Streaming is
   # generally backwards compatible, but it can be useful to know check
   # the data if issues occur, or when you require an object from a specific
   # sandpit.
   BSS.Stream.getMetadata("system.bss")
   {'bss_object': 'BioSimSpace._SireWrappers._system.System',
   'bss_version': '2023.2.2',
   'bss_revisionid': '9d46295',
   'sire_version': '2023.3.0.dev',
   'sire_revisionid': '90eef20',
   'sandpit': 'None'}

   # You can also get more extensive metadata regarding the underlying Sire
   # object.
   BSS.Stream.getSireMetadata("system.bss")
  {'Object type(s)': 'SireSystem::System',
   'Created by': 'lester',
   'Creation date': 'Thu Jun 29 09:38:21 2023',
   'Created on': 'porridge',
   'System': {'UNIX system': 'Linux',
    'UNIX release': '6.3.9-arch1-1',
    'UNIX version': '#1 SMP PREEMPT_DYNAMIC Wed, 21 Jun 2023 20:46:20 +0000',
    'UNIX machine': 'x86_64',
    'Compiler': 'GNU C++ (13.1.0)',
    'Wordsize': '64 bit',
    'ByteOrder': 'Little endian',
    'Qt runtime version': '5.15.8',
    'Qt compile version': '5.15.8',
    'Sire compile version': '2023.3.0',
    'Compile flags': '-std=c++17 -DSIRE_HAS_CPP_17 -DSIRE_HAS_CPP_14 -DSIRE_HAS_CPP_1Y -DSIRE_HAS_CPP_11 -Wall -Wno-attributes -pipe -DSIRE_ALWAYS_INLINE=inline -Wno-strict-aliasing -DSIRE_VISIBILITY_AVAILABLE -fvisibility=hidden -fvisibility-inlines-hidden -fopenmp-simd -mavx -DSIRE_USE_AVX  -std=c++17 -DSIRE_HAS_CPP_17 -DSIRE_HAS_CPP_14 -DSIRE_HAS_CPP_1Y -DSIRE_HAS_CPP_11 -O3 -ffast-math -fomit-frame-pointer',
    'Link flags': '-rdynamic -Wl,--no-undefined -lpthread'},
   'Repository': 'git@github.com:OpenBioSim/sire.git | devel',
   'Packed size': '317.173 kB',
   'Unpacked size': '3872.57 kB',
   'Country': 'United Kingdom',
   'Language': 'English'}
"""

import sire as _sr

_sr.use_new_api()
del _sr

from ._stream import *
