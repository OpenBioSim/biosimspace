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
.. currentmodule:: BioSimSpace.Notebook

Functions
=========

.. autosummary::
    :toctree: generated/

    plot
    plotContour
    plotOverlapMatrix

Examples
========

Generate a line graph using two lists of data.

.. image:: ../../../doc/source/_static/plot_01.png
   :width: 800px
   :align: center

If no argument is passed for the ``x`` data then each ``y`` data value is
plotted against its list index.

.. image:: ../../../doc/source/_static/plot_02.png
   :width: 800px
   :align: center

Use the ``xlabel`` and ``ylabel`` arguments to add labels to your plots.

.. image:: ../../../doc/source/_static/plot_03.png
   :width: 800px
   :align: center

Error bars can be added using ``xerror`` and ``yerror``.

.. image:: ../../../doc/source/_static/plot_04.png
   :width: 800px
   :align: center

It is possible to generate plots from the output of a real-time simulation.
Where functions return time-series data containing :ref:`ref-Types`, then
axis labels will be automatically generated. (The ``xlabel`` and ``ylabel``
still take precedence.)

.. image:: ../../../doc/source/_static/plot_05.png
   :width: 800px
   :align: center

Classes
=======

.. autosummary::
    :toctree: generated/

    View

Examples
--------

Load and visualise a molecular system.

.. image:: ../../../doc/source/_static/view_01.png
   :width: 800px
   :align: center

Attach a :class:`View <BioSimSpace.Notebook.View>` to a running molecular
dynamics `Process <BioSimSpace.Process>` and visualise the first molecule in
the system from the latest configuration in real-time.

.. image:: ../../../doc/source/_static/view_02.png
   :width: 800px
   :align: center
"""

import sire as _sr

_sr.use_new_api()
del _sr

from ._plot import *
from ._view import *
