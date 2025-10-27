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

"""Tools for plotting data."""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["plot", "plotContour", "plotOverlapMatrix"]


from os import environ as _environ

from .. import _is_notebook

# Check to see if DISPLAY is set.
if "DISPLAY" in _environ:
    _display = _environ.get("DISPLAY")
else:
    _display = None
del _environ

if _display is not None:
    _has_display = True
    try:
        import matplotlib.pyplot as _plt
        import matplotlib.colors as _colors

        _has_matplotlib = True
    except ImportError:
        _has_matplotlib = False
else:
    if _is_notebook:
        try:
            import matplotlib.pyplot as _plt
            import matplotlib.colors as _colors

            _has_matplotlib = True
        except ImportError:
            _has_matplotlib = False
    else:
        _has_matplotlib = False
        _has_display = False
        # _warn("The DISPLAY environment variable is unset. Plotting functionality disabled!")

del _display

if _has_matplotlib:
    # Define font sizes.
    _SMALL_SIZE = 14
    _MEDIUM_SIZE = 16
    _BIGGER_SIZE = 18

    # Set font sizes.
    _plt.rc("font", size=_SMALL_SIZE)  # controls default text sizes
    _plt.rc("axes", titlesize=_SMALL_SIZE)  # fontsize of the axes title
    _plt.rc("axes", labelsize=_MEDIUM_SIZE)  # fontsize of the x and y labels
    _plt.rc("xtick", labelsize=_SMALL_SIZE)  # fontsize of the tick labels
    _plt.rc("ytick", labelsize=_SMALL_SIZE)  # fontsize of the tick labels
    _plt.rc("legend", fontsize=_SMALL_SIZE)  # legend fontsize
    _plt.rc("figure", titlesize=_BIGGER_SIZE)  # fontsize of the figure title


def plot(
    x=None,
    y=None,
    xerr=None,
    yerr=None,
    xlabel=None,
    ylabel=None,
    logx=False,
    logy=False,
):
    """
    A simple function to create x/y plots with matplotlib.

    Parameters
    ----------

    x : list
        A list of x data values.

    y : list
        A list of y data values.

    xerr : list
        A list of error values for the x data.

    yerr : list
        A list of error values for the y data.

    xlabel : str
        The x axis label string.

    ylabel : str
        The y axis label string.

    logx : bool
        Whether the x axis is logarithmic.

    logy : bool
        Whether the y axis is logarithmic.
    """
    from ..Types._type import Type as _Type
    from warnings import warn as _warn
    from .. import _is_interactive

    # Make sure were running interactively.
    if not _is_interactive:
        _warn("You can only use BioSimSpace.Notebook.plot when running interactively.")
        return None

    # Matplotlib failed to import.
    if not _has_matplotlib and _has_display:
        _warn(
            "BioSimSpace.Notebook.plot is disabled as matplotlib failed "
            "to load. Please check your matplotlib installation."
        )
        return None

    if not isinstance(logx, bool):
        raise TypeError("'logx' must be of type 'bool'.")
    if not isinstance(logy, bool):
        raise TypeError("'logy' must be of type 'bool'.")

    # Whether we need to convert the x and y data to floats.
    is_unit_x = False
    is_unit_y = False

    if x is None:
        if y is None:
            raise ValueError("'y' data must be defined!")

        # No x data, use array index as value.
        x = [x for x in range(0, len(y))]

    else:
        # No y data, we assume that the user wants to plot the x
        # data as a series.
        if y is None:
            y = x
            x = [x for x in range(0, len(y))]

    # The x argument must be a list or tuple of data records.
    if not isinstance(x, (list, tuple)):
        raise TypeError("'x' must be of type 'list'")

    else:
        # Make sure all records are of the same type. Missing data will be
        # None, so find the first unit type.
        for idx, xx in enumerate(x):
            if xx is not None:
                _type = type(xx)
                break
        if not all(isinstance(xx, (_type, type(None))) for xx in x):
            raise TypeError("All 'x' data values must be of same type")

        # Convert int to float.
        if _type is int:
            x = [float(xx) for xx in x]
            _type = float

            try:
                xerr = [float(xx) for xx in xerr]
            except:
                pass

        # Make sure any associated error has the same unit.
        if xerr is not None:
            if not all(isinstance(xx, (_type, type(None))) for xx in xerr):
                raise TypeError("All 'xerr' values must be of same type as x data")

        # Does this type have units?
        if isinstance(x[idx], _Type):
            is_unit_x = True

    # The y argument must be a list or tuple of data records.
    if not isinstance(y, (list, tuple)):
        raise TypeError("'y' must be of type 'list'")

    else:
        # Make sure all records are of the same type. Missing data will be
        # None, so find the first unit type.
        for idx, yy in enumerate(y):
            if yy is not None:
                _type = type(yy)
                break
        if not all(isinstance(yy, (_type, type(None))) for yy in y):
            raise TypeError("All 'y' data values must be of same type")

        # Convert int to float.
        if _type is int:
            y = [float(yy) for yy in y]
            _type = float

            try:
                yerr = [float(yy) for yy in yerr]
            except:
                pass

        # Make sure any associated error has the same unit.
        if yerr is not None:
            if not all(isinstance(yy, (_type, type(None))) for yy in yerr):
                raise TypeError("All 'yerr' values must be of same type as y data")

        # Does this type have units?
        if isinstance(y[idx], _Type):
            is_unit_y = True

    # Strip any missing values.

    # x-dimension
    idx = [i for i, v in enumerate(x) if v is not None]
    x = list(filter(lambda v: v is not None, x))
    y = [y[i] for i in idx]
    if xerr is not None:
        xerr = [xerr[i] for i in idx]
    if yerr is not None:
        yerr = [yerr[i] for i in idx]

    # y-dimension
    idx = [i for i, v in enumerate(y) if v is not None]
    y = list(filter(lambda v: v is not None, y))
    x = [x[i] for i in idx]
    if xerr is not None:
        xerr = [xerr[i] for i in idx]
    if yerr is not None:
        yerr = [yerr[i] for i in idx]

    # Lists must contain the same number of records.
    # Truncate the longer list to the length of the shortest.
    if len(x) != len(y):
        _warn("Mismatch in list sizes: len(x) = %d, len(y) = %d" % (len(x), len(y)))

        len_x = len(x)
        len_y = len(y)

        if len_x < len_y:
            y = y[:len_x]
        else:
            x = x[:len_y]

        if xerr is not None:
            xerr = xerr[: len(x)]
        if yerr is not None:
            yerr = yerr[: len(y)]

    if xlabel is not None:
        if not isinstance(xlabel, str):
            raise TypeError("'xlabel' must be of type 'str'")
    else:
        if isinstance(x[0], _Type):
            xlabel = (
                x[0].__class__.__qualname__
                + " ("
                + x[0]._print_format[x[0].unit()]
                + ")"
            )

    if ylabel is not None:
        if not isinstance(ylabel, str):
            raise TypeError("'ylabel' must be of type 'str'")
    else:
        if isinstance(y[0], _Type):
            ylabel = (
                y[0].__class__.__qualname__
                + " ("
                + y[0]._print_format[y[0].unit()]
                + ")"
            )

    # Convert the x and y values to floats.
    if is_unit_x:
        x = [x.value() for x in x]
        if xerr is not None:
            xerr = [x.value() for x in xerr]
    if is_unit_y:
        y = [y.value() for y in y]
        if yerr is not None:
            yerr = [y.value() for y in yerr]

    # Set the figure size.
    _plt.figure(figsize=(8, 6))

    # Create the plot.
    if xerr is None and yerr is None:
        _plt.plot(x, y, "-bo")
    else:
        if xerr is None:
            _plt.errorbar(x, y, yerr=yerr, fmt="-bo")
        else:
            if yerr is None:
                _plt.errorbar(x, y, xerr=xerr, fmt="-bo")
            else:
                _plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt="-bo")

    # Add axis labels.
    if xlabel is not None:
        _plt.xlabel(xlabel)
    if ylabel is not None:
        _plt.ylabel(ylabel)

    # Scale the axes.
    if logx:
        _plt.xscale("log")
    if logy:
        _plt.yscale("log")

    # Turn on grid.
    _plt.grid()

    return _plt.show()


def plotContour(x, y, z, xlabel=None, ylabel=None, zlabel=None):
    """
    A simple function to create two-dimensional contour plots with matplotlib.

    Parameters
    ----------

    x : list
        A list of x data values.

    y : list
        A list of y data values.

    z : list
        A list of z data values.

    xlabel : str
        The x axis label string.

    ylabel : str
        The y axis label string.

    zlabel : str
        The z axis label string.
    """
    import numpy as _np
    from warnings import warn as _warn
    from ..Types._type import Type as _Type
    from .. import _is_interactive

    import numpy as _np
    import scipy.interpolate as _interp

    from mpl_toolkits.axes_grid1 import make_axes_locatable as _make_axes_locatable

    # Make sure were running interactively.
    if not _is_interactive:
        _warn("You can only use BioSimSpace.Notebook.plot when running interactively.")
        return None

    # Matplotlib failed to import.
    if not _has_matplotlib and _has_display:
        _warn(
            "BioSimSpace.Notebook.plot is disabled as matplotlib failed "
            "to load. Please check your matplotlib installation."
        )
        return None

    # Whether we need to convert the x, y, and z data to floats.
    is_unit_x = False
    is_unit_y = False
    is_unit_z = False

    # The x argument must be a list or tuple of data records.
    if not isinstance(x, (list, tuple)):
        raise TypeError("'x' must be of type 'list'")

    else:
        # Make sure all records are of the same type.
        _type = type(x[0])
        if not all(isinstance(xx, _type) for xx in x):
            raise TypeError("All 'x' data values must be of same type")

        # Convert int to float.
        if _type is int:
            x = [float(xx) for xx in x]
            _type = float

        # Does this type have units?
        if isinstance(x[0], _Type):
            is_unit_x = True

    # The y argument must be a list or tuple of data records.
    if not isinstance(y, (list, tuple)):
        raise TypeError("'y' must be of type 'list'")

    else:
        # Make sure all records are of the same type.
        _type = type(y[0])
        if not all(isinstance(yy, _type) for yy in y):
            raise TypeError("All 'y' data values must be of same type")

        # Convert int to float.
        if _type is int:
            y = [float(yy) for yy in y]
            _type = float

        # Does this type have units?
        if isinstance(y[0], _Type):
            is_unit_y = True

    if not isinstance(z, (list, tuple)):
        raise TypeError("'z' must be of type 'list'")

    else:
        # Make sure all records are of the same type.
        _type = type(z[0])
        if not all(isinstance(zz, _type) for zz in z):
            raise TypeError("All 'z' data values must be of same type")

        # Convert int to float.
        if _type is int:
            z = [float(zz) for zz in z]
            _type = float

        # Does this type have units?
        if isinstance(z[0], _Type):
            is_unit_z = True

    # Lists must contain the same number of records.
    # Truncate the longer list to the length of the shortest.
    if len(x) != len(y) or len(x) != len(z) or len(y) != len(z):
        _warn(
            "Mismatch in list sizes: len(x) = %d, len(y) = %d, len(z) = %d"
            % (len(x), len(y), len(z))
        )

        lens = [len(x), len(y), len(z)]
        min_len = min(lens)

        x = x[:min_len]
        y = y[:min_len]
        z = z[:min_len]

    if xlabel is not None:
        if not isinstance(xlabel, str):
            raise TypeError("'xlabel' must be of type 'str'")
    else:
        if isinstance(x[0], _Type):
            xlabel = (
                x[0].__class__.__qualname__
                + " ("
                + x[0]._print_format[x[0].unit()]
                + ")"
            )

    if ylabel is not None:
        if not isinstance(ylabel, str):
            raise TypeError("'ylabel' must be of type 'str'")
    else:
        if isinstance(y[0], _Type):
            ylabel = (
                y[0].__class__.__qualname__
                + " ("
                + y[0]._print_format[y[0].unit()]
                + ")"
            )

    if zlabel is not None:
        if not isinstance(zlabel, str):
            raise TypeError("'zlabel' must be of type 'str'")
    else:
        if isinstance(z[0], _Type):
            zlabel = (
                z[0].__class__.__qualname__
                + " ("
                + z[0]._print_format[z[0].unit()]
                + ")"
            )

    # Convert the x and y values to floats.
    if is_unit_x:
        x = [x.value() for x in x]
    if is_unit_y:
        y = [y.value() for y in y]
    if is_unit_z:
        z = [z.value() for z in z]

    # Convert to two-dimensional arrays. We don't assume the data is on a grid,
    # so we interpolate the z values.
    try:
        (X, Y) = _np.meshgrid(
            _np.linspace(_np.min(x), _np.max(x), 1000),
            _np.linspace(_np.min(y), _np.max(y), 1000),
        )
        Z = _interp.griddata((x, y), z, (X, Y), method="linear")
    except:
        raise ValueError("Unable to interpolate x, y, and z data to a grid.")

    # Set the figure size.
    _plt.figure(figsize=(8, 8))

    # Create the contour plot.
    cp = _plt.contourf(X, Y, Z)

    # Add axis labels.
    if xlabel is not None:
        _plt.xlabel(xlabel)
    if ylabel is not None:
        _plt.ylabel(ylabel)

    # Get the current axes.
    ax = _plt.gca()

    # Make sure the axes are equal.
    ax.set_aspect("equal", adjustable="box")

    # Make sure the colour bar matches size of the axes.
    divider = _make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)

    # Add a colour bar and label it.
    cbar = _plt.colorbar(cp, cax=cax)
    if zlabel is not None:
        cbar.set_label(zlabel)

    return _plt.show()


def plotOverlapMatrix(
    overlap, continuous_cbar=False, color_bar_cutoffs=[0.03, 0.1, 0.3]
):
    """
    Plot the overlap matrix from a free-energy perturbation analysis.

    Parameters
    ----------

    overlap : List of List of float, or 2D numpy array of float
        The overlap matrix.

    continuous_cbar : bool, optional, default=False
        If True, use a continuous colour bar. Otherwise, use a discrete
        set of values defined by the 'color_bar_cutoffs' argument to
        assign a colour to each element in the matrix.

    color_bar_cutoffs : List of float, optional, default=[0.03, 0.1, 0.3]
        The cutoffs to use when assigning a colour to each element in the
        matrix. This is used for both the continuous and discrete color bars.
        Can not contain more than 3 elements.
    """
    import numpy as _np
    import matplotlib.colors as _colors
    from warnings import warn as _warn
    from .. import _is_interactive

    # Make sure were running interactively.
    if not _is_interactive:
        _warn("You can only use BioSimSpace.Notebook.plot when running interactively.")
        return None

    # Matplotlib failed to import.
    if not _has_matplotlib and _has_display:
        _warn(
            "BioSimSpace.Notebook.plot is disabled as matplotlib failed "
            "to load. Please check your matplotlib installation."
        )
        return None

    # Validate the input
    if not isinstance(overlap, (list, tuple, _np.ndarray)):
        raise TypeError(
            "The 'overlap' matrix must be a list of list types, or a numpy array!"
        )

    # Try converting to a NumPy array.
    try:
        overlap = _np.array(overlap)
    except:
        raise TypeError(
            "'overlap' must be of type 'np.matrix',  'np.ndarray', or a list of lists."
        )

    # Store the number of rows.
    num_rows = len(overlap)

    # Check the data in each row.
    for row in overlap:
        if not isinstance(row, (list, tuple, _np.ndarray)):
            raise TypeError("The 'overlap' matrix must be a list of list types!")
        if len(row) != num_rows:
            raise ValueError("The 'overlap' matrix must be square!")
        if not all(isinstance(x, float) for x in row):
            raise TypeError("The 'overlap' matrix must contain 'float' types!")

    # Check the colour bar options
    if not isinstance(continuous_cbar, bool):
        raise TypeError("The 'continuous_cbar' option must be a boolean!")
    if not isinstance(color_bar_cutoffs, (list, tuple, _np.ndarray)):
        raise TypeError(
            "The 'color_bar_cutoffs' option must be a list of floats "
            " or a numpy array when 'continuous_cbar' is False!"
        )
    if not all(isinstance(x, float) for x in color_bar_cutoffs):
        raise TypeError("The 'color_bar_cutoffs' option must be a list of floats!")
    if len(color_bar_cutoffs) > 3:
        raise ValueError(
            "The 'color_bar_cutoffs' option must contain no more than 3 elements!"
        )

    # Add 0 and 1 to the colour bar cutoffs.
    if color_bar_cutoffs is not None:
        color_bounds = [0] + color_bar_cutoffs + [1]

    # Tuple of colours and associated font colours.
    # The last and first colours are for the top and bottom of the scale
    # for the continuous colour bar, but are ignored for the discrete bar.
    all_colors = (
        ("#FBE8EB", "black"),  # Lighter pink
        ("#FFD3E0", "black"),
        ("#88CCEE", "black"),
        ("#78C592", "black"),
        ("#117733", "white"),
        ("#004D00", "white"),
    )  # Darker green

    # Set the colour map.
    if continuous_cbar:
        # Create a color map using the extended palette and positions
        box_colors = [all_colors[i][0] for i in range(len(color_bounds) + 1)]
        cmap = _colors.LinearSegmentedColormap.from_list(
            "CustomMap", list(zip(color_bounds, box_colors))
        )

        # Normalise the same way each time so that plots are always comparable.
        norm = _colors.Normalize(vmin=0, vmax=1)
    else:
        # Throw away the first and last colours.
        box_colors = [colors[0] for colors in all_colors[1:-1]]
        cmap = _colors.ListedColormap(
            [box_colors[i] for i in range(len(color_bounds) - 1)]
        )
        norm = _colors.BoundaryNorm(color_bounds, cmap.N)

    # Create the figure and axis. Use a default size for fewer than 16 windows,
    # otherwise scale the figure size to the number of windows.
    if num_rows < 16:
        fig, ax = _plt.subplots(figsize=(8, 8), dpi=300)
    else:
        fig, ax = _plt.subplots(figsize=(num_rows / 2, num_rows / 2), dpi=300)

    # Create the heatmap. Separate the cells with white lines.
    im = ax.imshow(overlap, cmap=cmap, norm=norm)
    for i in range(num_rows - 1):
        for j in range(num_rows - 1):
            # Make sure these are on the edges of the cells.
            ax.axhline(i + 0.5, color="white", linewidth=0.5)
            ax.axvline(j + 0.5, color="white", linewidth=0.5)

    # Label each cell with the overlap value.
    for i in range(num_rows):
        for j in range(num_rows):
            # Get the text colour based on the overlap value.
            overlap_val = overlap[i][j]
            # Get the index of first color bound greater than the overlap value.
            for idx, bound in enumerate(color_bounds):
                if bound > overlap_val:
                    break
            text_color = all_colors[1:-1][idx - 1][1]
            ax.text(
                j,
                i,
                "{:.2f}".format(overlap[i][j]),
                ha="center",
                va="center",
                fontsize=10,
                color=text_color,
            )

    # Create a colorbar. Reduce the height of the colorbar to match the figure and remove the border.
    if continuous_cbar:
        cbar = ax.figure.colorbar(im, ax=ax, cmap=cmap, norm=norm, shrink=0.7)
    else:
        cbar = ax.figure.colorbar(
            im,
            ax=ax,
            cmap=cmap,
            norm=norm,
            boundaries=color_bounds,
            ticks=color_bounds,
            shrink=0.7,
        )
    cbar.outline.set_visible(False)

    # Set the axis labels.
    # Set the x axis at the top of the plot.
    _plt.xlabel(r"$\lambda$ Index")
    ax.xaxis.set_label_position("top")
    _plt.ylabel(r"$\lambda$ Index")

    ticks = [x for x in range(0, num_rows)]

    # Set ticks every lambda window.
    _plt.xticks(ticks)
    ax.xaxis.tick_top()
    _plt.yticks(ticks)

    # Remove the borders.
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Create a tight layout to trim whitespace.
    fig.tight_layout()

    return _plt.show()
