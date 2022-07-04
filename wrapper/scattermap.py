# https://gist.github.com/lukauskas/f2f43aad6078a8b5d71b986174487b8c

from seaborn.matrix import _HeatMapper
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

# from seaborn.external.six import string_types

from seaborn.utils import despine, axis_ticklabels_overlap, relative_luminance, to_utf8

class _ScatterMapper(_HeatMapper):
    """
    Draw a scattermap plot, similar to heatmap plot, but use scatter dots instead of heatmap
    """

    def __init__(self, data,
                 marker, marker_size,
                 vmin, vmax, cmap, center, robust, cbar, cbar_kws,
                 xticklabels=True, yticklabels=True, mask=None):

        super(_ScatterMapper, self).__init__(
            data, vmin, vmax, cmap, center, robust, cbar=cbar, cbar_kws=cbar_kws,
            xticklabels=xticklabels, yticklabels=yticklabels, mask=mask,
            # Don't support annotation
            annot=False, fmt=None, annot_kws=None,
        )

        self.marker = marker

        if isinstance(marker_size, float) or isinstance(marker_size, int):
            self.marker_size = marker_size
        elif isinstance(marker_size, pd.DataFrame):
            self.marker_size = marker_size.loc[self.data.index, self.data.columns].values
        else:
            self.marker_size = marker_size
        
    def plot(self, ax, cax, kws):
        """Draw the scattermap on the provided Axes."""
        # Remove all the Axes spines
        despine(ax=ax, left=True, bottom=True)

        # Draw the heatmap
        data = self.plot_data

        range_y = np.arange(data.shape[0], dtype=int) + 0.5
        range_x = np.arange(data.shape[1], dtype=int) + 0.5
        x, y = np.meshgrid(range_x, range_y)

        hmap = ax.scatter(x, y,
                          c=data,
                          marker=self.marker,
                          cmap=self.cmap,
                          vmin=self.vmin, vmax=self.vmax,
                          s=self.marker_size, **kws)

        # Set the axis limits
        ax.set(xlim=(0, self.data.shape[1]), ylim=(0, self.data.shape[0]))

        # Possibly add a colorbar
        if self.cbar:
            cb = ax.figure.colorbar(hmap, cax, ax, **self.cbar_kws)
            cb.outline.set_linewidth(0)
            # If rasterized is passed to pcolormesh, also rasterize the
            # colorbar to avoid white lines on the PDF rendering
            if kws.get('rasterized', False):
                cb.solids.set_rasterized(True)

        # Add row and column labels
        if isinstance(self.xticks, str) and self.xticks == "auto":
            xticks, xticklabels = self._s(ax, self.xticklabels, 0)
        else:
            xticks, xticklabels = self.xticks, self.xticklabels

        if isinstance(self.yticks, str) and self.yticks == "auto":
            yticks, yticklabels = self._auto_ticks(ax, self.yticklabels, 1)
        else:
            yticks, yticklabels = self.yticks, self.yticklabels

        ax.set(xticks=xticks, yticks=yticks)
        xtl = ax.set_xticklabels(xticklabels)
        ytl = ax.set_yticklabels(yticklabels, rotation="vertical")

        # Possibly rotate them if they overlap
        ax.figure.draw(ax.figure.canvas.get_renderer())
        if axis_ticklabels_overlap(xtl):
            plt.setp(xtl, rotation="vertical")
        if axis_ticklabels_overlap(ytl):
            plt.setp(ytl, rotation="horizontal")

        # Add the axis labels
        ax.set(xlabel=self.xlabel, ylabel=self.ylabel)

        # Annotate the cells with the formatted values
        if self.annot:
            self._annotate_heatmap(ax, hmap)

        # Invert the y axis to show the plot in matrix form
        ax.invert_yaxis()

def scattermap(data,
               marker='o',
               marker_size=100,
               vmin=None, vmax=None, cmap=None, center=None, robust=False,
               linewidths=0, linecolor="white",
               cbar=True, cbar_kws=None, cbar_ax=None,
               square=False, xticklabels="auto", yticklabels="auto",
               mask=None, ax=None, **kwargs):
    """Plot rectangular data as a color-encoded matrix.
    This function is similar to `sns.heatmap`, as it is an Axes-level function that will draw the
    heatmap into the currently-active Axes if none is provided to the ``ax`` argument.
    The main difference is that instead of drawing an actual heatmap with filled squares,
    this function will use the `plt.scatter` behind the scenes to draw a scatterplot-heatmap.
    The default is set to plot a grid of circles, however this can be changed via `marker`
    parameter.
    Parameters
    ----------
    data : rectangular dataset
        2D dataset that can be coerced into an ndarray. If a Pandas DataFrame
        is provided, the index/column information will be used to label the
        columns and rows.
    marker: string, optional
        Marker to use: any marker that `pyplot.scatter` supports. Defaults to circle.
    marker_size: int or rectangular dataset
        Either an integer to set the marker size of all data points to,
        or a 2D dataset (like in `data`) that sets individual point sizes.
        Defaults to 100.
    vmin, vmax : floats, optional
        Values to anchor the colormap, otherwise they are inferred from the
        data and other keyword arguments.
    cmap : matplotlib colormap name or object, or list of colors, optional
        The mapping from data values to color space. If not provided, the
        default will depend on whether ``center`` is set.
    center : float, optional
        The value at which to center the colormap when plotting divergant data.
        Using this parameter will change the default ``cmap`` if none is
        specified.
    robust : bool, optional
        If True and ``vmin`` or ``vmax`` are absent, the colormap range is
        computed with robust quantiles instead of the extreme values.
    linewidths : float, optional
        Width of the border lines that will surround the markers
    linecolor : color, optional
        Color of the border lines to the markers
    cbar : boolean, optional
        Whether to draw a colorbar.
    cbar_kws : dict of key, value mappings, optional
        Keyword arguments for `fig.colorbar`.
    cbar_ax : matplotlib Axes, optional
        Axes in which to draw the colorbar, otherwise take space from the
        main Axes.
    square : boolean, optional
        If True, set the Axes aspect to "equal" so each cell will be
        square-shaped.
    xticklabels, yticklabels : "auto", bool, list-like, or int, optional
        If True, plot the column names of the dataframe. If False, don't plot
        the column names. If list-like, plot these alternate labels as the
        xticklabels. If an integer, use the column names but plot only every
        n label. If "auto", try to densely plot non-overlapping labels.
    mask : boolean array or DataFrame, optional
        If passed, data will not be shown in cells where ``mask`` is True.
        Cells with missing values are automatically masked.
    ax : matplotlib Axes, optional
        Axes in which to draw the plot, otherwise use the currently-active
        Axes.
    kwargs : other keyword arguments
        All other keyword arguments are passed to ``ax.pcolormesh``.
    Returns
    -------
    ax : matplotlib Axes
        Axes object with the heatmap.
    See also
    --------
    clustermap : Plot a matrix using hierachical clustering to arrange the
                 rows and columns.
    Examples
    --------
    Plot a scattermap for a numpy array:
    .. plot::
        :context: close-figs
        >>> import numpy as np; np.random.seed(0)
        >>> import seaborn as sns; sns.set()
        >>> uniform_data = np.random.rand(10, 12)
        >>> ax = scattermap(uniform_data)
    Draw on white axes
    .. plot::
        :context: close-figs
        >>> uniform_data = np.random.rand(10, 12)
        >>> with sns.axes_style("white"):
        ...     ax = scattermap(uniform_data)
    Change the limits of the scattermap:
    .. plot::
        :context: close-figs
        >>> ax = scattermap(uniform_data, vmin=0, vmax=1)
    Plot a scattermap for data centered on 0 with a diverging colormap:
    .. plot::
        :context: close-figs
        >>> normal_data = np.random.randn(10, 12)
        >>> ax = scattermap(normal_data, center=0)
    Plot a dataframe with meaningful row and column labels:
    .. plot::
        :context: close-figs
        >>> flights = sns.load_dataset("flights")
        >>> flights = flights.pivot("month", "year", "passengers")
        >>> ax = scattermap(flights)
    Add border lines around each glyph:
    .. plot::
        :context: close-figs
        >>> ax = scattermap(flights, linewidths=1, linecolor='black')
    Use a different colormap:
    .. plot::
        :context: close-figs
        >>> ax = scattermap(flights, cmap="YlGnBu")
    Center the colormap at a specific value:
    .. plot::
        :context: close-figs
        >>> ax = scattermap(flights, center=flights.loc["January", 1955])
    Plot every other column label and don't plot row labels:
    .. plot::
        :context: close-figs
        >>> data = np.random.randn(50, 20)
        >>> ax = scattermap(data, xticklabels=2, yticklabels=False)
    Don't draw a colorbar:
    .. plot::
        :context: close-figs
        >>> ax = scattermap(flights, cbar=False)
    Use different axes for the colorbar:
    .. plot::
        :context: close-figs
        >>> grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
        >>> f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws)
        >>> ax = scattermap(flights, ax=ax,
        ...                  cbar_ax=cbar_ax,
        ...                  cbar_kws={"orientation": "horizontal"})
    Use a mask to plot only part of a matrix
    .. plot::
        :context: close-figs
        >>> corr = np.corrcoef(np.random.randn(10, 200))
        >>> mask = np.zeros_like(corr)
        >>> mask[np.triu_indices_from(mask)] = True
        >>> with sns.axes_style("white"):
        ...     ax = scattermap(corr, mask=mask, vmax=.3, square=True)
     Change glyph, plot stars instead of circles
    .. plot::
        :context: close-figs
        >>> ax = scattermap(corr, vmax=.3, square=True, marker='*')
    Plot multiple markers on the same plot
    >>> corr = np.corrcoef(np.random.randn(10, 200))
    >>> mask = np.zeros_like(corr)
    >>> mask[np.triu_indices_from(mask)] = True
    >>> with sns.axes_style("white"):
    ...     ax = scattermap(corr, mask=mask, vmax=.3, square=True)
    ...     ax = scattermap(corr, mask=mask.T, vmax=.3, square=True, ax=ax, marker='*', cbar=False)
    Specify size for points
    .. plot::
        :context: close-figs
    >>> with sns.axes_style("white"):
    ...     ax = scattermap(corr, vmax=.3, square=True, marker_size=np.abs(corr)*300)
    """
    # Initialize the plotter object
    plotter = _ScatterMapper(data,
                             marker, marker_size,
                             vmin, vmax, cmap, center, robust,
                             cbar, cbar_kws, xticklabels,
                             yticklabels, mask)

    # Add the pcolormesh kwargs here
    kwargs["linewidths"] = linewidths
    kwargs["edgecolor"] = linecolor

    # Draw the plot and return the Axes
    if ax is None:
        ax = plt.gca()
    if square:
        ax.set_aspect("equal")
    plotter.plot(ax, cbar_ax, kwargs)
    return ax

        