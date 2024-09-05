import numpy as np


def _plot2D_root(
    seg,
    pts=None,
    txt=None,
    plotfile=None,
    point_kw=None,
    box_kw=None,
    axs=None,
    subplot_kw=None,
    gridspec_kw=None,
    fig_kw=None,
    save_kw=None,
    title=None,
    xlabel="x",
    ylabel="y",
    label_kw=None,
):
    r"""Plot a 2D kd-tree.

    Args:
        seg (list of np.ndarray): Line segments to plot defining box edges.
        pts (np.ndarray, optional): Points contained by the kdtree. Defaults to
            None if not provided and points are not plotted.
        txt (list of tuples, optional): Each tuple contains the (x, y, string)
            information for text labels to be added to the boxes. Defaults to
            None and text is not added.
        plotfile (:obj:`str`, optional): Full path to file where the plot
            should be saved. If None, the plot is displayed. Defaults to None
        point_kw (:obj:`dict`, optional): Keywords passed directly to
            :func:`matplotlib.pyplot.scatter` for drawing the points. Defaults
            to empty dict.
        box_kw (:obj:`dict`, optional): Keywords passed directly to
            :class:`matplotlib.collections.LineCollection` for drawing the
            leaf boxes. Defaults to empty dict.

        axs (:obj:`matplotlib.pyplot.Axes`, optional): Axes that should be used
            for plotting. Defaults to None and new axes are created.
        subplot_kw (:obj:`dict`, optional): Keywords passed directly to
            :meth:`matplotlib.figure.Figure.add_subplot`. Defaults to {}.
        gridspec_kw (:obj:`dict`, optional): Keywords passed directly to
            :class:`matplotlib.gridspec.GridSpec`. Defaults to empty dict.
        fig_kw (:obj:`dict`, optional): Keywords passed directly to
            :func:`matplotlib.pyplot.figure`. Defaults to empty dict.
        save_kw (:obj:`dict`, optional): Keywords passed directly to
            :func:`matplotlib.pyplot.savefig`. Defaults to empty dict.

        title (:obj:`str`, optional): Title that the plot should be given.
            Defaults to None and no title is displayed.
        xlabel (:obj:`str`, optional): Label for the x-axis. Defaults to 'x'.
        ylabel (:obj:`str`, optional): Label for the y-axis. Defaults to 'y'.
        label_kw (:obj:`dict`, optional): Keywords passed directly to
            :class:`matplotlib.text.Text` when creating box labels. Defaults
            to empty dict.

    Returns:
        :obj:`matplotlib.pyplot.Axes`: Axes containing the plot.

    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    if point_kw is None:
        point_kw = {}
    if box_kw is None:
        box_kw = {}
    if subplot_kw is None:
        subplot_kw = {}
    if gridspec_kw is None:
        gridspec_kw = {}
    if fig_kw is None:
        fig_kw = {}
    if save_kw is None:
        save_kw = {}
    if label_kw is None:
        label_kw = {}
    # Axes creation
    if axs is None:
        plt.close("all")
        fig, axs = plt.subplots(
            subplot_kw=subplot_kw, gridspec_kw=gridspec_kw, **fig_kw
        )

    # Labels
    if title is not None:
        axs.set_title(title)
    axs.set_xlabel(xlabel, **label_kw)
    axs.set_ylabel(ylabel, **label_kw)

    # Plot points
    if isinstance(pts, list):
        for p in pts:
            if p is not None:
                axs.scatter(p[:, 0], p[:, 1], **point_kw)
    elif pts is not None:
        axs.scatter(pts[:, 0], pts[:, 1], **point_kw)

    # Plot boxes
    lc = LineCollection(seg, **box_kw)
    axs.add_collection(lc)

    # Labels
    if txt is not None:
        # label_kw.setdefault('axes', axs)
        label_kw.setdefault("verticalalignment", "bottom")
        label_kw.setdefault("horizontalalignment", "left")
        for t in txt:
            plt.text(*t, **label_kw)

    axs.autoscale()
    axs.margins(0.1)

    # Save
    if plotfile is not None:
        plt.savefig(plotfile, **save_kw)
    else:
        plt.show()

    # Return axes
    return axs


def plot2D_serial(tree, pts=None, label_boxes=False, **kwargs):
    r"""Plot a 2D kd-tree constructed in serial.

    Parameters
    ----------

    tree: :class:`cykdtree.kdtree.PyKDTree`
        kd-tree class.
    pts: np.ndarray, optional
        Points contained by the kdtree.
    label_boxes: bool
        If True, leaves in the tree are labeled with their index. Defaults to False.

    Additional keywords are passed to :func:`cykdtree.plot._plot2D_root`.

    Returns
    -------

    :obj:`matplotlib.pyplot.Axes`: Axes containing the plot.

    """
    # Box edges
    seg = []
    for leaf in tree.leaves:
        le = leaf.left_edge
        re = leaf.right_edge
        # Top
        seg.append(np.array([[le[0], re[1]], [re[0], re[1]]], "float"))
        # Bottom
        seg.append(np.array([[le[0], le[1]], [re[0], le[1]]], "float"))
        # Left
        seg.append(np.array([[le[0], le[1]], [le[0], re[1]]], "float"))
        # Right
        seg.append(np.array([[re[0], le[1]], [re[0], re[1]]], "float"))

    # Labels
    txt = None
    if label_boxes:
        txt = []
        for leaf in tree.leaves:
            txt.append((leaf.left_edge[0], leaf.left_edge[1], "%d" % leaf.id))

    # Return axes
    return _plot2D_root(seg, pts=pts, txt=txt, **kwargs)
