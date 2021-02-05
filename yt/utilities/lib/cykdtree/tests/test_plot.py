import os

from yt.utilities.lib.cykdtree.kdtree import PyKDTree
from yt.utilities.lib.cykdtree.plot import plot2D_serial
from yt.utilities.lib.cykdtree.tests import make_points


def test_plot2D_serial():
    fname_test = "test_plot2D_serial.png"
    pts, le, re, ls = make_points(100, 2)
    tree = PyKDTree(pts, le, re, leafsize=ls)
    axs = plot2D_serial(
        tree, pts, title="Serial Test", plotfile=fname_test, label_boxes=True
    )
    os.remove(fname_test)
    # plot2D_serial(tree, pts, axs=axs)
    del axs
