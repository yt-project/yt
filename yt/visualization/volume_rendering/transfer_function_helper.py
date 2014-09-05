"""
A helper class to build, display, and modify transfer functions for volume
rendering.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.funcs import mylog
from yt.data_objects.profiles import BinnedProfile1D
from yt.visualization.volume_rendering.api import ColorTransferFunction
from yt.visualization._mpl_imports import FigureCanvasAgg
from matplotlib.figure import Figure
from yt.extern.six.moves import StringIO
import numpy as np


class TransferFunctionHelper(object):

    profiles = None

    def __init__(self, ds):
        r"""A transfer function helper.

        This attempts to help set up a good transfer function by finding
        bounds, handling linear/log options, and displaying the transfer
        function combined with 1D profiles of rendering quantity.

        Parameters
        ----------
        ds: A Dataset instance
            A static output that is currently being rendered. This is used to
            help set up data bounds.

        Notes
        -----
        """
        self.ds = ds
        self.field = None
        self.log = False
        self.tf = None
        self.bounds = None
        self.grey_opacity = True
        self.profiles = {}

    def set_bounds(self, bounds=None):
        """
        Set the bounds of the transfer function.

        Parameters
        ----------
        bounds: array-like, length 2, optional
            A length 2 list/array in the form [min, max]. These should be the
            raw values and not the logarithm of the min and max. If bounds is
            None, the bounds of the data are calculated from all of the data
            in the dataset.  This can be slow for very large datasets.
        """
        if bounds is None:
            bounds = self.ds.all_data().quantities.extrema(self.field)
        self.bounds = bounds

        # Do some error checking.
        assert(len(self.bounds) == 2)
        if self.log:
            assert(self.bounds[0] > 0.0)
            assert(self.bounds[1] > 0.0)
        return

    def set_field(self, field):
        """
        Set the field to be rendered

        Parameters
        ----------
        field: string
            The field to be rendered.
        """
        self.field = field

    def set_log(self, log):
        """
        Set whether or not the transfer function should be in log or linear
        space. Also modifies the ds.field_info[field].take_log attribute to
        stay in sync with this setting.

        Parameters
        ----------
        log: boolean
            Sets whether the transfer function should use log or linear space.
        """
        self.log = log
        self.ds.index
        self.ds._get_field_info(self.field).take_log = log

    def build_transfer_function(self):
        """
        Builds the transfer function according to the current state of the
        TransferFunctionHelper.

        Parameters
        ----------
        None

        Returns
        -------

        A ColorTransferFunction object.

        """
        if self.bounds is None:
            mylog.info('Calculating data bounds. This may take a while.' +
                       '  Set the .bounds to avoid this.')
            self.set_bounds()

        if self.log:
            mi, ma = np.log10(self.bounds[0]), np.log10(self.bounds[1])
        else:
            mi, ma = self.bounds
        self.tf = ColorTransferFunction((mi, ma),
                                        grey_opacity=self.grey_opacity,
                                        nbins=512)
        return self.tf

    def plot(self, fn=None, profile_field=None, profile_weight=None):
        """
        Save the current transfer function to a bitmap, or display
        it inline.

        Parameters
        ----------
        fn: string, optional
            Filename to save the image to. If None, the returns an image
            to an IPython session.

        Returns
        -------

        If fn is None, will return an image to an IPython notebook.

        """
        if self.tf is None:
            self.build_transfer_function()
        tf = self.tf
        if self.log:
            xfunc = np.logspace
            xmi, xma = np.log10(self.bounds[0]), np.log10(self.bounds[1])
        else:
            xfunc = np.linspace
            # Need to strip units off of the bounds to avoid a recursion error
            # in matplotlib 1.3.1
            xmi, xma = [np.float64(b) for b in self.bounds]

        x = xfunc(xmi, xma, tf.nbins)
        y = tf.funcs[3].y
        w = np.append(x[1:]-x[:-1], x[-1]-x[-2])
        colors = np.array([tf.funcs[0].y, tf.funcs[1].y, tf.funcs[2].y,
                           np.ones_like(x)]).T

        fig = Figure(figsize=[6, 3])
        canvas = FigureCanvasAgg(fig)
        ax = fig.add_axes([0.2, 0.2, 0.75, 0.75])
        ax.bar(x, tf.funcs[3].y, w, edgecolor=[0.0, 0.0, 0.0, 0.0],
               log=self.log, color=colors, bottom=[0])

        if profile_field is not None:
            try:
                prof = self.profiles[self.field]
            except KeyError:
                self.setup_profile(profile_field, profile_weight)
                prof = self.profiles[self.field]
            if profile_field not in prof.keys():
                prof.add_fields([profile_field], fractional=False,
                                weight=profile_weight)
            # Strip units, if any, for matplotlib 1.3.1
            xplot = np.array(prof[self.field])
            yplot = np.array(prof[profile_field]*tf.funcs[3].y.max() /
                             prof[profile_field].max())
            ax.plot(xplot, yplot, color='w', linewidth=3)
            ax.plot(xplot, yplot, color='k')

        ax.set_xscale({True: 'log', False: 'linear'}[self.log])
        ax.set_xlim(x.min(), x.max())
        ax.set_xlabel(self.ds._get_field_info(self.field).get_label())
        ax.set_ylabel(r'$\mathrm{alpha}$')
        ax.set_ylim(y.max()*1.0e-3, y.max()*2)

        if fn is None:
            from IPython.core.display import Image
            f = StringIO()
            canvas.print_figure(f)
            f.seek(0)
            img = f.read()
            return Image(img)
        else:
            fig.savefig(fn)

    def setup_profile(self, profile_field=None, profile_weight=None):
        if profile_field is None:
            profile_field = 'cell_volume'
        prof = BinnedProfile1D(self.ds.all_data(), 128, self.field,
                               self.bounds[0], self.bounds[1],
                               log_space=self.log,
                               end_collect=False)
        prof.add_fields([profile_field], fractional=False,
                        weight=profile_weight)
        self.profiles[self.field] = prof
        return
