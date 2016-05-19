"""
Import the components of the volume rendering extension



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
from yt.data_objects.construction_data_containers import YTStreamline
from yt.funcs import get_pbar
from yt.units.yt_array import YTArray
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface, parallel_passthrough
from yt.utilities.amr_kdtree.api import AMRKDTree

def sanitize_length(length, ds):
    # Ensure that lengths passed in with units are returned as code_length
    # magnitudes without units
    if isinstance(length, YTArray):
        return ds.arr(length).in_units('code_length').d
    else:
        return length

class Streamlines(ParallelAnalysisInterface):
    r"""A collection of streamlines that flow through the volume

    The Streamlines object contains a collection of streamlines
    defined as paths that are parallel to a specified vector field.  

    Parameters
    ----------
    ds : `~yt.data_objects.Dataset`
        This is the dataset to streamline
    pos : array_like
        An array of initial starting positions of the streamlines.
    xfield: field, optional
        The x component of the vector field to be streamlined.
        Default:'velocity_x'
    yfield: field, optional
        The y component of the vector field to be streamlined.
        Default:'velocity_y'
    zfield: field, optional
        The z component of the vector field to be streamlined.
        Default:'velocity_z'
    volume : `yt.extensions.volume_rendering.AMRKDTree`, optional
        The volume to be streamlined.  Can be specified for
        finer-grained control, but otherwise will be automatically
        generated.  At this point it must use the AMRKDTree. 
        Default: None
    dx : float, optional
        Optionally specify the step size during the integration.
        Default: minimum dx
    length : float, optional
        Optionally specify the length of integration.  
        Default: np.max(self.ds.domain_right_edge-self.ds.domain_left_edge)
    direction : real, optional
        Specifies the direction of integration.  The magnitude of this
        value has no effect, only the sign.
    get_magnitude: bool, optional
        Specifies whether the Streamlines.magnitudes array should be
        filled with the magnitude of the vector field at each point in
        the streamline.  This seems to be a ~10% hit to performance.
        Default: False
    
    Examples
    --------
    >>> import yt
    >>> import numpy as np
    >>> import matplotlib.pylab as pl
    >>>
    >>> from yt.visualization.api import Streamlines
    >>> from mpl_toolkits.mplot3d import Axes3D
    >>>
    >>> # Load the dataset and set some parameters
    >>> ds = load('IsolatedGalaxy/galaxy0030/galaxy0030')
    >>> c = np.array([0.5]*3)
    >>> N = 100
    >>> scale = 1.0
    >>> pos_dx = np.random.random((N,3))*scale-scale/2.
    >>> pos = c+pos_dx
    >>>
    >>> # Define and construct streamlines
    >>> streamlines = Streamlines(
            ds,pos, 'velocity_x', 'velocity_y', 'velocity_z', length=1.0) 
    >>> streamlines.integrate_through_volume()
    >>>
    >>> # Make a 3D plot of the streamlines and save it to disk
    >>> fig=pl.figure() 
    >>> ax = Axes3D(fig)
    >>> for stream in streamlines.streamlines:
    >>>     stream = stream[np.all(stream != 0.0, axis=1)]
    >>>     ax.plot3D(stream[:,0], stream[:,1], stream[:,2], alpha=0.1)
    >>> pl.savefig('streamlines.png')
    """
    def __init__(self, ds, positions, xfield='velocity_x', yfield='velocity_x',
                 zfield='velocity_x', volume=None,
                 dx=None, length=None, direction=1,
                 get_magnitude=False):
        ParallelAnalysisInterface.__init__(self)
        self.ds = ds
        self.start_positions = sanitize_length(positions, ds)
        self.N = self.start_positions.shape[0]
        # I need a data object to resolve the field names to field tuples
        # via _determine_fields()
        ad = self.ds.all_data()
        self.xfield = ad._determine_fields(xfield)[0]
        self.yfield = ad._determine_fields(yfield)[0]
        self.zfield = ad._determine_fields(zfield)[0]
        self.get_magnitude=get_magnitude
        self.direction = np.sign(direction)
        if volume is None:
            volume = AMRKDTree(self.ds)
            volume.set_fields([self.xfield,self.yfield,self.zfield],
                              [False,False,False],
                              False)
            volume.join_parallel_trees()
        self.volume = volume
        if dx is None:
            dx = self.ds.index.get_smallest_dx()
        self.dx = sanitize_length(dx, ds)
        if length is None:
            length = np.max(self.ds.domain_right_edge-self.ds.domain_left_edge)
        self.length = sanitize_length(length, ds)
        self.steps = int(length/dx)+1
        # Fix up the dx.
        self.dx = 1.0*self.length/self.steps
        self.streamlines = np.zeros((self.N,self.steps,3), dtype='float64')
        self.magnitudes = None
        if self.get_magnitude:
            self.magnitudes = np.zeros((self.N,self.steps), dtype='float64')
        
    def integrate_through_volume(self):
        nprocs = self.comm.size
        my_rank = self.comm.rank
        self.streamlines[my_rank::nprocs,0,:] = \
            self.start_positions[my_rank::nprocs]

        pbar = get_pbar("Streamlining", self.N)
        for i,stream in enumerate(self.streamlines[my_rank::nprocs]):
            thismag = None
            if self.get_magnitude:
                thismag = self.magnitudes[i,:]
            step = self.steps
            while (step > 1):
                this_node = self.volume.locate_node(stream[-step,:])
                step = self._integrate_through_brick(
                    this_node, stream, step, mag=thismag)
            pbar.update(i)
        pbar.finish()
        
        self._finalize_parallel(None)
        self.streamlines = self.ds.arr(self.streamlines, 'code_length')
        if self.get_magnitude:
            self.magnitudes = self.ds.arr(
                self.magnitudes, self.ds.field_info[self.xfield].units)
       
    @parallel_passthrough
    def _finalize_parallel(self,data):
        self.streamlines = self.comm.mpi_allreduce(self.streamlines, op='sum')
        if self.get_magnitude:
            self.magnitudes = self.comm.mpi_allreduce(
                self.magnitudes, op='sum')

    def _integrate_through_brick(self, node, stream, step,
                                 periodic=False, mag=None):
        LE = self.ds.domain_left_edge.d
        RE = self.ds.domain_right_edge.d
        while (step > 1):
            self.volume.get_brick_data(node)
            brick = node.data
            stream[-step+1] = stream[-step]
            if mag is None:
                brick.integrate_streamline(
                    stream[-step+1], self.direction*self.dx, None)
            else:
                marr = [mag]
                brick.integrate_streamline(
                    stream[-step+1], self.direction*self.dx, marr)
                mag[-step+1] = marr[0]

            cur_stream = stream[-step+1, :]
            if np.sum(np.logical_or(cur_stream < LE, cur_stream >= RE)):
                return 0

            nLE = node.get_left_edge()
            nRE = node.get_right_edge()
            if np.sum(np.logical_or(cur_stream < nLE, cur_stream >= nRE)):
                return step-1
            step -= 1
        return step

    def clean_streamlines(self):
        temp = np.empty(self.N, dtype='object')
        temp2 = np.empty(self.N, dtype='object')
        for i,stream in enumerate(self.streamlines):
            mask = np.all(stream != 0.0, axis=1)
            temp[i] = stream[mask]
            temp2[i] = self.magnitudes[i,mask]
        self.streamlines = temp
        self.magnitudes = temp2

    def path(self, streamline_id):
        """
        Returns an YTSelectionContainer1D object defined by a streamline.

        Parameters
        ----------
        streamline_id : int
            This defines which streamline from the Streamlines object
            that will define the YTSelectionContainer1D object.

        Returns
        -------
        An YTStreamline YTSelectionContainer1D object

        Examples
        --------

        >>> from yt.visualization.api import Streamlines
        >>> streamlines = Streamlines(ds, [0.5]*3) 
        >>> streamlines.integrate_through_volume()
        >>> stream = streamlines.path(0)
        >>> matplotlib.pylab.semilogy(stream['t'], stream['Density'], '-x')

        """
        return YTStreamline(self.streamlines[streamline_id], ds=self.ds,
                            length=self.length)
