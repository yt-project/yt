import numpy as na
from yt.funcs import *
from yt.utilities.parallel_tools.parallel_analysis_interface import \
    ParallelAnalysisInterface
from yt.utilities.amr_kdtree.api import AMRKDTree
from yt.config import ytcfg

my_rank = ytcfg.getint("yt", "__parallel_rank")
nprocs = ytcfg.getint("yt", "__parallel_size")

class StreamLines(ParallelAnalysisInterface):
    r"""A collection of streamlines that flow through the volume

        Examples
        --------

        >>> c = na.array([0.5]*3)
        >>> c = na.array([0.5]*3)
        >>> N = 100
        >>> scale = 1.0
        >>> pos_dx = na.random.random((N,3))*scale-scale/2.
        >>> pos = c+pos_dx

        >>> streamlines = StreamLines(pf,pos,'x-velocity', 'y-velocity', 'z-velocity', length=1.0) 
        >>> streamlines.integrate_through_volume()

        >>> import matplotlib.pylab as pl
        >>> from mpl_toolkits.mplot3d import Axes3D
        >>> fig=pl.figure() 
        >>> ax = Axes3D(fig)
        >>> for stream in streamlines.streamlines:
        >>>     stream = stream[na.all(stream != 0.0, axis=1)]
        >>>     ax.plot3D(stream[:,0], stream[:,1], stream[:,2], alpha=0.1)
        >>> pl.savefig('streamlines.png')
        """
    def __init__(self, pf, positions, xfield, yfield, zfield, volume=None,
                 dx=None, length=None):
        self.pf = pf
        self.start_positions = positions
        self.N = self.start_positions.shape[0]
        self.xfield = xfield
        self.yfield = yfield
        self.zfield = zfield
        if volume is None:
            volume = AMRKDTree(self.pf, fields=[self.xfield,self.yfield,self.zfield],
                            log_fields=[False,False,False], merge_trees=True)
        self.volume = volume
        if dx is None:
            dx = self.pf.h.get_smallest_dx()
        self.dx = dx
        if length is None:
            length = na.max(self.pf.domain_right_edge-self.pf.domain_left_edge)
        self.length = length
        self.steps = int(length/dx)
        self.streamlines = na.zeros((self.N,self.steps,3), dtype='float64')

    def integrate_through_volume(self):
        my_count = 0
        # Initialize streamlines
        for i in range(self.N):
            if i%nprocs != my_rank:
                continue
            self.streamlines[my_count,0,:] = self.start_positions[i]
            my_count+=1
        self.streamlines = self.streamlines[:my_count,:,:]
        
        pbar = get_pbar("Streamlining", self.N)
        for i,stream in enumerate(self.streamlines):
            step = self.steps
            while (step > 1):
                this_brick = self.volume.locate_brick(stream[-step,:])
                step = self._integrate_through_brick(this_brick, stream, step)
            pbar.update(i)
        pbar.finish()

        self._finalize_parallel()
            
    def _finalize_parallel(self):
        if nprocs <= 1: return
        self.streamlines = self._mpi_cat2d_double_array([
            self.streamlines, [self.N,self.steps,3]])

    def _integrate_through_brick(self, node, stream, step, periodic=False):
        my_dx = node.grid.dds[0]
        while (step > 1):
            self.volume.get_brick_data(node)
            brick = node.brick
            stream[-step+1] = stream[-step]
            brick.integrate_streamline(stream[-step+1], my_dx)
            if na.any(stream[-step+1,:] <= self.pf.domain_left_edge) | \
                   na.any(stream[-step+1,:] >= self.pf.domain_right_edge):
                return 0

            if na.any(stream[-step+1,:] < node.l_corner) | \
                   na.any(stream[-step+1,:] >= node.r_corner):
                return step-1
            step -= 1
        return step

    


        
