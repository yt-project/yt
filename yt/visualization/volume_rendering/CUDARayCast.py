"""
An attempt at putting the ray-casting operation into CUDA

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys;sys.path.insert(0,'.')

from yt.mods import *
import yt.extensions.HierarchySubset as hs
import numpy as np
import h5py, time

import matplotlib;matplotlib.use("Agg");import pylab

from yt.extensions.volume_rendering.TransferFunction import ColorTransferFunction

if __name__ == "__main__":

    # This is boilerplate code for setting up pycuda
    import pycuda.driver as cuda
    import pycuda.compiler as compiler
    import pycuda.autoinit
    import pycuda.gpuarray as gpuarray
    cuda.init()
    assert (cuda.Device.count() >= 1)

    print "Extracting hierarchy."
    opf = load("/u/ki/mturk/ki05/MSM96-SIM3-restart-J64/DataDump0081.dir/DataDump0081")
    pf = hs.ExtractedParameterFile(opf, 20)

    cpu = {}
    gpu = {}

    print "Reading data."
    #fn = "DataDump0081_partitioned.h5"
    fn = "RedshiftOutput0005_partitioned.h5"
    f = h5py.File("/u/ki/mturk/ki05/%s" % fn)
    cpu['grid_data'] = f["/PGrids/Data"][:].astype("float32")
    cpu['dims'] = f["/PGrids/Dims"][:].astype("int32") - 1
    cpu['left_edge'] = f["/PGrids/LeftEdges"][:].astype("float32")
    cpu['right_edge'] = f["/PGrids/RightEdges"][:].astype("float32")

    print "Constructing transfer function."
    if "Data" in fn:
        mh = np.log10(1.67e-24)
        tf = ColorTransferFunction((7.5+mh, 14.0+mh))
        tf.add_gaussian( 8.25+mh, 0.002, [0.2, 0.2, 0.4, 0.1])
        tf.add_gaussian( 9.75+mh, 0.002, [0.0, 0.0, 0.3, 0.1])
        tf.add_gaussian(10.25+mh, 0.004, [0.0, 0.3, 0.0, 0.1])
        tf.add_gaussian(11.50+mh, 0.005, [1.0, 0.0, 0.0, 0.7])
        tf.add_gaussian(12.75+mh, 0.010, [1.0, 1.0, 1.0, 1.0])
    elif "Red" in fn:
        tf = ColorTransferFunction((-31, -27))
        tf.add_gaussian(-30.0, 0.05, [1.0, 0.0, 0.0, 0.1])
        tf.add_gaussian(-29.5, 0.03, [0.0, 1.0, 0.0, 0.3])
        tf.add_gaussian(-29.0, 0.05, [0.0, 0.0, 1.0, 0.5])
        tf.add_gaussian(-28.5, 0.05, [1.0, 1.0, 1.0, 1.0])
    else: raise RuntimeError

    cpu['ngrids'] = np.array([cpu['dims'].shape[0]], dtype='int32')
    cpu['tf_r'] = tf.red.y.astype("float32")
    cpu['tf_g'] = tf.green.y.astype("float32")
    cpu['tf_b'] = tf.blue.y.astype("float32")
    cpu['tf_a'] = tf.alpha.y.astype("float32")

    cpu['tf_bounds'] = np.array(tf.x_bounds, dtype='float32')

    cpu['v_dir'] = np.array([0.3, 0.5, 0.6], dtype='float32')

    c = np.array([0.47284317, 0.48062515, 0.58282089], dtype='float32')

    print "Getting cutting plane."
    cp = pf.h.cutting(cpu['v_dir'], c)

    W = 2000.0/pf['au']
    W = 0.25
    Nvec = 128
    back_c = c - cp._norm_vec * W
    front_c = c + cp._norm_vec * W

    px, py = np.mgrid[-W:W:Nvec*1j,-W:W:Nvec*1j]
    xv = cp._inv_mat[0,0]*px + cp._inv_mat[0,1]*py + cp.center[0]
    yv = cp._inv_mat[1,0]*px + cp._inv_mat[1,1]*py + cp.center[1]
    zv = cp._inv_mat[2,0]*px + cp._inv_mat[2,1]*py + cp.center[2]
    cpu['v_pos'] = np.array([xv, yv, zv], dtype='float32').transpose()

    cpu['image_r'] = np.zeros((Nvec, Nvec), dtype='float32').ravel()
    cpu['image_g'] = np.zeros((Nvec, Nvec), dtype='float32').ravel()
    cpu['image_b'] = np.zeros((Nvec, Nvec), dtype='float32').ravel()
    cpu['image_a'] = np.zeros((Nvec, Nvec), dtype='float32').ravel()

    print "Generating module"
    source = open("yt/extensions/volume_rendering/_cuda_caster.cu").read()
    mod = compiler.SourceModule(source)
    func = mod.get_function("ray_cast")

    for n, a in cpu.items():
        ss = a.size * a.dtype.itemsize
        print "Allocating %0.3e megabytes for %s" % (ss/(1024*1024.), n)
        gpu[n] = cuda.to_device(a.ravel('F'))
        #pycuda.autoinit.context.synchronize()

    BLOCK_SIZE = 8
    grid_size = Nvec / BLOCK_SIZE

    print "Running ray_cast function."
    t1 = time.time()
    ret = func(gpu['ngrids'],
               gpu['grid_data'],
               gpu['dims'],
               gpu['left_edge'],
               gpu['right_edge'],
               gpu['tf_r'],
               gpu['tf_g'],
               gpu['tf_b'],
               gpu['tf_a'],
               gpu['tf_bounds'],
               gpu['v_dir'],
               gpu['v_pos'],
               gpu['image_r'],
               gpu['image_g'],
               gpu['image_b'],
               gpu['image_a'],
         block=(BLOCK_SIZE,BLOCK_SIZE,1),
         grid=(grid_size, grid_size), time_kernel=True)
    t2 = time.time()
    print "BACK: %0.3e" % (t2-t1)

    mi, ma = 1e300, -1e300
    image = []
    for im in 'rgb':
        ii = 'image_%s' % im
        sh, dtype = cpu[ii].shape, cpu[ii].dtype
        del cpu[ii]
        cpu[ii] = cuda.from_device(gpu[ii], sh, dtype).reshape((Nvec,Nvec))
        mi, ma = min(cpu[ii].min(),mi), max(cpu[ii].max(), ma)
        image.append(cpu[ii])
        print "Min/max of %s %0.3e %0.3e" % (
                im, image[-1].min(), image[-1].max())
        pylab.clf()
        pylab.imshow(image[-1], interpolation='nearest')
        pylab.savefig("/u/ki/mturk/public_html/vr6/%s.png" % (ii))

    image = np.array(image).transpose()
    image = (image - mi) / (ma - mi)
    pylab.clf()
    pylab.imshow(image, interpolation='nearest')
    pylab.savefig("/u/ki/mturk/public_html/vr6/image_rgb.png")
