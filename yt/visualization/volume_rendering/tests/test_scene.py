
"""
Test for Volume Rendering Scene, and their movement.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import yt
from yt.testing import fake_random_ds
from yt.visualization.volume_rendering.api import Scene, \
    volume_render, Camera, VolumeSource
import numpy as np

def test_rotation():
    ds = fake_random_ds(64)
    ds2 = fake_random_ds(64)
    dd = ds.sphere(ds.domain_center, ds.domain_width[0] / 2)
    dd2 = ds2.sphere(ds2.domain_center, ds2.domain_width[0] / 2)
    
    im, sc = volume_render(dd, field=('gas', 'density'))
    im.write_png('test.png')
    
    vol = sc.get_source(0)
    tf = vol.transfer_function
    tf.clear()
    mi, ma = dd.quantities.extrema('density')
    mi = np.log10(mi)
    ma = np.log10(ma)
    mi_bound = ((ma-mi)*(0.10))+mi
    ma_bound = ((ma-mi)*(0.90))+mi
    tf.map_to_colormap(mi_bound, ma_bound, scale=0.01, colormap='Blues_r')
    
    vol2 = VolumeSource(dd2, field=('gas', 'density'))
    sc.add_source(vol2)
    
    tf = vol2.transfer_function
    tf.clear()
    mi, ma = dd2.quantities.extrema('density')
    mi = np.log10(mi)
    ma = np.log10(ma)
    mi_bound = ((ma-mi)*(0.10))+mi
    ma_bound = ((ma-mi)*(0.90))+mi
    tf.map_to_colormap(mi_bound, ma_bound,  scale=0.01, colormap='Reds_r')
    sc.render('test_scene.png', clip_ratio=6.0)
    
    nrot = 2 
    for i in range(nrot):
        sc.camera.pitch(2*np.pi/nrot)
        sc.render('test_rot_%04i.png' % i, clip_ratio=6.0)
