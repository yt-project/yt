from yt.visualization.volume_rendering.theia.scene import TheiaScene
from yt.visualization.volume_rendering.theia.algorithms.front_to_back import FrontToBackRaycaster

from yt.visualization.volume_rendering.transfer_functions import ColorTransferFunction
from yt.visualization.color_maps import *

import subprocess as sp

import numpy as np

def scale_func(v, mi, ma):
      return  np.minimum(1.0, (v-mi)/(ma-mi) + 0.0)

bolshoi = "/home/bogert/log_densities_1024.npy"
density_grid = np.load(bolshoi)

ts = TheiaScene(volume = density_grid, raycaster = FrontToBackRaycaster(size = (1920,1080) ))

mi, ma = 0.0, 3.6
bins = 5000
tf = ColorTransferFunction( (mi, ma), bins)
tf.map_to_colormap(0.5, ma, colormap="spring", scale_func = scale_func)
ts.source.raycaster.set_transfer(tf)

ts.source.raycaster.set_density_scale(0.03)
ts.source.raycaster.set_brightness(2.3)
ts.camera.zoom(30.0)

FFMPEG_BIN = "ffmpeg"

pipe = sp.Popen([ FFMPEG_BIN,
        '-y', # (optional) overwrite the output file if it already exists
        '-f', 'rawvideo',
        '-vcodec','rawvideo',
        '-s', '1920x1080', # size of one frame
        '-pix_fmt', 'rgba',
        '-r', '9.97', # frames per second
        '-i', '-', # The input comes from a pipe
        '-an', # Tells FFMPEG not to expect any audio
        #'-vcodec', 'dnxhd', '-b:v', '220M',
        '-vcodec', 'libx264', '-preset', 'ultrafast', '-qp', '0',# '-b:v', '220M',
        '-pix_fmt', 'yuv420p',
        'bolshoiplanck2.mkv' ],
        stdin=sp.PIPE,stdout=sp.PIPE)
		
		
#ts.source.surface.bounds = (1920,1080)
for k in range (0,50) :
    ts.update()
    ts.camera.rotateX(0.01)
    ts.camera.rotateZ(0.01)
    ts.camera.rotateY(0.01)
    ts.camera.zoom(0.01)

    array = ts.source.get_results()

    array.tofile(pipe.stdin)

pipe.terminate()
