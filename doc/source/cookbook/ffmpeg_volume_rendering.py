#This is an example of how to make videos of 
#uniform grid data using Theia and ffmpeg

#The Scene object to hold the ray caster and view camera
from yt.visualization.volume_rendering.theia.scene import TheiaScene

#GPU based raycasting algorithm to use 
from yt.visualization.volume_rendering.theia.algorithms.front_to_back import FrontToBackRaycaster

#These will be used to define how to color the data
from yt.visualization.volume_rendering.transfer_functions import ColorTransferFunction
from yt.visualization.color_maps import *

#This will be used to launch ffmpeg
import subprocess as sp

#Of course we need numpy for math magic
import numpy as np

#Opacity scaling function
def scale_func(v, mi, ma):
      return  np.minimum(1.0, (v-mi)/(ma-mi) + 0.0)

#load the uniform grid from a numpy array file
bolshoi = "/home/bogert/log_densities_1024.npy"
density_grid = np.load(bolshoi)

#Set the TheiaScene to use the density_grid and 
#setup the raycaster for a resulting 1080p image
ts = TheiaScene(volume = density_grid, raycaster = FrontToBackRaycaster(size = (1920,1080) ))

#the min and max values in the data to color
mi, ma = 0.0, 3.6

#setup colortransferfunction
bins = 5000
tf = ColorTransferFunction( (mi, ma), bins)
tf.map_to_colormap(0.5, ma, colormap="spring", scale_func = scale_func)

#pass the transfer function to the ray caster
ts.source.raycaster.set_transfer(tf)

#Initial configuration for start of video
#set initial opacity and brightness values
#then zoom into the center of the data 30%
ts.source.raycaster.set_opacity(0.03)
ts.source.raycaster.set_brightness(2.3)
ts.camera.zoom(30.0)

#path to ffmpeg executable
FFMPEG_BIN = "/usr/local/bin/ffmpeg"

pipe = sp.Popen([ FFMPEG_BIN,
        '-y', # (optional) overwrite the output file if it already exists
	#This must be set to rawvideo because the image is an array
        '-f', 'rawvideo', 
	#This must be set to rawvideo because the image is an array
        '-vcodec','rawvideo',
	#The size of the image array and resulting video
        '-s', '1920x1080', 
	#This must be rgba to match array format (uint32)
        '-pix_fmt', 'rgba',
	#frame rate of video
        '-r', '29.97', 
        #Indicate that the input to ffmpeg comes from a pipe
        '-i', '-', 
        # Tells FFMPEG not to expect any audio
        '-an', 
        #Setup video encoder
	#Use any encoder you life available from ffmpeg
        '-vcodec', 'libx264', '-preset', 'ultrafast', '-qp', '0',
        '-pix_fmt', 'yuv420p',
        #Name of the output
        'bolshoiplanck2.mkv' ],
        stdin=sp.PIPE,stdout=sp.PIPE)
		
		
#Now we loop and produce 500 frames
for k in range (0,500) :
    #update the scene resulting in a new image
    ts.update()

    #get the image array from the ray caster
    array = ts.source.get_results()

    #send the image array to ffmpeg
    array.tofile(pipe.stdin)

    #rotate the scene by 0.01 rads in x,y & z
    ts.camera.rotateX(0.01)
    ts.camera.rotateZ(0.01)
    ts.camera.rotateY(0.01)

    #zoom in 0.01% for a total of a 5% zoom
    ts.camera.zoom(0.01)


#Close the pipe to ffmpeg
pipe.terminate()
