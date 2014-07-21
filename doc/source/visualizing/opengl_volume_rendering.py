from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL.ARB.vertex_buffer_object import *

import sys, time
import numpy as np
import pycuda.driver as cuda_driver
import pycuda.gl as cuda_gl

from yt.visualization.volume_rendering.theia.scene import TheiaScene
from yt.visualization.volume_rendering.theia.algorithms.front_to_back import FrontToBackRaycaster
from yt.visualization.volume_rendering.transfer_functions import ColorTransferFunction
from yt.visualization.color_maps import *

import numexpr as ne

window = None     # Number of the glut window.
rot_enabled = True

#Theia Scene
ts = None

#RAY CASTING values
c_tbrightness = 1.0
c_tdensity = 0.05

output_texture = None # pointer to offscreen render target

leftButton = False
middleButton = False
rightButton = False

#Screen width and height
width = 1024
height = 1024

eyesep = 0.1

(pbo, pycuda_pbo) = [None]*2

def create_PBO(w, h):
    global pbo, pycuda_pbo
    num_texels = w*h
    array = np.zeros((w,h,3),np.uint32)

    pbo = glGenBuffers(1)
    glBindBuffer(GL_ARRAY_BUFFER, pbo)
    glBufferData(GL_ARRAY_BUFFER, array, GL_DYNAMIC_DRAW)
    glBindBuffer(GL_ARRAY_BUFFER, 0)
    pycuda_pbo = cuda_gl.RegisteredBuffer(long(pbo))

def destroy_PBO(self):
    global pbo, pycuda_pbo
    glBindBuffer(GL_ARRAY_BUFFER, long(pbo))
    glDeleteBuffers(1, long(pbo));
    glBindBuffer(GL_ARRAY_BUFFER, 0)
    pbo,pycuda_pbo = [None]*2

#consistent with C initPixelBuffer()
def create_texture(w,h):
    global output_texture
    output_texture = glGenTextures(1)
    glBindTexture(GL_TEXTURE_2D, output_texture)
    # set basic parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
    # buffer data
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
                 w, h, 0, GL_RGBA, GL_UNSIGNED_INT_8_8_8_8, None)

#consistent with C initPixelBuffer()
def destroy_texture():
    global output_texture
    glDeleteTextures(output_texture);
    output_texture = None

def init_gl(w = 512 , h = 512):
    Width, Height = (w, h)

    glClearColor(0.1, 0.1, 0.5, 1.0)
    glDisable(GL_DEPTH_TEST)

    #matrix functions
    glViewport(0, 0, Width, Height)
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    #matrix functions
    gluPerspective(60.0, Width/float(Height), 0.1, 10.0)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)

def resize(Width, Height):
    global width, height
    (width, height) = Width, Height
    glViewport(0, 0, Width, Height)        # Reset The Current Viewport And Perspective Transformation
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(60.0, Width/float(Height), 0.1, 10.0)


def do_tick():
    global time_of_last_titleupdate, frame_counter, frames_per_second
    if ((time.clock () * 1000.0) - time_of_last_titleupdate >= 1000.):
        frames_per_second = frame_counter                   # Save The FPS
        frame_counter = 0  # Reset The FPS Counter
        szTitle = "%d FPS" % (frames_per_second )
        glutSetWindowTitle ( szTitle )
        time_of_last_titleupdate = time.clock () * 1000.0
    frame_counter += 1

oldMousePos = [ 0, 0 ]
def mouseButton( button, mode, x, y ):
	"""Callback function (mouse button pressed or released).

	The current and old mouse positions are stored in
	a	global renderParam and a global list respectively"""

	global leftButton, middleButton, rightButton, oldMousePos

        if button == GLUT_LEFT_BUTTON:
	    if mode == GLUT_DOWN:
	        leftButton = True
            else:
		leftButton = False

        if button == GLUT_MIDDLE_BUTTON:
	    if mode == GLUT_DOWN:
	        middleButton = True
            else:
		middleButton = False

        if button == GLUT_RIGHT_BUTTON:
	    if mode == GLUT_DOWN:
	        rightButton = True
            else:
		rightButton = False

	oldMousePos[0], oldMousePos[1] = x, y
	glutPostRedisplay( )

def mouseMotion( x, y ):
	"""Callback function (mouse moved while button is pressed).

	The current and old mouse positions are stored in
	a	global renderParam and a global list respectively.
	The global translation vector is updated according to
	the movement of the mouse pointer."""

	global ts, leftButton, middleButton, rightButton, oldMousePos
	deltaX = x - oldMousePos[ 0 ]
	deltaY = y - oldMousePos[ 1 ]

	factor = 0.001

	if leftButton == True:
             ts.camera.rotateX( - deltaY * factor)
             ts.camera.rotateY( - deltaX * factor)
	if middleButton == True:
	     ts.camera.translateX( deltaX* 2.0 * factor)
	     ts.camera.translateY( - deltaY* 2.0 * factor)
	if rightButton == True:
	     ts.camera.scale += deltaY * factor

	oldMousePos[0], oldMousePos[1] = x, y
	glutPostRedisplay( )

def keyPressed(*args):
    global c_tbrightness, c_tdensity
    # If escape is pressed, kill everything.
    if args[0] == '\033':
        print 'Closing..'
        destroy_PBOs()
        destroy_texture()
        exit()

    #change the brightness of the scene
    elif args[0] == ']':
        c_tbrightness += 0.025
    elif args[0] == '[':
        c_tbrightness -= 0.025

    #change the density scale
    elif args[0] == ';':
        c_tdensity -= 0.001
    elif args[0] == '\'':
        c_tdensity += 0.001 

def idle():
    glutPostRedisplay()

def display():
    try:
        #process left eye
        process_image()
        display_image()

        glutSwapBuffers()

    except:
        from traceback import print_exc
        print_exc()
        from os import _exit
        _exit(0)

def process(eye = True):
    global ts, pycuda_pbo, eyesep, c_tbrightness, c_tdensity

    ts.get_raycaster().set_opacity(c_tdensity)
    ts.get_raycaster().set_brightness(c_tbrightness)

    dest_mapping = pycuda_pbo.map()
    (dev_ptr, size) = dest_mapping.device_ptr_and_size()
    ts.get_raycaster().surface.device_ptr = dev_ptr
    ts.update()
   # ts.get_raycaster().cast()
    dest_mapping.unmap()


def process_image():
    global output_texture, pbo, width, height
    """ copy image and process using CUDA """
    # run the Cuda kernel
    process()
    # download texture from PBO
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, np.uint64(pbo))
    glBindTexture(GL_TEXTURE_2D, output_texture)

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
                 width, height, 0, GL_RGBA, GL_UNSIGNED_INT_8_8_8_8_REV, None)

def display_image(eye = True):
    global width, height
    """ render a screen sized quad """
    glDisable(GL_DEPTH_TEST)
    glDisable(GL_LIGHTING)
    glEnable(GL_TEXTURE_2D)
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)

    #matix functions should be moved
    glMatrixMode(GL_PROJECTION)
    glPushMatrix()
    glLoadIdentity()
    glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)
    glMatrixMode( GL_MODELVIEW)
    glLoadIdentity()
    glViewport(0, 0, width, height)

    glBegin(GL_QUADS)
    glTexCoord2f(0.0, 0.0)
    glVertex3f(-1.0, -1.0, 0.5)
    glTexCoord2f(1.0, 0.0)
    glVertex3f(1.0, -1.0, 0.5)
    glTexCoord2f(1.0, 1.0)
    glVertex3f(1.0, 1.0, 0.5)
    glTexCoord2f(0.0, 1.0)
    glVertex3f(-1.0, 1.0, 0.5)
    glEnd()

    glMatrixMode(GL_PROJECTION)
    glPopMatrix()

    glDisable(GL_TEXTURE_2D)
    glBindTexture(GL_TEXTURE_2D, 0)
    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0)
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0)


#note we may need to init cuda_gl here and pass it to camera
def main():
    global window, ts, width, height
    (width, height) = (1024, 1024)

    glutInit(sys.argv)
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH )
    glutInitWindowSize(width, height)
    glutInitWindowPosition(0, 0)
    window = glutCreateWindow("Stereo Volume Rendering")


    glutDisplayFunc(display)
    glutIdleFunc(idle)
    glutReshapeFunc(resize)
    glutMouseFunc( mouseButton )
    glutMotionFunc( mouseMotion )
    glutKeyboardFunc(keyPressed)
    init_gl(width, height)

    # create texture for blitting to screen
    create_texture(width, height)

    import pycuda.gl.autoinit
    import pycuda.gl
    cuda_gl = pycuda.gl

    create_PBO(width, height)
    # ----- Load and Set Volume Data -----

    density_grid = np.load("/home/bogert/dd150_log_densities.npy")

    mi, ma= 21.5, 24.5
    bins = 5000
    tf = ColorTransferFunction( (mi, ma), bins)
    tf.map_to_colormap(mi, ma, colormap="algae", scale_func = scale_func)

    ts = TheiaScene(volume = density_grid, raycaster = FrontToBackRaycaster(size = (width, height), tf = tf))

    ts.get_raycaster().set_sample_size(0.01)
    ts.get_raycaster().set_max_samples(5000)
    ts.update()

    glutMainLoop()

def scale_func(v, mi, ma):
    return  np.minimum(1.0, np.abs((v)-ma)/np.abs(mi-ma) + 0.0)

# Print message to console, and kick off the main to get it rolling.
if __name__ == "__main__":
    print "Hit ESC key to quit, 'a' to toggle animation, and 'e' to toggle cuda"
    main()
