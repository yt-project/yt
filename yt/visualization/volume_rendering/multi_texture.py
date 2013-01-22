"""
A first-pass at a VisVis-based interactive viewer.

Hi Matt,

You could define your MultiTexture object (which would be a Wobject) and use
multiple visvis.textures.TextureObject instances to represent the sub-textures.
This TextureObject class is not a Wobject, but simply of a wrapper around an
OpenGl texture that handles uploading, updating, etc. 

The MultiTexture object you could derive from BaseTexture or maybe even
Texture3D if your data is 3D. Or maybe just from Wobject if this makes more
sense. Also, the visvis.textures.TextureObjectToVisualize might be used instead
of the TextureObject class.

I mentioned a lot of texture classes, let's make a small list:
  * TextureObject - wraps an OpenGl texture
  * TextureObjectToVisualize(TextureObject) - implements intensity scaling
    (what is changed with the clim property)
  * BaseTexture(Wobject) - the base visvis Wobject class to represent textures.
    It *contains* one TextureObjectToVisualize and has several properties to
    influence it's appearance.  * Texture2D(BaseTexture) - Implementation for 2D
    textures.
  * Texture3D(BaseTexture) - Implementation for 3D textures.
  * YourNewMultiTexture(BaseTexture or ...) - contains multiple
    TextureObjectToVisualize or TextureObject instances.

This is I think the most sensible way to go about this. 

Now in the OnDraw method, you can draw the textures one by one using the same
settings. Or if you want a single draw in which you combine the information
from all textures to determine the final color, you should write a new glsl
shader to which you pass all textures.

I hope this helps,
  Almar
"""
from yt.mods import *
from yt.funcs import *

import visvis as vv
import visvis.textures as vvt

import OpenGL.GL as gl
import OpenGL.GL.ARB.shader_objects as gla
import OpenGL.GLU as glu

import numpy as np

class MultipleTexture(vv.Wobject):
    def __init__(self, parent, data, global_size, renderStyle='mip'):
        vv.Wobject.__init__(self, parent)

        self._global_size = global_size
        
        # create colormap
        self._colormap = vvt.Colormap()
        
        # create glsl program for this texture...
        self._program1 = program =  vvt.GlslProgram()
        
        # scale and translation transforms
        self._trafo_scale = vv.Transform_Scale()
        self._trafo_trans = vv.Transform_Translate()
        self.transformations.append(self._trafo_trans)
        self.transformations.append(self._trafo_scale)        

        data = ensure_list(data)
        self._textures = []
        self._quads = {}
        for obj in data:
            tex = vvt.TextureObjectToVisualize(3, obj) 
            self._textures.append(tex)
            tex.SetData(obj)
        self._program1.SetVertexShader(vvt.vshaders['calculateray'])
        self._program1.SetFragmentShader(vvt.fshaders['mip'])
        self._renderStyle = ''
        self.renderStyle = renderStyle
        if not self._renderStyle:
            self.renderStyle = 'mip'
        self._quadlists = None

    def OnDrawShape(self, clr):
        gl.glColor(clr[0], clr[1], clr[2], 1.0)
        self._DrawQuads()

    def OnDraw(self, fast=False):
        # Prepare by setting things to their defaults. This might release some
        # memory so result in a bigger chance that the shader is run in 
        # hardware mode. On ATI, the line and point smoothing should be off
        # if you want to use gl_FragCoord. (Yeah, I do not see the connection
        # either...)
        gl.glPointSize(1)
        gl.glLineWidth(1)
        gl.glDisable(gl.GL_LINE_STIPPLE)
        gl.glDisable(gl.GL_LINE_SMOOTH)
        gl.glDisable(gl.GL_POINT_SMOOTH)
        
        # only draw front-facing parts
        gl.glEnable(gl.GL_CULL_FACE)
        gl.glCullFace(gl.GL_BACK)
        gl.glBlendFunc(gl.GL_ONE, gl.GL_ONE)
        gl.glBlendEquation(gl.GL_MAX)
        gl.glDisable(gl.GL_DEPTH_TEST)
        
        if self._program1.IsUsable():
            self._program1.Enable()

        if fast:
            self._program1.SetUniformf('stepRatio', [0.4])
        else:
            self._program1.SetUniformf('stepRatio', [1.0])

        self._program1.SetUniformi('texture', [0])        
        self._colormap.Enable(1)
        self._program1.SetUniformi('colormap', [1])

        # enable this texture
        t1 = time.time()
        for i,tex in enumerate(self._textures):

            tex.Enable(0)
            
            if not tex._shape:
                continue

            # fragment shader on
                
            # bind texture- and help-textures (create if it does not exist)
            
            # set uniforms: parameters
            shape = tex._shape # as in opengl
            self._program1.SetUniformf('shape',reversed(list(shape)) )
            
            # do the actual drawing

            self._DrawQuads(tex, i)
            tex.Disable()
        
        # clean up
        gl.glFlush()
        t2 = time.time()
        print "Rendering: %0.3e" % (t2-t1)
        self._colormap.Disable()
        self._program1.Disable()
        #
        gl.glBlendEquation(gl.GL_FUNC_ADD)
        gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
        gl.glDisable(gl.GL_CULL_FACE)
        gl.glEnable(gl.GL_LINE_SMOOTH)
        gl.glEnable(gl.GL_POINT_SMOOTH)
        gl.glEnable(gl.GL_DEPTH_TEST)
        
        
    def _DrawQuads(self, tex, tex_id):
        """ Draw the quads of the texture. 
        This is done in a seperate method to reuse code in 
        OnDraw() and OnDrawShape(). """        
        
        # should we draw?
        if not tex._shape:
            return 
        
        # should we create quads?
        if tex_id not in self._quads:
            self._CreateQuads(tex, tex_id)
        
        # get data
        tex_coord, ver_coord, ind = self._quads[tex_id]
        
        # init vertex and texture array
        gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
        gl.glEnableClientState(gl.GL_TEXTURE_COORD_ARRAY)
        gl.glVertexPointerf(ver_coord.data)
        gl.glTexCoordPointerf(tex_coord.data)
        
        # draw
        gl.glDrawElements(gl.GL_QUADS, len(ind), gl.GL_UNSIGNED_BYTE, ind)
        
        # disable vertex array        
        gl.glDisableClientState(gl.GL_VERTEX_ARRAY)
        gl.glDisableClientState(gl.GL_TEXTURE_COORD_ARRAY)
    
    
    def _GetLimits(self):
        """ _GetLimits()
        Get the limits in world coordinates between which the object exists.
        """
        
        x1, x2 = -0.5, -0.5 + self._global_size[2]
        y1, y2 = -0.5, -0.5 + self._global_size[1]
        z1, z2 = -0.5, -0.5 + self._global_size[0]

        return vv.Wobject._GetLimits(self, x1, x2, y1, y2, z1, z2)
    
    
    def _CreateQuads(self, tex, tex_id):
        axes = self.GetAxes()
        if not axes:
            return
        
        # Store daspect so we can detect it changing
        self._daspectStored = axes.daspect
        
        # Note that we could determine the world coordinates and use
        # them directly here. However, the way that we do it now (using
        # the transformations) is to be preferred, because that way the
        # transformations are applied via the ModelView matrix stack,
        # and can easily be made undone in the raycaster.
        # The -0.5 offset is to center pixels/voxels. This works correctly
        # for anisotropic data.

        dr = tex._dataRef

        ss = [s - 1 for s in tex._shape]

        x0 = dr.origin[2] * self._global_size[2] - 0.5
        x1 = x0 + ss[2] * dr.sampling[2] * self._global_size[2]

        y0 = dr.origin[1] * self._global_size[1]
        y1 = y0 + ss[1] * dr.sampling[1] * self._global_size[1]

        z0 = dr.origin[0] * self._global_size[0]
        z1 = z0 + ss[0] * dr.sampling[0] * self._global_size[0]
        
        # prepare texture coordinates
        t0, t1 = 0, 1
        
        # if any axis are flipped, make sure the correct polygons are front
        # facing
        tmp = 1
        for i in axes.daspect:
            if i<0:
                tmp*=-1        
        if tmp==1:
            t0, t1 = t1, t0
            x0, x1 = x1, x0
            y0, y1 = y1, y0
            z0, z1 = z1, z0

        # using glTexCoord* is the same as glMultiTexCoord*(GL_TEXTURE0)
        # Therefore we need to bind the base texture to 0.
        
        # draw. So we draw the six planes of the cube (well not a cube,
        # a 3d rectangle thingy). The inside is only rendered if the 
        # vertex is facing front, so only 3 planes are rendered at a        
        # time...                
        
        tex_coord, ver_coord = vvt.points.Pointset(3), vvt.points.Pointset(3)
        indices = [0,1,2,3, 4,5,6,7, 3,2,6,5, 0,4,7,1, 0,3,5,4, 1,7,6,2]
        
        # bottom
        tex_coord.Append((t0,t0,t0)); ver_coord.Append((x0, y0, z0)) # 0
        tex_coord.Append((t1,t0,t0)); ver_coord.Append((x1, y0, z0)) # 1
        tex_coord.Append((t1,t1,t0)); ver_coord.Append((x1, y1, z0)) # 2
        tex_coord.Append((t0,t1,t0)); ver_coord.Append((x0, y1, z0)) # 3
        # top
        tex_coord.Append((t0,t0,t1)); ver_coord.Append((x0, y0, z1)) # 4    
        tex_coord.Append((t0,t1,t1)); ver_coord.Append((x0, y1, z1)) # 5
        tex_coord.Append((t1,t1,t1)); ver_coord.Append((x1, y1, z1)) # 6
        tex_coord.Append((t1,t0,t1)); ver_coord.Append((x1, y0, z1)) # 7
        
        # Store quads
        self._quads[tex_id] = (tex_coord, ver_coord, np.array(indices,dtype=np.uint8))

def visvis_plot(vp):
    """
    This function accepts a volume rendering object, which it then tosses into
    visvis for plotting.
    """
        
    vp.partition_grids()
    gs = vp.bricks

    mi = min((g.my_data[0].min() for g in gs))
    ma = max((g.my_data[0].max() for g in gs))

    texes = []

    ax = vv.gca()

    for i,g in enumerate(gs):
        ss = ((g.RightEdge - g.LeftEdge) / (np.array(g.my_data[0].shape)-1)).tolist()
        origin = g.LeftEdge.astype("float32").tolist()
        dd = (g.my_data[0].astype("float32") - mi)/(ma - mi)
        dd = np.clip(dd, 0.0, 1.0)
        print ss
        texes.append(vv.Aarray(dd, origin = origin, sampling = ss))

    mtex = MultipleTexture(ax, texes, global_size=vp.pf.domain_dimensions)

    ax.daspectAuto = False
    ax.SetLimits()
    ax.bgcolor = (0,0,0)

    # set camera
    ax.cameraType = '3d'

    # done
    ax.Draw()

    return mtex, ax
