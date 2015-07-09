cimport numpy as np
cimport cython
cimport pyembree.rtcore as rtc 
from mesh_traversal cimport YTEmbreeScene
cimport pyembree.rtcore_geometry as rtcg
cimport pyembree.rtcore_ray as rtcr
cimport pyembree.rtcore_geometry_user as rtcgu
from yt.utilities.lib.element_mappings import Q1Sampler3D
from filter_feedback_functions cimport \
    sample_hex
from pyembree.rtcore cimport \
    Vertex, \
    Triangle, \
    Vec3f
from libc.stdlib cimport malloc, free
import numpy as np

cdef extern from "mesh_construction.h":
    enum:
        MAX_NUM_TRI
        
    int HEX_NV
    int HEX_NT
    int TETRA_NV
    int TETRA_NT
    int triangulate_hex[MAX_NUM_TRI][3]
    int triangulate_tetra[MAX_NUM_TRI][3]


cdef class TriangleMesh:
    r'''

    This class constructs a polygon mesh with triangular elements and 
    adds it to the scene. 

    Parameters
    ----------

    scene : YTEmbreeScene
        This is the scene to which the constructed polygons will be
        added.
    vertices : a np.ndarray of floats. 
        This specifies the x, y, and z coordinates of the vertices in 
        the polygon mesh. This should either have the shape 
        (num_triangles, 3, 3), or the shape (num_vertices, 3), depending
        on the value of the `indices` parameter.
    indices : either None, or a np.ndarray of ints
        If None, then vertices must have the shape (num_triangles, 3, 3).
        In this case, `vertices` specifices the coordinates of each
        vertex of each triangle in the mesh, with vertices being 
        duplicated if they are shared between triangles. For example,
        if indices is None, then vertices[2][1][0] should give you 
        the x-coordinate of the 2nd vertex of the 3rd triangle.
        If indices is a np.ndarray, then it must have the shape
        (num_triangles, 3), and `vertices` must have the shape
        (num_vertices, 3). In this case, indices[2][1] tells you 
        the index of the 2nd vertex of the 3rd triangle in `indices`,
        while vertices[5][2] tells you the z-coordinate of the 6th
        vertex in the mesh. Note that the indexing is assumed to be
        zero-based. In this setup, vertices can be shared between
        triangles, and the number of vertices can be less than 3 times
        the number of triangles.
            
    '''

    cdef Vertex* vertices
    cdef Triangle* indices
    cdef unsigned int mesh
    cdef double* field_data
    cdef rtcg.RTCFilterFunc filter_func
    cdef int tpe, vpe
    cdef int[MAX_NUM_TRI][3] tri_array
    cdef long* element_indices
    cdef UserData user_data

    def __init__(self, YTEmbreeScene scene,
                 np.ndarray vertices,
                 np.ndarray indices = None):

        if indices is None:
            self._build_from_flat(scene, vertices)
        else:
            self._build_from_indices(scene, vertices, indices)

    cdef void _build_from_flat(self, YTEmbreeScene scene, 
                               np.ndarray tri_vertices):
        cdef int i, j
        cdef int nt = tri_vertices.shape[0]
        # In this scheme, we don't share any vertices.  This leads to cracks,
        # but also means we have exactly three times as many vertices as
        # triangles.
        cdef unsigned int mesh = rtcg.rtcNewTriangleMesh(scene.scene_i,
                        rtcg.RTC_GEOMETRY_STATIC, nt, nt*3, 1) 
        
        cdef Vertex* vertices = <Vertex*> rtcg.rtcMapBuffer(scene.scene_i, mesh,
                        rtcg.RTC_VERTEX_BUFFER)

        for i in range(nt):
            for j in range(3):
                vertices[i*3 + j].x = tri_vertices[i,j,0]
                vertices[i*3 + j].y = tri_vertices[i,j,1]
                vertices[i*3 + j].z = tri_vertices[i,j,2]
        rtcg.rtcUnmapBuffer(scene.scene_i, mesh, rtcg.RTC_VERTEX_BUFFER)

        cdef Triangle* triangles = <Triangle*> rtcg.rtcMapBuffer(scene.scene_i,
                        mesh, rtcg.RTC_INDEX_BUFFER)
        for i in range(nt):
            triangles[i].v0 = i*3 + 0
            triangles[i].v1 = i*3 + 1
            triangles[i].v2 = i*3 + 2

        rtcg.rtcUnmapBuffer(scene.scene_i, mesh, rtcg.RTC_INDEX_BUFFER)
        self.vertices = vertices
        self.indices = triangles
        self.mesh = mesh

    cdef void _build_from_indices(self, YTEmbreeScene scene,
                                  np.ndarray tri_vertices,
                                  np.ndarray tri_indices):
        cdef int i
        cdef int nv = tri_vertices.shape[0]
        cdef int nt = tri_indices.shape[0]

        cdef unsigned int mesh = rtcg.rtcNewTriangleMesh(scene.scene_i,
                                        rtcg.RTC_GEOMETRY_STATIC, nt, nv, 1)

        # set up vertex and triangle arrays. In this case, we just read
        # them directly from the inputs
        cdef Vertex* vertices = <Vertex*> rtcg.rtcMapBuffer(scene.scene_i, mesh,
                                                    rtcg.RTC_VERTEX_BUFFER)

        for i in range(nv):
                vertices[i].x = tri_vertices[i, 0]
                vertices[i].y = tri_vertices[i, 1]
                vertices[i].z = tri_vertices[i, 2]

        rtcg.rtcUnmapBuffer(scene.scene_i, mesh, rtcg.RTC_VERTEX_BUFFER)

        cdef Triangle* triangles = <Triangle*> rtcg.rtcMapBuffer(scene.scene_i,
                        mesh, rtcg.RTC_INDEX_BUFFER)

        for i in range(nt):
            triangles[i].v0 = tri_indices[i][0]
            triangles[i].v1 = tri_indices[i][1]
            triangles[i].v2 = tri_indices[i][2]

        rtcg.rtcUnmapBuffer(scene.scene_i, mesh, rtcg.RTC_INDEX_BUFFER)

        self.vertices = vertices
        self.indices = triangles
        self.mesh = mesh


cdef class ElementMesh(TriangleMesh):
    r'''

    Currently, we handle non-triangular mesh types by converting them 
    to triangular meshes. This class performs this transformation.
    Currently, this is implemented for hexahedral and tetrahedral
    meshes.

    Parameters
    ----------

    scene : EmbreeScene
        This is the scene to which the constructed polygons will be
        added.
    vertices : a np.ndarray of floats. 
        This specifies the x, y, and z coordinates of the vertices in 
        the polygon mesh. This should either have the shape 
        (num_vertices, 3). For example, vertices[2][1] should give the 
        y-coordinate of the 3rd vertex in the mesh.
    indices : a np.ndarray of ints
        This should either have the shape (num_elements, 4) or 
        (num_elements, 8) for tetrahedral and hexahedral meshes, 
        respectively. For tetrahedral meshes, each element will 
        be represented by four triangles in the scene. For hex meshes,
        each element will be represented by 12 triangles, 2 for each 
        face. For hex meshes, we assume that the node ordering is as
        defined here: 
        http://homepages.cae.wisc.edu/~tautges/papers/cnmev3.pdf
            
    '''

    def __init__(self, YTEmbreeScene scene,
                 np.ndarray vertices, 
                 np.ndarray indices,
                 np.ndarray data,
                 sampler_type):

        # We need now to figure out if we've been handed quads or tetrahedra.
        # If it's quads, we can build the mesh slightly differently.
        # http://stackoverflow.com/questions/23723993/converting-quadriladerals-in-an-obj-file-into-triangles

        if indices.shape[1] == 8:
            self.vpe = HEX_NV
            self.tpe = HEX_NT
            self.tri_array = triangulate_hex
        elif indices.shape[1] == 4:
            self.vpe = TETRA_NV
            self.tpe = TETRA_NT
            self.tri_array = triangulate_tetra
        else:
            raise NotImplementedError

        self.field_data = NULL
        self.element_indices = NULL
        self._build_from_indices(scene, vertices, indices)
        self._set_field_data(scene, data)
        self._set_sampler_type(scene, sampler_type)

    cdef void _build_from_indices(self, YTEmbreeScene scene,
                                  np.ndarray vertices_in,
                                  np.ndarray indices_in):
        cdef int i, j, ind
        cdef int nv = vertices_in.shape[0]
        cdef int ne = indices_in.shape[0]
        cdef int nt = self.tpe*ne

        cdef unsigned int mesh = rtcg.rtcNewTriangleMesh(scene.scene_i,
                    rtcg.RTC_GEOMETRY_STATIC, nt, nv, 1)

        # first just copy over the vertices
        cdef Vertex* vertices = <Vertex*> rtcg.rtcMapBuffer(scene.scene_i, mesh,
                        rtcg.RTC_VERTEX_BUFFER)

        for i in range(nv):
            vertices[i].x = vertices_in[i, 0]
            vertices[i].y = vertices_in[i, 1]
            vertices[i].z = vertices_in[i, 2]
        rtcg.rtcUnmapBuffer(scene.scene_i, mesh, rtcg.RTC_VERTEX_BUFFER)

        # now build up the triangles
        cdef Triangle* triangles = <Triangle*> rtcg.rtcMapBuffer(scene.scene_i,
                        mesh, rtcg.RTC_INDEX_BUFFER)

        for i in range(ne):
            for j in range(self.tpe):
                triangles[self.tpe*i+j].v0 = indices_in[i][self.tri_array[j][0]]
                triangles[self.tpe*i+j].v1 = indices_in[i][self.tri_array[j][1]]
                triangles[self.tpe*i+j].v2 = indices_in[i][self.tri_array[j][2]]

        rtcg.rtcUnmapBuffer(scene.scene_i, mesh, rtcg.RTC_INDEX_BUFFER)

        cdef long* element_indices = <long *> malloc(ne * self.vpe * sizeof(long))
    
        for i in range(ne):
            for j in range(self.vpe):
                element_indices[i*self.vpe + j] = indices_in[i][j]

        self.element_indices = element_indices
        self.vertices = vertices
        self.indices = triangles
        self.mesh = mesh

    cdef void _set_field_data(self, YTEmbreeScene scene,
                              np.ndarray data_in):

        cdef int ne = data_in.shape[0]
        cdef double* field_data = <double *> malloc(ne * self.vpe * sizeof(double))

        for i in range(ne):
            for j in range(self.vpe):
                field_data[self.vpe*i+j] = data_in[i][j]

        self.field_data = field_data

        cdef UserData user_data
        user_data.vertices = self.vertices
        user_data.indices = self.indices
        user_data.field_data = self.field_data
        user_data.element_indices = self.element_indices
        user_data.tpe = self.tpe
        user_data.vpe = self.vpe
        self.user_data = user_data
        
        rtcg.rtcSetUserData(scene.scene_i, self.mesh, &self.user_data)

    cdef void _set_sampler_type(self, YTEmbreeScene scene, sampler_type):
        if sampler_type == 'surface':
            self.filter_func = <rtcg.RTCFilterFunc> sample_hex
        else:
            print "Error - sampler type not implemented."
            raise NotImplementedError

        rtcg.rtcSetIntersectionFilterFunction(scene.scene_i,
                                              self.mesh,
                                              self.filter_func)

    # @cython.boundscheck(False)
    # @cython.wraparound(False)
    # @cython.cdivision(True)
    # @cython.initializedcheck(False)
    # def sample_at_point(self, double u, double v, int primID):
        
    #     cdef int elemID, faceID
    #     position = self._get_hit_position(u, v, primID)
    #     vertices = np.empty((8, 3), dtype=np.float64)
        
    #     elemID = primID / self.tpe
    #     # faceID = (primID % self.tpe) / 2
        
    #     # faces = np.array([[0, 1, 2, 3],
    #     #                   [4, 5, 6, 7],
    #     #                   [0, 1, 5, 4],
    #     #                   [1, 2, 6, 5],
    #     #                   [0, 3, 7, 4],
    #     #                   [3, 2, 6, 7]])

    #     # locs = faces[faceID]

    #     element_indices = self.element_indices[elemID]
    #     field_data = np.asarray(self.field_data[elemID], dtype=np.float64)

    #     for i in range(8):
    #         vertices[i][0] = self.vertices[element_indices[i]].x
    #         vertices[i][1] = self.vertices[element_indices[i]].y
    #         vertices[i][2] = self.vertices[element_indices[i]].z    
                             
    #     sampler = Q1Sampler3D()
    #     result = sampler.sample_at_real_point(position, vertices, field_data)

    #     return result

    # @cython.boundscheck(False)
    # @cython.wraparound(False)
    # @cython.cdivision(True)
    # @cython.initializedcheck(False)
    # cdef np.ndarray _get_hit_position(self, double u, double v, int primID):

    #     cdef Triangle tri
    #     cdef Vertex v0, v1, v2
    #     cdef int i
    #     position = np.empty(3, dtype=np.float64)
    #     vertices = np.empty((3, 3), dtype=np.float64)
        
    #     tri = self.indices[primID]
    #     v0 = self.vertices[tri.v0]
    #     v1 = self.vertices[tri.v1]
    #     v2 = self.vertices[tri.v2]

    #     vertices[0][0] = v0.x
    #     vertices[0][1] = v0.y
    #     vertices[0][2] = v0.z
        
    #     vertices[1][0] = v1.x
    #     vertices[1][1] = v1.y
    #     vertices[1][2] = v1.z

    #     vertices[2][0] = v2.x
    #     vertices[2][1] = v2.y
    #     vertices[2][2] = v2.z

    #     for i in range(3):
    #         position[i] = vertices[0][i]*(1.0 - u - v) + vertices[1][i]*u + vertices[2][i]*v

    #     return position

    # @cython.boundscheck(False)
    # @cython.wraparound(False)
    # @cython.cdivision(True)
    # @cython.initializedcheck(False)
    # def sample_triangular(self, double u, double v, int primID):

    #     cdef int i, j
    #     cdef double d0, d1, d2

    #     i = primID / self.tpe
    #     j = primID % self.tpe

    #     d0 = self.field_data[i][self.tri_array[j][0]]
    #     d1 = self.field_data[i][self.tri_array[j][1]]
    #     d2 = self.field_data[i][self.tri_array[j][2]]

    #     return d0*(1.0 - u - v) + d1*u + d2*v
        
    def __dealloc__(self):
        if self.field_data is not NULL:
            free(self.field_data)
        if self.element_indices is not NULL:
            free(self.element_indices)
