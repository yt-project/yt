"""
This file contains the ElementMesh, which represents the target that the 
rays will be cast at when rendering finite element data. This class handles
the interface between the internal representation of the mesh and the pyembree
representation.


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport numpy as np
cimport cython
cimport pyembree.rtcore as rtc 
from mesh_traversal cimport YTEmbreeScene
cimport pyembree.rtcore_geometry as rtcg
cimport pyembree.rtcore_ray as rtcr
cimport pyembree.rtcore_geometry_user as rtcgu
from mesh_samplers cimport sample_hex
from mesh_samplers cimport sample_tetra
from mesh_samplers cimport sample_element
from pyembree.rtcore cimport \
    Vertex, \
    Triangle, \
    Vec3f
from mesh_intersection cimport \
    patchIntersectFunc, \
    patchBoundsFunc
from libc.stdlib cimport malloc, free
from libc.math cimport fmax, sqrt
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
    int hex20_faces[6][8]

cdef class ElementMesh:
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

    cdef Vertex* vertices
    cdef Triangle* indices
    cdef unsigned int mesh
    cdef double* field_data
    cdef rtcg.RTCFilterFunc filter_func
    # triangles per element, vertices per element, and field points per 
    # element, respectively:
    cdef int tpe, vpe, fpe
    cdef int[MAX_NUM_TRI][3] tri_array
    cdef int* element_indices
    cdef MeshDataContainer datac

    def __init__(self, YTEmbreeScene scene,
                 np.ndarray vertices, 
                 np.ndarray indices,
                 np.ndarray data):

        # We need now to figure out if we've been handed quads or tetrahedra.
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

        self._build_from_indices(scene, vertices, indices)
        self._set_field_data(scene, data)
        self._set_sampler_type(scene)

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
        cdef Vertex* vertices = <Vertex*> malloc(nv * sizeof(Vertex))
        for i in range(nv):
            vertices[i].x = vertices_in[i, 0]
            vertices[i].y = vertices_in[i, 1]
            vertices[i].z = vertices_in[i, 2]       
        rtcg.rtcSetBuffer(scene.scene_i, mesh, rtcg.RTC_VERTEX_BUFFER,
                          vertices, 0, sizeof(Vertex))

        # now build up the triangles
        cdef Triangle* triangles = <Triangle*> malloc(nt * sizeof(Triangle))
        for i in range(ne):
            for j in range(self.tpe):
                triangles[self.tpe*i+j].v0 = indices_in[i][self.tri_array[j][0]]
                triangles[self.tpe*i+j].v1 = indices_in[i][self.tri_array[j][1]]
                triangles[self.tpe*i+j].v2 = indices_in[i][self.tri_array[j][2]]
        rtcg.rtcSetBuffer(scene.scene_i, mesh, rtcg.RTC_INDEX_BUFFER,
                          triangles, 0, sizeof(Triangle))

        cdef int* element_indices = <int *> malloc(ne * self.vpe * sizeof(int))    
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
        self.fpe = data_in.shape[1]

        cdef double* field_data = <double *> malloc(ne * self.fpe * sizeof(double))

        for i in range(ne):
            for j in range(self.fpe):
                field_data[i*self.fpe+j] = data_in[i][j]

        self.field_data = field_data

        cdef MeshDataContainer datac
        datac.vertices = self.vertices
        datac.indices = self.indices
        datac.field_data = self.field_data
        datac.tpe = self.tpe
        datac.vpe = self.vpe
        self.datac = datac
        
        rtcg.rtcSetUserData(scene.scene_i, self.mesh, &self.datac)

    cdef void _set_sampler_type(self, YTEmbreeScene scene):

        if self.fpe == 1:
            self.filter_func = <rtcg.RTCFilterFunc> sample_element
        elif self.fpe == 4:
            self.filter_func = <rtcg.RTCFilterFunc> sample_tetra
        elif self.fpe == 8:
            self.filter_func = <rtcg.RTCFilterFunc> sample_hex
        else:
            print "Error - sampler type not implemented."
            raise NotImplementedError

        rtcg.rtcSetIntersectionFilterFunction(scene.scene_i,
                                              self.mesh,
                                              self.filter_func)
        
    def __dealloc__(self):
        free(self.field_data)
        free(self.element_indices)
        free(self.vertices)
        free(self.indices)


cdef class Order2ElementMesh:
    r'''

    Currently, we handle non-triangular mesh types by converting them 
    to triangular meshes. This class performs this transformation.
    Currently, this is implemented for 20-point hexahedral meshes only.

    Parameters
    ----------

    scene : EmbreeScene
        This is the scene to which the constructed polygons will be
        added.
    vertices : a np.ndarray of floats. 
        This specifies the x, y, and z coordinates of the vertices in 
        the mesh. This should either have the shape 
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

    cdef Patch* patches
    cdef unsigned int mesh
    # patches per element, vertices per element, and field points per 
    # element, respectively:
    cdef int ppe, vpe, fpe

    def __init__(self, YTEmbreeScene scene,
                 np.ndarray vertices, 
                 np.ndarray indices,
                 np.ndarray data):

        # only 20-point hexes are supported right now.
        if indices.shape[1] == 20:
            self.vpe = 20
        else:
            raise NotImplementedError

        self._build_from_indices(scene, vertices, indices)

    cdef void _build_from_indices(self, YTEmbreeScene scene,
                                  np.ndarray vertices_in,
                                  np.ndarray indices_in):
        cdef int i, j, ind, idim
        cdef int nv = vertices_in.shape[0]
        cdef int ne = indices_in.shape[0]
        cdef int np = 6*ne;

        cdef unsigned int mesh = rtcgu.rtcNewUserGeometry(scene.scene_i, np)
        
        cdef Patch* patches = <Patch*> malloc(np * sizeof(Patch));
        cdef Patch* patch
        for i in range(ne):  # for each element
            for j in range(6):  # for each face
                patch = &(patches[i*6+j])
                patch.geomID = mesh
                for k in range(8):  # for each vertex
                    ind = hex20_faces[j][k]
                    for idim in range(3):  # for each spatial dimension (yikes)
                        patch.v[k][idim] = vertices_in[indices_in[i]][ind][idim]
                self._set_bounding_sphere(patch)

        self.patches = patches
        self.mesh = mesh

        rtcg.rtcSetUserData(scene.scene_i, self.mesh, self.patches)
        rtcgu.rtcSetBoundsFunction(scene.scene_i, self.mesh,
                                   <rtcgu.RTCBoundsFunc> patchBoundsFunc)
        rtcgu.rtcSetIntersectFunction(scene.scene_i, self.mesh, 
                                      <rtcgu.RTCIntersectFunc> patchIntersectFunc)

    cdef void _set_bounding_sphere(self, Patch* patch):

        # set the center to be the centroid of the patch vertices
        cdef int i, j
        for j in range(8):
            for i in range(3):
                patch.center[i] += patch.v[j][i]
        for i in range(3):
            patch.center[i] /= 8.0

        # set the radius to be slightly larger than the distance between
        # the center and the farthest vertex
        cdef float r = 0.0
        cdef float d
        for j in range(8):
            d = 0.0
            for i in range(3):
                d += (patch.v[j][i] - patch.center[i])**2
            d = sqrt(d)
            r = fmax(r, d)
        patch.radius = 1.05*r

    def __dealloc__(self):
        free(self.patches)

