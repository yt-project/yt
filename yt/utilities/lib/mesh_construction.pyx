"""
This file contains the ElementMesh classes, which represent the target that the 
rays will be cast at when rendering finite element data. This class handles
the interface between the internal representation of the mesh and the pyembree
representation.

Note - this file is only used for the Embree-accelerated ray-tracer.

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np
cimport cython
from libc.stdlib cimport malloc, free
from libc.math cimport fmax, sqrt
cimport numpy as np

cimport pyembree.rtcore as rtc 
cimport pyembree.rtcore_geometry as rtcg
cimport pyembree.rtcore_ray as rtcr
cimport pyembree.rtcore_geometry_user as rtcgu
from pyembree.rtcore cimport \
    Vertex, \
    Triangle, \
    Vec3f

from mesh_traversal cimport YTEmbreeScene
from mesh_samplers cimport \
    sample_hex, \
    sample_tetra, \
    sample_wedge
from mesh_intersection cimport \
    patchIntersectFunc, \
    patchBoundsFunc
from yt.utilities.exceptions import YTElementTypeNotRecognized

cdef extern from "mesh_triangulation.h":
    enum:
        MAX_NUM_TRI
        
    int HEX_NV
    int HEX_NT
    int TETRA_NV
    int TETRA_NT
    int WEDGE_NV
    int WEDGE_NT
    int triangulate_hex[MAX_NUM_TRI][3]
    int triangulate_tetra[MAX_NUM_TRI][3]
    int triangulate_wedge[MAX_NUM_TRI][3]
    int hex20_faces[6][8]


cdef class LinearElementMesh:
    r'''

    This creates a 1st-order mesh to be ray-traced with embree.
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

        # We need to figure out what kind of elements we've been handed.
        if indices.shape[1] == 8:
            self.vpe = HEX_NV
            self.tpe = HEX_NT
            self.tri_array = triangulate_hex
        elif indices.shape[1] == 6:
            self.vpe = WEDGE_NV
            self.tpe = WEDGE_NT
            self.tri_array = triangulate_wedge
        elif indices.shape[1] == 4:
            self.vpe = TETRA_NV
            self.tpe = TETRA_NT
            self.tri_array = triangulate_tetra
        else:
            raise YTElementTypeNotRecognized(vertices.shape[1], 
                                             indices.shape[1])

        self._build_from_indices(scene, vertices, indices)
        self._set_field_data(scene, data)
        self._set_sampler_type(scene)

    cdef void _build_from_indices(self, YTEmbreeScene scene,
                                  np.ndarray vertices_in,
                                  np.ndarray indices_in):
        cdef int i, j
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
        datac.element_indices = self.element_indices
        datac.tpe = self.tpe
        datac.vpe = self.vpe
        datac.fpe = self.fpe
        self.datac = datac
        
        rtcg.rtcSetUserData(scene.scene_i, self.mesh, &self.datac)

    cdef void _set_sampler_type(self, YTEmbreeScene scene):

        if self.vpe == 4:
            self.filter_func = <rtcg.RTCFilterFunc> sample_tetra
        elif self.vpe == 6:
            self.filter_func = <rtcg.RTCFilterFunc> sample_wedge
        elif self.vpe == 8:
            self.filter_func = <rtcg.RTCFilterFunc> sample_hex
        else:
            raise NotImplementedError("Sampler type not implemented.")

        rtcg.rtcSetIntersectionFilterFunction(scene.scene_i,
                                              self.mesh,
                                              self.filter_func)
        
    def __dealloc__(self):
        free(self.field_data)
        free(self.element_indices)
        free(self.vertices)
        free(self.indices)


cdef class QuadraticElementMesh:
    r'''

    This creates a mesh of quadratic patches corresponding to the faces
    of 2nd-order Lagrange elements for direct rendering via embree.
    Currently, this is implemented for 20-point hexahedral meshes only.

    Parameters
    ----------

    scene : EmbreeScene
        This is the scene to which the constructed patches will be
        added.
    vertices : a np.ndarray of floats. 
        This specifies the x, y, and z coordinates of the vertices in 
        the mesh. This should either have the shape 
        (num_vertices, 3). For example, vertices[2][1] should give the 
        y-coordinate of the 3rd vertex in the mesh.
    indices : a np.ndarray of ints
        This should have the shape (num_elements, 20). Each hex will be
        represented in the scene by 6 bi-quadratic patches. We assume that 
        the node ordering is as defined here: 
        http://homepages.cae.wisc.edu/~tautges/papers/cnmev3.pdf
            
    '''

    cdef Patch* patches
    cdef np.float64_t* vertices
    cdef np.float64_t* field_data
    cdef unsigned int mesh
    # patches per element, vertices per element, and field points per 
    # element, respectively:
    cdef int ppe, vpe, fpe

    def __init__(self, YTEmbreeScene scene,
                 np.ndarray vertices, 
                 np.ndarray indices,
                 np.ndarray field_data):

        # only 20-point hexes are supported right now.
        if indices.shape[1] == 20:
            self.vpe = 20
        else:
            raise NotImplementedError

        self._build_from_indices(scene, vertices, indices, field_data)

    cdef void _build_from_indices(self, YTEmbreeScene scene,
                                  np.ndarray vertices_in,
                                  np.ndarray indices_in,
                                  np.ndarray field_data):
        cdef int i, j, ind, idim
        cdef int ne = indices_in.shape[0]
        cdef int nv = vertices_in.shape[0]
        cdef int npatch = 6*ne;

        cdef unsigned int mesh = rtcgu.rtcNewUserGeometry(scene.scene_i, npatch)
        cdef np.ndarray[np.float64_t, ndim=2] element_vertices
        cdef Patch* patches = <Patch*> malloc(npatch * sizeof(Patch))
        self.vertices = <np.float64_t*> malloc(20 * ne * 3 * sizeof(np.float64_t))
        self.field_data = <np.float64_t*> malloc(20 * ne * sizeof(np.float64_t))

        for i in range(ne):
            element_vertices = vertices_in[indices_in[i]]
            for j in range(20):
                self.field_data[i*20 + j] = field_data[i][j]
                for k in range(3):
                    self.vertices[i*20*3 + j*3 + k] = element_vertices[j][k]

        cdef Patch* patch
        for i in range(ne):  # for each element
            element_vertices = vertices_in[indices_in[i]]
            for j in range(6):  # for each face
                patch = &(patches[i*6+j])
                patch.geomID = mesh
                for k in range(8):  # for each vertex
                    ind = hex20_faces[j][k]
                    for idim in range(3):  # for each spatial dimension (yikes)
                        patch.v[k][idim] = element_vertices[ind][idim]
                patch.vertices = self.vertices + i*20*3
                patch.field_data = self.field_data + i*20

        self.patches = patches
        self.mesh = mesh

        rtcg.rtcSetUserData(scene.scene_i, self.mesh, self.patches)
        rtcgu.rtcSetBoundsFunction(scene.scene_i, self.mesh,
                                   <rtcgu.RTCBoundsFunc> patchBoundsFunc)
        rtcgu.rtcSetIntersectFunction(scene.scene_i, self.mesh,
                                      <rtcgu.RTCIntersectFunc> patchIntersectFunc)

    def __dealloc__(self):
        free(self.patches)
        free(self.vertices)
        free(self.field_data)
