from pyembree.rtcore cimport \
    Vertex, \
    Triangle, \
    Vec3f
cimport numpy as np

ctypedef struct MeshDataContainer:
    Vertex* vertices       # array of triangle vertices
    Triangle* indices      # which vertices belong to which triangles
    double* field_data     # the field values at the vertices
    int* element_indices   # which vertices belong to which *element*
    int tpe                # the number of triangles per element
    int vpe                # the number of vertices per element
    int fpe                # the number of field values per element

ctypedef struct Patch:
    float[8][3] v
    unsigned int geomID
    np.float64_t* vertices
    np.float64_t* field_data
