from pyembree.rtcore cimport \
    Vertex, \
    Triangle, \
    Vec3f

ctypedef struct MeshDataContainer:
    Vertex* vertices       # array of triangle vertices
    Triangle* indices      # which vertices belong to which triangles
    double* field_data     # the field values at the vertices
    int* element_indices   # which vertices belong to which *element*
    int tpe                # the number of triangles per element
    int vpe                # the number of vertices per element

ctypedef struct Patch:
    Vec3f v[8]
    Vec3f center
    double radius
    unsigned int geomID
