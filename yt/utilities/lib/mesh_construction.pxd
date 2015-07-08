from pyembree.rtcore cimport \
    Vertex, \
    Triangle, \
    Vec3f

ctypedef struct UserData:
    Vertex* vertices
    Triangle* indices
    double[:,:] field_data
    long[:,:] element_indices
    int tpe
    int vpe
