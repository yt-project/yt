cimport numpy as np

cdef void _cartesian_bounds_of_spherical_element(np.float64_t r_i,
                                            np.float64_t theta_i,
                                            np.float64_t phi_i,
                                            np.float64_t dr_i,
                                            np.float64_t dtheta_i,
                                            np.float64_t dphi_i,
                                            np.float64_t xyz_i[3],
                                            np.float64_t dxyz_i[3]
                                            ) noexcept nogil