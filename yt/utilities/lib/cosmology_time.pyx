# distutils: language = c++
# distutils: libraries = STD_LIBS

cimport numpy as np

import numpy as np

from libc.math cimport sqrt
import cython



@cython.cdivision(True)
cdef inline double _a_dot(double a, double h0, double om_m, double om_l) noexcept:
    om_k = 1.0 - om_m - om_l
    return h0 * a * sqrt(om_m * (a ** -3) + om_k * (a ** -2) + om_l)


@cython.cdivision(True)
cpdef double _a_dot_recip(double a, double h0, double om_m, double om_l):
    return 1. / _a_dot(a, h0, om_m, om_l)


cdef inline double _da_dtau(double a, double h0, double om_m, double om_l) noexcept:
    return a**2 * _a_dot(a, h0, om_m, om_l)

@cython.cdivision(True)
cpdef double _da_dtau_recip(double a, double h0, double om_m, double om_l) noexcept:
    return 1. / _da_dtau(a, h0, om_m, om_l)



def t_frw(ds, z):
    from scipy.integrate import quad
    aexp = 1 / (1 + z)

    h0 = ds.hubble_constant
    om_m = ds.omega_matter
    om_l = ds.omega_lambda
    conv = ds.quan(0.01, "Mpc/km*s").to("Gyr")

    if isinstance(z, (int, float)):
        return ds.quan(
            quad(_a_dot_recip, 0, aexp, args=(h0, om_m, om_l))[0],
            units=conv,
        )

    return ds.arr(
        [quad(_a_dot_recip, 0, a, args=(h0, om_m, om_l))[0] for a in aexp],
        units=conv,
    )

def tau_frw(ds, z):
    from scipy.integrate import quad
    aexp = 1 / (1 + z)

    h0 = ds.hubble_constant
    om_m = ds.omega_matter
    om_l = ds.omega_lambda

    if isinstance(z, (int, float)):
        return quad(_da_dtau_recip, 1, aexp, args=(h0, om_m, om_l))[0]

    return np.asarray(
        [quad(_da_dtau_recip, 1, a, args=(h0, om_m, om_l))[0] for a in aexp],
    )
