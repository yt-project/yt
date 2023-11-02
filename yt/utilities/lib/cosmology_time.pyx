# distutils: language = c++
# distutils: libraries = STD_LIBS

cimport numpy as np

import numpy as np

from libc.math cimport sqrt
import cython
from scipy.integrate import quad



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


cdef double dadtau(double aexp_tau,double O_mat_0,double O_vac_0,double O_k_0) noexcept:
    cdef double aexp_tau3 = aexp_tau * aexp_tau * aexp_tau
    return sqrt( aexp_tau3 * (O_mat_0 + O_vac_0*aexp_tau3 + O_k_0*aexp_tau) )

@cython.cdivision(True)
cdef double dadt(double aexp_t,double O_mat_0,double O_vac_0,double O_k_0) noexcept:
    cdef double aexp_t3 = aexp_t * aexp_t * aexp_t
    return sqrt( (1./aexp_t)*(O_mat_0 + O_vac_0*aexp_t3 + O_k_0*aexp_t) )


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

@cython.cdivision(True)
cdef void step_cosmo(
    const double alpha,
    double &tau,
    double &aexp_tau,
    double &t,
    double &aexp_t,
    const double O_mat_0,
    const double O_vac_0,
    const double O_k_0,
) noexcept:
    dtau = alpha * aexp_tau / dadtau(aexp_tau,O_mat_0,O_vac_0,O_k_0)
    aexp_tau_pre = aexp_tau - dadtau(aexp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.0
    aexp_tau = aexp_tau - dadtau(aexp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
    tau = tau - dtau

    dt = alpha * aexp_t / dadt(aexp_t,O_mat_0,O_vac_0,O_k_0)
    aexp_t_pre = aexp_t - dadt(aexp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.0
    aexp_t = aexp_t - dadt(aexp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
    t = t - dt

cpdef friedman(double O_mat_0,double O_vac_0,double O_k_0):
    cdef double alpha=1e-4 ,aexp_min=1e-3, aexp_tau=1, aexp_t=1, tau=0, t=0
    cdef int nstep=0, ntable=1000, n_out
    cdef np.ndarray[double,mode='c'] t_out = np.zeros([ntable+1])
    cdef np.ndarray[double,mode='c'] tau_out = np.zeros([ntable+1])
    cdef double age_tot, delta_tau, next_tau

    while aexp_tau >= aexp_min or aexp_t >= aexp_min:
        nstep = nstep + 1

        step_cosmo(
            alpha,
            tau,
            aexp_tau,
            t,
            aexp_t,
            O_mat_0,
            O_vac_0,
            O_k_0,
        )

    age_tot=-t
    if nstep < ntable :
        ntable = nstep
        alpha = alpha / 2.

    delta_tau = 20.*tau/ntable/11.

    aexp_tau = 1.
    aexp_t = 1.
    tau = 0.
    t = 0.

    n_out = 0
    t_out[n_out] = t
    tau_out[n_out] = tau

    next_tau = tau + delta_tau/10.

    while n_out < ntable/2 :
        step_cosmo(
            alpha,
            tau,
            aexp_tau,
            t,
            aexp_t,
            O_mat_0,
            O_vac_0,
            O_k_0,
        )

        if tau < next_tau:
            n_out = n_out + 1
            t_out[n_out] = t
            tau_out[n_out] = tau
            next_tau = next_tau + delta_tau/10.

    while aexp_tau >= aexp_min or aexp_t >= aexp_min:
        step_cosmo(
            alpha,
            tau,
            aexp_tau,
            t,
            aexp_t,
            O_mat_0,
            O_vac_0,
            O_k_0,
        )

        if tau < next_tau:
            n_out = n_out + 1
            t_out[n_out] = t
            tau_out[n_out] = tau
            next_tau = next_tau + delta_tau

    n_out = ntable
    t_out[n_out] = t
    tau_out[n_out] = tau

    return tau_out,t_out,delta_tau,ntable,age_tot
