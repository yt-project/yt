cimport numpy as np
import numpy as np


cdef double dadtau(double aexp_tau,double O_mat_0,double O_vac_0,double O_k_0):
    return ( aexp_tau**3 * (O_mat_0 + O_vac_0*aexp_tau**3 + O_k_0*aexp_tau) )**0.5

cdef double dadt(double aexp_t,double O_mat_0,double O_vac_0,double O_k_0):
    return ( (1./aexp_t)*(O_mat_0 + O_vac_0*aexp_t**3 + O_k_0*aexp_t) )**0.5


cdef step_cosmo(double alpha,double tau,double aexp_tau,double t,double aexp_t,double O_mat_0,double O_vac_0,double O_k_0):
    cdef double dtau,aexp_tau_pre,dt,aexp_t_pre

    dtau = alpha * aexp_tau / dadtau(aexp_tau,O_mat_0,O_vac_0,O_k_0)
    aexp_tau_pre = aexp_tau - dadtau(aexp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.0
    aexp_tau = aexp_tau - dadtau(aexp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
    tau = tau - dtau

    dt = alpha * aexp_t / dadt(aexp_t,O_mat_0,O_vac_0,O_k_0)
    aexp_t_pre = aexp_t - dadt(aexp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.0
    aexp_t = aexp_t - dadt(aexp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
    t = t - dt

    return tau,aexp_tau,t,aexp_t


cpdef friedman(double O_mat_0,double O_vac_0,double O_k_0):
    cdef double alpha=1.e-5,aexp_min=1.e-3,aexp_tau=1.,aexp_t=1.,tau=0.,t=0.
    cdef int nstep=0,ntable=1000,n_out
    cdef np.ndarray[double,mode='c'] t_out=np.zeros([ntable+1]),tau_out=np.zeros([ntable+1])
    cdef double age_tot,delta_tau,next_tau

    while aexp_tau >= aexp_min or aexp_t >= aexp_min:
       nstep = nstep + 1
       tau,aexp_tau,t,aexp_t = step_cosmo(alpha,tau,aexp_tau,t,aexp_t,O_mat_0,O_vac_0,O_k_0)

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
        tau,aexp_tau,t,aexp_t = step_cosmo(alpha,tau,aexp_tau,t,aexp_t,O_mat_0,O_vac_0,O_k_0)

        if tau < next_tau:
            n_out = n_out + 1
            t_out[n_out] = t
            tau_out[n_out] = tau
            next_tau = next_tau + delta_tau/10.

    while aexp_tau >= aexp_min or aexp_t >= aexp_min:
        tau,aexp_tau,t,aexp_t = step_cosmo(alpha,tau,aexp_tau,t,aexp_t,O_mat_0,O_vac_0,O_k_0)

        if tau < next_tau:
            n_out = n_out + 1
            t_out[n_out] = t
            tau_out[n_out] = tau
            next_tau = next_tau + delta_tau

    n_out = ntable
    t_out[n_out] = t
    tau_out[n_out] = tau

    return tau_out,t_out,delta_tau,ntable,age_tot

cpdef get_ramses_ages(np.ndarray[double,mode='c'] tf, 
                     np.ndarray[double,mode='c'] tauf,  
                     double dtau, 
                     double tsim, 
                     double t_scale, 
                     np.ndarray[double,mode='c'] ages, 
                     int nOver2, 
                     int ntot):

    cdef np.ndarray[double,mode='c'] t
    cdef np.ndarray[double,mode='c'] dage
    cdef np.ndarray[int,mode='c'] iage

    dage = 1 + (10*ages/dtau)
    dage = np.minimum(dage, nOver2 + (dage - nOver2)/10.)
    iage = np.array(dage,dtype=np.int32)

    t = (tf[iage]*(ages - tauf[iage - 1]) / (tauf[iage] - tauf[iage - 1]))
    t = t + (tf[iage-1]*(ages-tauf[iage]) / (tauf[iage-1]-tauf[iage]))
    return  (tsim - t)*t_scale
 
