# distutils: include_dirs = LIB_DIR
# distutils: libraries = STD_LIBS
cimport cython
cimport libc.math as math
cimport numpy as np
import numpy as np

np.import_array()

cdef np.float64_t h_eos4(np.float64_t e):
    cdef np.float64_t x
    x = 2.25 * e * e
    return 2.5 * e + x / (1.0 + math.sqrt(x + 1.0))

cdef np.float64_t _lorentz_factor_eos4(
        np.float64_t rho, 
        np.float64_t mx,
        np.float64_t my,
        np.float64_t mz,
        np.float64_t e,
        np.float64_t c
    ):
    cdef np.float64_t u2, c2, vx, vy, vz
    cdef np.float64_t htilde, fac;

    c2 = c * c
    fac = (1.0 / (rho * (h_eos4(e) + 1.0))) ** 2
    u2 = (mx * mx + my * my + mz * mz) * fac
    return math.sqrt(1.0 + u2 / c2)

def gamma_eos4(temp):
    cdef np.float64_t[:] kT = temp.ravel()
    
    out = np.empty_like(kT)
    cdef np.float64_t[:] outp = out.ravel()
    cdef np.float64_t x, c_p, c_v
    cdef int i

    for i in range(outp.shape[0]):
        x = 2.25 * kT[i] / math.sqrt(2.25 * kT[i] * kT[i] + 1.0)
        c_p = 2.5 + x
        c_v = 1.5 + x
        outp[i] = c_p / c_v
    return out

cdef np.float64_t cs_eos4(np.float64_t kT, np.float64_t c):
    cdef np.float64_t hp, cs2;
    hp = h_eos4(kT) + 1.0
    cs2 = kT / (3.0 * hp)
    cs2 *= (5.0 * hp - 8.0 * kT) / (hp - kT)
    return c * math.sqrt(cs2)

def sound_speed_eos4(temp, np.float64_t c):
    cdef np.float64_t[:] kT = temp.ravel()
    out = np.empty_like(kT)
    cdef np.float64_t[:] outp = out.ravel()
    cdef np.float64_t hp, cs2;

    cdef int i
    for i in range(outp.shape[0]):
        outp[i] = cs_eos4(kT[i], c)
    return out

cdef np.float64_t _four_vel(
    np.float64_t rho, 
    np.float64_t kT,
    np.float64_t mi
):
    return mi / (rho * (h_eos4(kT) + 1.0))

def four_velocity_xyz_eos4(mom, dens, temp, np.float64_t c):
    cdef np.float64_t[:] rho = dens.ravel()
    cdef np.float64_t[:] mi = mom.ravel()
    cdef np.float64_t[:] kT = temp.ravel()
    
    out = np.empty_like(dens)
    cdef np.float64_t[:] outp = out.ravel()
    cdef np.float64_t[:] ui
    cdef int i

    for i in range(outp.shape[0]):
        outp[i] = _four_vel(rho[i], kT[i], mi[i])
    return out

def lorentz_factor_eos4(dens, momx, momy, momz, temp, np.float64_t c):
    cdef np.float64_t[:] rho = dens.ravel()
    
    cdef np.float64_t[:] mx = momx.ravel()
    cdef np.float64_t[:] my = momy.ravel()
    cdef np.float64_t[:] mz = momz.ravel()
    cdef np.float64_t[:] kT = temp.ravel()
    
    out = np.empty_like(dens)
    cdef np.float64_t[:] outp = out.ravel()
    cdef np.float64_t lf
    cdef int i

    for i in range(outp.shape[0]):
        outp[i] = _lorentz_factor_eos4(rho[i], mx[i], my[i], mz[i], kT[i], c)
    return out

def velocity_xyz_eos4(dens, momx, momy, momz, temp, momi, np.float64_t c):
    cdef np.float64_t[:] rho = dens.ravel()
    cdef np.float64_t[:] mx = momx.ravel()
    cdef np.float64_t[:] my = momy.ravel()
    cdef np.float64_t[:] mz = momz.ravel()
    cdef np.float64_t[:] kT = temp.ravel()
    cdef np.float64_t[:] mi = momi.ravel()
    
    out = np.empty_like(dens)
    cdef np.float64_t[:] outp = out.ravel()
    cdef np.float64_t lf
    cdef int i

    for i in range(outp.shape[0]):
        lf = _lorentz_factor_eos4(rho[i], mx[i], my[i], mz[i], kT[i], c)
        outp[i] = _four_vel(rho[i], kT[i], mi[i]) / lf
    return out

def density_eos4(dens, momx, momy, momz, temp, np.float64_t c):
    cdef np.float64_t[:] rho = dens.ravel()
    
    cdef np.float64_t[:] mx = momx.ravel()
    cdef np.float64_t[:] my = momy.ravel()
    cdef np.float64_t[:] mz = momz.ravel()
    cdef np.float64_t[:] kT = temp.ravel()
    
    out = np.empty_like(dens)
    cdef np.float64_t[:] outp = out.ravel()
    cdef np.float64_t lf
    cdef int i

    for i in range(outp.shape[0]):
        lf = _lorentz_factor_eos4(rho[i], mx[i], my[i], mz[i], kT[i], c)
        outp[i] = rho[i] / lf
    return out

def pressure_eos4(dens, momx, momy, momz, temp, np.float64_t c):
    cdef np.float64_t[:] rho = dens.ravel()
    cdef np.float64_t[:] mx = momx.ravel()
    cdef np.float64_t[:] my = momy.ravel()
    cdef np.float64_t[:] mz = momz.ravel()
    cdef np.float64_t[:] kT = temp.ravel()
    
    out = np.empty_like(dens)
    cdef np.float64_t[:] outp = out.ravel()
    cdef np.float64_t lf, c2
    cdef int i

    c2 = c * c
    for i in range(outp.shape[0]):
        lf = _lorentz_factor_eos4(rho[i], mx[i], my[i], mz[i], kT[i], c)
        outp[i] = rho[i] / lf * c2 * kT[i]
    return out

def specific_thermal_energy_eos4(dens, momx, momy, momz, temp, np.float64_t c):
    cdef np.float64_t[:] rho = dens.ravel()
    cdef np.float64_t[:] mx = momx.ravel()
    cdef np.float64_t[:] my = momy.ravel()
    cdef np.float64_t[:] mz = momz.ravel()
    cdef np.float64_t[:] kT = temp.ravel()
    
    out = np.empty_like(dens)
    cdef np.float64_t[:] outp = out.ravel()
    cdef np.float64_t lf, p, c2, ht
    cdef int i
        
    c2 = c * c
    for i in range(outp.shape[0]):
        outp[i] = c2 * (h_eos4(kT[i])- kT[i])
    return out

def kinetic_energy_density_eos4(dens, momx, momy, momz, temp, np.float64_t c):
    cdef np.float64_t[:] rho = dens.ravel()
    cdef np.float64_t[:] mx = momx.ravel()
    cdef np.float64_t[:] my = momy.ravel()
    cdef np.float64_t[:] mz = momz.ravel()
    cdef np.float64_t[:] kT = temp.ravel()
    
    out = np.empty_like(dens)
    cdef np.float64_t[:] outp = out.ravel()
    cdef np.float64_t lf, u2, hp, ux, uy, uz, c2
    cdef int i

    c2 = c * c
    for i in range(outp.shape[0]):
        hp = h_eos4(kT[i]) + 1.0
        ux = _four_vel(rho[i], kT[i], mx[i])
        uy = _four_vel(rho[i], kT[i], my[i])
        uz = _four_vel(rho[i], kT[i], mz[i])
        u2 = ux**2 + uy**2 + uz**2
        lf = _lorentz_factor_eos4(rho[i], mx[i], my[i], mz[i], kT[i], c)
        gm1 = u2 / c2 / (lf + 1.0)
        p = rho[i] / lf * c2 * kT[i]
        outp[i] = gm1 * (rho[i] * hp * c2 + p)
    return out

def mach_number_eos4(dens, momx, momy, momz, temp, np.float64_t c):
    cdef np.float64_t[:] rho = dens.ravel()
    cdef np.float64_t[:] mx = momx.ravel()
    cdef np.float64_t[:] my = momy.ravel()
    cdef np.float64_t[:] mz = momz.ravel()
    cdef np.float64_t[:] kT = temp.ravel()
    
    out = np.empty_like(dens)
    cdef np.float64_t[:] outp = out.ravel()
    cdef np.float64_t lf, u2, hp, cs
    cdef int i
    
    for i in range(outp.shape[0]):
        cs = cs_eos4(kT[i], c)
        us = cs / math.sqrt(1.0 - cs**2 / c**2)
        u = math.sqrt(mx[i]**2 + my[i]**2 + mz[i]**2) / (rho[i] * (h_eos4(kT[i]) + 1.0))
        outp[i] = u / us
    return out
