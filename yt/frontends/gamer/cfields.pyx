# distutils: include_dirs = LIB_DIR
# distutils: libraries = STD_LIBS
cimport cython
cimport libc.math as math
cimport numpy as np

import numpy as np


cdef np.float64_t gamma_eos(np.float64_t kT, np.float64_t g) noexcept nogil:
    return g

cdef np.float64_t gamma_eos_tb(np.float64_t kT, np.float64_t g) noexcept nogil:
    cdef np.float64_t x, c_p, c_v
    x = 2.25 * kT / ( math.sqrt(2.25 * kT * kT + 1.0) + 1.0 )
    c_p = 2.5 + x
    c_v = 1.5 + x
    return c_p / c_v

cdef np.float64_t cs_eos_tb(np.float64_t kT, np.float64_t h, np.float64_t g) noexcept nogil:
    cdef np.float64_t hp, cs2
    hp = h + 1.0
    cs2 = kT / (3.0 * hp)
    cs2 *= (5.0 * hp - 8.0 * kT) / (hp - kT)
    return math.sqrt(cs2)

cdef np.float64_t cs_eos(np.float64_t kT, np.float64_t h, np.float64_t g) noexcept nogil:
    cdef np.float64_t hp, cs2
    hp = h + 1.0
    cs2 = g / hp * kT
    return math.sqrt(cs2)


ctypedef np.float64_t (*f2_type)(np.float64_t, np.float64_t) noexcept nogil
ctypedef np.float64_t (*f3_type)(np.float64_t, np.float64_t, np.float64_t) noexcept nogil


cdef class SRHDFields:
    cdef f2_type gamma
    cdef f3_type cs
    cdef np.float64_t _gamma

    def __init__(self, int eos, np.float64_t gamma):
        self._gamma = gamma
        # Select aux functions based on eos no.
        if eos == 1:
            self.gamma = gamma_eos
            self.cs = cs_eos
        else:
            self.gamma = gamma_eos_tb
            self.cs = cs_eos_tb

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def gamma_field(self, temp):
        cdef np.float64_t[:] kT = temp.ravel()

        out = np.empty_like(kT)
        cdef np.float64_t[:] outp = out.ravel()
        cdef int i

        for i in range(outp.shape[0]):
            outp[i] = self.gamma(kT[i], self._gamma)
        return out

    cdef np.float64_t _lorentz_factor(
            self,
            np.float64_t rho,
            np.float64_t mx,
            np.float64_t my,
            np.float64_t mz,
            np.float64_t h,
        ) noexcept nogil:
        cdef np.float64_t u2
        cdef np.float64_t fac

        fac = (1.0 / (rho * (h + 1.0))) ** 2
        u2 = (mx * mx + my * my + mz * mz) * fac
        return math.sqrt(1.0 + u2)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def lorentz_factor(self, dens, momx, momy, momz, enth):
        cdef np.float64_t[:] rho = dens.ravel()

        cdef np.float64_t[:] mx = momx.ravel()
        cdef np.float64_t[:] my = momy.ravel()
        cdef np.float64_t[:] mz = momz.ravel()
        cdef np.float64_t[:] h = enth.ravel()

        out = np.empty_like(dens)
        cdef np.float64_t[:] outp = out.ravel()
        cdef int i

        for i in range(outp.shape[0]):
            outp[i] = self._lorentz_factor(rho[i], mx[i], my[i], mz[i], h[i])
        return out

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def sound_speed(self, temp, enth):
        cdef np.float64_t[:] kT = temp.ravel()
        cdef np.float64_t[:] h = enth.ravel()
        out = np.empty_like(kT)
        cdef np.float64_t[:] outp = out.ravel()

        cdef int i
        for i in range(outp.shape[0]):
            outp[i] = self.cs(kT[i], h[i], self._gamma)
        return out

    cdef np.float64_t _four_vel(
            self,
            np.float64_t rho,
            np.float64_t mi,
            np.float64_t h,
        ) noexcept nogil:
        return mi / (rho * (h + 1.0))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def four_velocity_xyz(self, dens, mom, enth):
        cdef np.float64_t[:] rho = dens.ravel()
        cdef np.float64_t[:] mi = mom.ravel()
        cdef np.float64_t[:] h = enth.ravel()

        out = np.empty_like(dens)
        cdef np.float64_t[:] outp = out.ravel()
        cdef int i

        for i in range(outp.shape[0]):
            outp[i] = self._four_vel(rho[i], mi[i], h[i])
        return out

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def kinetic_energy_density(self, dens, momx, momy, momz, temp, enth):
        cdef np.float64_t[:] rho = dens.ravel()
        cdef np.float64_t[:] mx = momx.ravel()
        cdef np.float64_t[:] my = momy.ravel()
        cdef np.float64_t[:] mz = momz.ravel()
        cdef np.float64_t[:] kT = temp.ravel()
        cdef np.float64_t[:] h = enth.ravel()

        out = np.empty_like(dens)
        cdef np.float64_t[:] outp = out.ravel()
        cdef np.float64_t lf, u2, ux, uy, uz
        cdef int i

        for i in range(outp.shape[0]):
            ux = self._four_vel(rho[i], mx[i], h[i])
            uy = self._four_vel(rho[i], my[i], h[i])
            uz = self._four_vel(rho[i], mz[i], h[i])
            u2 = ux**2 + uy**2 + uz**2
            lf = self._lorentz_factor(rho[i], mx[i], my[i], mz[i], h[i])
            gm1 = u2 / (lf + 1.0)
            p = rho[i] / lf * kT[i]
            outp[i] = gm1 * (rho[i] * (h[i] + 1.0) + p)
        return out

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def mach_number(self, dens, momx, momy, momz, temp, enth):
        cdef np.float64_t[:] rho = dens.ravel()
        cdef np.float64_t[:] mx = momx.ravel()
        cdef np.float64_t[:] my = momy.ravel()
        cdef np.float64_t[:] mz = momz.ravel()
        cdef np.float64_t[:] kT = temp.ravel()
        cdef np.float64_t[:] h = enth.ravel()

        out = np.empty_like(dens)
        cdef np.float64_t[:] outp = out.ravel()
        cdef np.float64_t cs, us, u
        cdef int i

        for i in range(outp.shape[0]):
            cs = self.cs(kT[i], h[i], self._gamma)
            us = cs / math.sqrt(1.0 - cs**2)
            u = math.sqrt(mx[i]**2 + my[i]**2 + mz[i]**2) / (rho[i] * (h[i] + 1.0))
            outp[i] = u / us
        return out
