# distutils: include_dirs = LIB_DIR
# distutils: libraries = STD_LIBS
cimport cython
cimport libc.math as math
cimport numpy as np

import numpy as np


cdef np.float64_t h_eos4(np.float64_t kT, np.float64_t g) noexcept nogil:
    cdef np.float64_t x
    x = 2.25 * kT * kT
    return 2.5 * kT + x / (1.0 + math.sqrt(x + 1.0))

cdef np.float64_t h_eos(np.float64_t kT, np.float64_t g) noexcept nogil:
    return g * kT / (g - 1.0)

cdef np.float64_t gamma_eos(np.float64_t kT, np.float64_t g) noexcept nogil:
    return g

cdef np.float64_t gamma_eos4(np.float64_t kT, np.float64_t g) noexcept nogil:
    cdef np.float64_t x, c_p, c_v
    x = 2.25 * kT / math.sqrt(2.25 * kT * kT + 1.0)
    c_p = 2.5 + x
    c_v = 1.5 + x
    return c_p / c_v

cdef np.float64_t cs_eos4(np.float64_t kT, np.float64_t c, np.float64_t g) noexcept nogil:
    cdef np.float64_t hp, cs2
    hp = h_eos4(kT, 0.0) + 1.0
    cs2 = kT / (3.0 * hp)
    cs2 *= (5.0 * hp - 8.0 * kT) / (hp - kT)
    return c * math.sqrt(cs2)

cdef np.float64_t cs_eos(np.float64_t kT, np.float64_t c, np.float64_t g) noexcept nogil:
    cdef np.float64_t hp, cs2
    hp = h_eos(kT, g) + 1.0
    cs2 = g / hp * kT
    return c * math.sqrt(cs2)

ctypedef np.float64_t (*f2_type)(np.float64_t, np.float64_t) noexcept nogil
ctypedef np.float64_t (*f3_type)(np.float64_t, np.float64_t, np.float64_t) noexcept nogil

cdef class SRHDFields:
    cdef f2_type h
    cdef f2_type gamma
    cdef f3_type cs
    cdef np.float64_t _gamma, _c, _c2

    def __init__(self, int eos, np.float64_t gamma, np.float64_t clight):
        self._gamma = gamma
        self._c = clight
        self._c2 = clight * clight
        # Select aux functions based on eos no.
        if (eos == 4):
            self.h = h_eos4
            self.gamma = gamma_eos4
            self.cs = cs_eos4
        else:
            self.h = h_eos
            self.gamma = gamma_eos
            self.cs = cs_eos

    cdef np.float64_t _lorentz_factor(
            self,
            np.float64_t rho,
            np.float64_t mx,
            np.float64_t my,
            np.float64_t mz,
            np.float64_t kT,
        ) noexcept nogil:
        cdef np.float64_t u2
        cdef np.float64_t fac

        fac = (1.0 / (rho * (self.h(kT, self._gamma) + 1.0))) ** 2
        u2 = (mx * mx + my * my + mz * mz) * fac
        return math.sqrt(1.0 + u2 / self._c2)

    cdef np.float64_t _four_vel(
        self,
        np.float64_t rho,
        np.float64_t kT,
        np.float64_t mi
    ) noexcept nogil:
        return mi / (rho * (self.h(kT, self._gamma) + 1.0))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def lorentz_factor(self, dens, momx, momy, momz, temp):
        cdef np.float64_t[:] rho = dens.ravel()

        cdef np.float64_t[:] mx = momx.ravel()
        cdef np.float64_t[:] my = momy.ravel()
        cdef np.float64_t[:] mz = momz.ravel()
        cdef np.float64_t[:] kT = temp.ravel()

        out = np.empty_like(dens)
        cdef np.float64_t[:] outp = out.ravel()
        cdef int i

        for i in range(outp.shape[0]):
            outp[i] = self._lorentz_factor(rho[i], mx[i], my[i], mz[i], kT[i])
        return out

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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def sound_speed(self, temp):
        cdef np.float64_t[:] kT = temp.ravel()
        out = np.empty_like(kT)
        cdef np.float64_t[:] outp = out.ravel()

        cdef int i
        for i in range(outp.shape[0]):
            outp[i] = self.cs(kT[i], self._c, self._gamma)
        return out

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def four_velocity_xyz(self, mom, dens, temp):
        cdef np.float64_t[:] rho = dens.ravel()
        cdef np.float64_t[:] mi = mom.ravel()
        cdef np.float64_t[:] kT = temp.ravel()

        out = np.empty_like(dens)
        cdef np.float64_t[:] outp = out.ravel()
        cdef int i

        for i in range(outp.shape[0]):
            outp[i] = self._four_vel(rho[i], kT[i], mi[i])
        return out

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def velocity_xyz(self, dens, momx, momy, momz, temp, momi):
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
            lf = self._lorentz_factor(rho[i], mx[i], my[i], mz[i], kT[i])
            outp[i] = self._four_vel(rho[i], kT[i], mi[i]) / lf
        return out

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def density(self, dens, momx, momy, momz, temp):
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
            lf = self._lorentz_factor(rho[i], mx[i], my[i], mz[i], kT[i])
            outp[i] = rho[i] / lf
        return out

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def pressure(self, dens, momx, momy, momz, temp):
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
            lf = self._lorentz_factor(rho[i], mx[i], my[i], mz[i], kT[i])
            outp[i] = rho[i] / lf * self._c2 * kT[i]
        return out

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def specific_thermal_energy(self, dens, momx, momy, momz, temp):
        cdef np.float64_t[:] kT = temp.ravel()

        out = np.empty_like(dens)
        cdef np.float64_t[:] outp = out.ravel()
        cdef int i

        for i in range(outp.shape[0]):
            outp[i] = self._c2 * (self.h(kT[i], self._gamma)- kT[i])
        return out

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def kinetic_energy_density(self, dens, momx, momy, momz, temp):
        cdef np.float64_t[:] rho = dens.ravel()
        cdef np.float64_t[:] mx = momx.ravel()
        cdef np.float64_t[:] my = momy.ravel()
        cdef np.float64_t[:] mz = momz.ravel()
        cdef np.float64_t[:] kT = temp.ravel()

        out = np.empty_like(dens)
        cdef np.float64_t[:] outp = out.ravel()
        cdef np.float64_t lf, u2, hp, ux, uy, uz
        cdef int i

        for i in range(outp.shape[0]):
            hp = self.h(kT[i], self._gamma) + 1.0
            ux = self._four_vel(rho[i], kT[i], mx[i])
            uy = self._four_vel(rho[i], kT[i], my[i])
            uz = self._four_vel(rho[i], kT[i], mz[i])
            u2 = ux**2 + uy**2 + uz**2
            lf = self._lorentz_factor(rho[i], mx[i], my[i], mz[i], kT[i])
            gm1 = u2 / self._c2 / (lf + 1.0)
            p = rho[i] / lf * self._c2 * kT[i]
            outp[i] = gm1 * (rho[i] * hp * self._c2 + p)
        return out

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    def mach_number(self, dens, momx, momy, momz, temp):
        cdef np.float64_t[:] rho = dens.ravel()
        cdef np.float64_t[:] mx = momx.ravel()
        cdef np.float64_t[:] my = momy.ravel()
        cdef np.float64_t[:] mz = momz.ravel()
        cdef np.float64_t[:] kT = temp.ravel()

        out = np.empty_like(dens)
        cdef np.float64_t[:] outp = out.ravel()
        cdef np.float64_t cs
        cdef int i

        for i in range(outp.shape[0]):
            cs = self.cs(kT[i], self._c, self._gamma)
            us = cs / math.sqrt(1.0 - cs**2 / self._c2)
            u = math.sqrt(mx[i]**2 + my[i]**2 + mz[i]**2) / (rho[i] * (self.h(kT[i], self._gamma) + 1.0))
            outp[i] = u / us
        return out
