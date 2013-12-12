#
# Physical Constants and Units Conversion Factors
#
# Values for these constants, unless otherwise noted, are drawn from IAU,
# IUPAC, and NIST data, whichever is newer.
# http://maia.usno.navy.mil/NSFA/IAU2009_consts.html
# http://goldbook.iupac.org/list_goldbook_phys_constants_defs.html
# http://physics.nist.gov/cuu/Constants/index.html

# Elementary masses
mass_electron_grams = 9.10938291e-28
amu_grams         = 1.660538921e-24

# Solar values (see Mamajek 2012)
# https://sites.google.com/site/mamajeksstarnotes/bc-scale
mass_sun_grams = 1.98841586e33
temp_sun_kelvin = 5870.0
luminosity_sun_ergs_per_sec = 3.8270e33
metallicity_sun = 0.02041

# Conversion Factors:  X au * mpc_per_au = Y mpc
# length
mpc_per_mpc   = 1e0
mpc_per_kpc   = 1e-3
mpc_per_pc    = 1e-6
mpc_per_au    = 4.84813682e-12
mpc_per_rsun  = 2.253962e-14
mpc_per_miles = 5.21552871e-20
mpc_per_km    = 3.24077929e-20
mpc_per_cm    = 3.24077929e-25
kpc_per_cm    = mpc_per_cm / mpc_per_kpc
pc_per_cm     = mpc_per_cm / mpc_per_pc
km_per_pc     = 3.08567758e13
km_per_m      = 1e-3
km_per_cm     = 1e-5
ly_per_cm     = 1.05702341e-18
rsun_per_cm   = 1.4378145e-11
au_per_cm     = 6.68458712e-14

m_per_fpc     = 0.0324077929

kpc_per_mpc   = 1.0 / mpc_per_kpc
pc_per_mpc    = 1.0 / mpc_per_pc
au_per_mpc    = 1.0 / mpc_per_au
rsun_per_mpc  = 1.0 / mpc_per_rsun
miles_per_mpc = 1.0 / mpc_per_miles
km_per_mpc    = 1.0 / mpc_per_km
cm_per_mpc    = 1.0 / mpc_per_cm
cm_per_kpc    = 1.0 / kpc_per_cm
cm_per_km     = 1.0 / km_per_cm
pc_per_km     = 1.0 / km_per_pc
cm_per_pc     = 1.0 / pc_per_cm
cm_per_ly     = 1.0 / ly_per_cm
cm_per_rsun   = 1.0 / rsun_per_cm
cm_per_au     = 1.0 / au_per_cm

# time
# "IAU Style Manual" by G.A. Wilkins, Comm. 5, in IAU Transactions XXB (1989)
sec_per_Gyr  = 31.5576e15
sec_per_Myr  = 31.5576e12
sec_per_kyr  = 31.5576e9
sec_per_year = 31.5576e6
sec_per_day  = 86400.0
sec_per_hr   = 3600.0
sec_per_min  = 60.0
day_per_year = 365.25

# temperature / energy
boltzmann_constant_erg_per_K = 1.3806488e-16
erg_per_eV = 1.602176562e-12
erg_per_keV = erg_per_eV * 1.0e3
K_per_keV = erg_per_keV / boltzmann_constant_erg_per_K
keV_per_K = 1.0 / K_per_keV

# Cosmological constants
rho_crit_g_cm3_h2 = 1.8788e-29
hubble_constant_hertz = 2.19724836e-18 # Planck 2013

# Misc. Approximations
mass_mean_atomic_cosmology = 1.22
mass_mean_atomic_galactic = 2.3

# Miscellaneous
HUGE = 1.0e90
TINY = 1.0e-40
