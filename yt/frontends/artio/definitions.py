yt_to_art = {
    "Density": "HVAR_GAS_DENSITY",
    "TotalEnergy": "HVAR_GAS_ENERGY",
    "GasEnergy": "HVAR_INTERNAL_ENERGY",
    "Pressure": "HVAR_PRESSURE",
    "XMomentumDensity": "HVAR_MOMENTUM_X",
    "YMomentumDensity": "HVAR_MOMENTUM_Y",
    "ZMomentumDensity": "HVAR_MOMENTUM_Z",
    "Gamma": "HVAR_GAMMA",
    "MetalDensitySNIa": "HVAR_METAL_DENSITY_Ia",
    "MetalDensitySNII": "HVAR_METAL_DENSITY_II",
    "Potential": "VAR_POTENTIAL",
    "PotentialHydro": "VAR_POTENTIAL_HYDRO",
    "particle_position_x": "POSITION_X",
    "particle_position_y": "POSITION_Y",
    "particle_position_z": "POSITION_Z",
    "particle_velocity_x": "VELOCITY_X",
    "particle_velocity_y": "VELOCITY_Y",
    "particle_velocity_z": "VELOCITY_Z",
    "particle_mass": "MASS",
    "particle_index": "PID",
    "particle_species": "SPECIES",
    "creation_time": "BIRTH_TIME",
    "particle_mass_initial": "INITIAL_MASS",
    "particle_metallicity1": "METALLICITY_SNIa",
    "particle_metallicity2": "METALLICITY_SNII",
    "stars": "STAR",
    "nbody": "N-BODY",
}
art_to_yt = dict(zip(yt_to_art.values(), yt_to_art.keys()))


class ARTIOconstants:
    def __init__(self):
        self.yr = 365.25 * 86400
        self.Myr = 1.0e6 * self.yr
        self.Gyr = 1.0e9 * self.yr

        self.pc = 3.0856775813e18
        self.kpc = 1.0e3 * self.pc
        self.Mpc = 1.0e6 * self.pc

        self.kms = 1.0e5

        self.mp = 1.672621637e-24
        self.k = 1.3806504e-16
        self.G = 6.67428e-8
        self.c = 2.99792458e10

        self.eV = 1.602176487e-12
        self.amu = 1.660538782e-24
        self.mH = 1.007825 * self.amu
        self.mHe = 4.002602 * self.amu

        self.Msun = 1.32712440018e26 / self.G
        self.Zsun = 0.0199

        self.Yp = 0.24
        self.wmu = 4.0 / (8.0 - 5.0 * self.Yp)
        self.wmu_e = 1.0 / (1.0 - 0.5 * self.Yp)
        self.XH = 1.0 - self.Yp
        self.XHe = 0.25 * self.Yp
        self.gamma = 5.0 / 3.0

        self.sigmaT = 6.6524e-25
