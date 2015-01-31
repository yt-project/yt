#!/usr/bin/env python

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('analysis_modules', parent_package, top_path)
    config.make_config_py()  # installs __config__.py
    config.add_subpackage("absorption_spectrum")
    config.add_subpackage("cosmological_observation")
    config.add_subpackage("halo_analysis")
    config.add_subpackage("halo_finding")
    config.add_subpackage("halo_mass_function")
    config.add_subpackage("level_sets")
    config.add_subpackage("particle_trajectories")
    config.add_subpackage("photon_simulator")
    config.add_subpackage("spectral_integrator")
    config.add_subpackage("star_analysis")
    config.add_subpackage("two_point_functions")
    config.add_subpackage("radmc3d_export")
    config.add_subpackage("sunrise_export")
    config.add_subpackage("sunyaev_zeldovich")
    config.add_subpackage("particle_trajectories")
    config.add_subpackage("photon_simulator")
    config.add_subpackage("ppv_cube")
    return config
