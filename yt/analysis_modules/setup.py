#!/usr/bin/env python
import setuptools


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('analysis_modules', parent_package, top_path)
    config.make_config_py()  # installs __config__.py
    #config.make_svn_version_py()
    config.add_subpackage("absorption_spectrum")
    config.add_subpackage("coordinate_transformation")
    config.add_subpackage("cosmology")
    config.add_subpackage("halo_finding")
    config.add_subpackage("halo_mass_function")
    config.add_subpackage("halo_merger_tree")
    config.add_subpackage("halo_profiler")
    config.add_subpackage("hierarchy_subset")
    config.add_subpackage("level_sets")
    config.add_subpackage("radial_column_density")
    config.add_subpackage("spectral_integrator")
    config.add_subpackage("star_analysis")
    config.add_subpackage("two_point_functions")
    return config
