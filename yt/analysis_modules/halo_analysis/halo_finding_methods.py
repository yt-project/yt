"""
Halo Finding methods



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np

from yt.analysis_modules.halo_finding.halo_objects import \
    FOFHaloFinder, HOPHaloFinder
from yt.frontends.halo_catalogs.halo_catalog.data_structures import \
    HaloCatalogDataset
from yt.frontends.stream.data_structures import \
    load_particles

from .operator_registry import \
    finding_method_registry


def add_finding_method(name, function):
    finding_method_registry[name] = HaloFindingMethod(function)
    
class HaloFindingMethod(object):
    r"""
    A halo finding method is a callback that performs halo finding on a 
    dataset and returns a new dataset that is the loaded halo finder output.
    """
    def __init__(self, function, args=None, kwargs=None):
        self.function = function
        self.args = args
        if self.args is None: self.args = []
        self.kwargs = kwargs
        if self.kwargs is None: self.kwargs = {}

    def __call__(self, ds):
        return self.function(ds, *self.args, **self.kwargs)

def _hop_method(pf, **finder_kwargs):
    r"""
    Run the Hop halo finding method.
    """
    
    halo_list = HOPHaloFinder(pf, **finder_kwargs)
    halos_pf = _parse_old_halo_list(pf, halo_list)
    return halos_pf
add_finding_method("hop", _hop_method)

def _fof_method(pf, **finder_kwargs):
    r"""
    Run the FoF halo finding method.
    """

    halo_list = FOFHaloFinder(pf, **finder_kwargs)
    halos_pf = _parse_old_halo_list(pf, halo_list)
    return halos_pf
add_finding_method("fof", _fof_method)

def _rockstar_method(pf, **finder_kwargs):
    r"""
    Run the Rockstar halo finding method.
    """

    from yt.frontends.halo_catalogs.rockstar.data_structures import \
     RockstarDataset
    from yt.analysis_modules.halo_finding.rockstar.api import \
     RockstarHaloFinder
    
    rh = RockstarHaloFinder(pf, **finder_kwargs)
    rh.run()


    halos_pf = RockstarDataset("rockstar_halos/halos_0.0.bin")
    try:
        halos_pf.create_field_info()
    except ValueError:
        return None

    return halos_pf
add_finding_method("rockstar", _rockstar_method)

def _parse_old_halo_list(data_pf, halo_list):
    r"""
    Convert the halo list into a loaded dataset.
    """

    num_halos = len(halo_list)

    if num_halos == 0: return None

    # Set up fields that we want to pull from identified halos and their units
    new_fields = ['particle_identifier', 'particle_mass', 'particle_position_x', 
        'particle_position_y','particle_position_z',
        'virial_radius']
    new_units = [ '', 'g', 'cm', 'cm','cm','cm']

    # Set up a dictionary based on those fields 
    # with empty arrays where we will fill in their values
    halo_properties = { f : (np.zeros(num_halos),unit) \
        for f, unit in zip(new_fields,new_units)}

    # Iterate through the halos pulling out their positions and virial quantities
    # and filling in the properties dictionary
    for i,halo in enumerate(halo_list):
        halo_properties['particle_identifier'][0][i] = i
        halo_properties['particle_mass'][0][i] = halo.virial_mass().in_cgs()
        halo_properties['virial_radius'][0][i] = halo.virial_radius().in_cgs()

        com = halo.center_of_mass().in_cgs()
        halo_properties['particle_position_x'][0][i] = com[0]
        halo_properties['particle_position_y'][0][i] = com[1]
        halo_properties['particle_position_z'][0][i] = com[2]

    # Define a bounding box based on original data pf
    bbox = np.array([data_pf.domain_left_edge.in_cgs(),
            data_pf.domain_right_edge.in_cgs()]).T

    # Create a pf with the halos as particles
    particle_pf = load_particles(halo_properties, 
            bbox=bbox, length_unit = 1, mass_unit=1)

    # Create the field info dictionary so we can reference those fields
    particle_pf.create_field_info()

    for attr in ["current_redshift", "current_time",
                 "domain_dimensions",
                 "cosmological_simulation", "omega_lambda",
                 "omega_matter", "hubble_constant"]:
        attr_val = getattr(data_pf, attr)
        setattr(particle_pf, attr, attr_val)
    particle_pf.current_time = particle_pf.current_time.in_cgs()
    
    return particle_pf
