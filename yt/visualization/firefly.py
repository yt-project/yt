"""
wrapper for firefly_api to output dataset to Firefly


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import numpy as np 
import os

from firefly_api.reader import Reader
from firefly_api.particlegroup import ParticleGroup
from firefly_api.errors import FireflyError,FireflyWarning,warnings

def create_firefly_object(
    region,
    path_to_firefly,
    fields_to_include= [],
    fields_units = [],
    default_decimation_factor = 100,
    velocity_units = 'km/s',
    coordinate_units = 'kpc',
    show_unused_fields=0):
    r"""This function links a region of data stored in a yt dataset
        to the Python frontend API for Firefly, a browser-based 
        particle visualization platform. 

        Example usage: 

            ramses_ds = yt.load(
                "/Users/agurvich/Desktop/yt_workshop/"+
                "DICEGalaxyDisk_nonCosmological/output_00002/info_00002.txt")

            region = ramses_ds.sphere(ramses_ds.domain_center,(1000,'kpc'))

            reader = create_firefly_object(
                region,
                path_to_firefly="/Users/agurvich/research/repos/Firefly",
                fields_to_include=[
                    'particle_extra_field_1',
                    'particle_extra_field_2'],
                fields_units = ['dimensionless','dimensionless'])

            reader.options['color']['io']=[1,1,0,1]
            reader.particleGroups[0].decimation_factor=100
            reader.dumpToJSON()

        Parameters
        ----------
        region : YTRegion object
            A region that includes particle data.

        path_to_firefly : string
            The (ideally) absolute path to the direction containing the index.html
            file of Firefly. 

        fields_to_include : array_like of strings
            A list of fields that you want to include in your 
            Firefly visualization for on-the-fly filtering and
            colormapping. 
        
        default_decimation_factor : integer
            The factor by which you want to decimate each particle group
            by (e.g. if there are 1e7 total particles in your simulation
            you might want to set this to 100 at first). Randomly samples
            your data like `shuffled_data[::decimation_factor]` so as to 
            not overtax a system. This is adjustable on a per particle group
            basis by changing the returned reader's 
            `reader.particleGroup[i].decimation_factor` before calling 
            `reader.dumpToJSON()`.
        
        velocity_units : string
            The units that the velocity should be converted to in order to 
            show streamlines in Firefly. Defaults to km/s. 
        
        coordinate_units: string
            The units that the coordinates should be converted to. Defaults to 
            kpc. 
        
        show_unused_fields: boolean
            A flag to optionally print the fields that are available, in the 
            dataset but were not explicitly requested to be tracked.

        Returns
        -------
        reader : firefly_api.reader.Reader object
            A reader object from the firefly_api, configured 
            to output
        """ 
    
    try:
        assert len(fields_units) == len(fields_to_include)
    except AssertionError:
        raise FireflyError("Each requested field must have units.")
    
    ## for safety, in case someone passes a float just cast it
    default_decimation_factor = int(default_decimation_factor)
    
    ## initialize a firefly reader instance
    reader = Reader(
        JSONdir=os.path.join(path_to_firefly,'data','yt'),
        write_startup=1,
        prefix='ytData',
        clean_JSONdir=True,
        )

    ## create a ParticleGroup object that contains *every* field
    for ptype in sorted(region.ds.particle_types_raw):
        ## skip this particle type if it has no particles in this region
        if region[ptype,'relative_particle_position'].shape[0]==0:
            continue
        
        if show_unused_fields:
            ## read the available extra fields from yt
            this_ptype_fields = region.ds.particle_fields_by_type[ptype]

            ## load the extra fields
            for field in this_ptype_fields:
                if 'position' in field or 'velocity' in field:
                    continue
                if field not in fields_to_include:
                    warnings.warn(FireflyWarning(
                        'detected (but did not request) {} {}'.format(ptype,field)))
                    continue

        ## you must have velocities (and they must be named "Velocities")
        tracked_arrays = [
            region[ptype,'relative_particle_velocity'].convert_to_units(velocity_units)]
        tracked_names = ['Velocities']

        ## explicitly go after the fields we want
        for field,units in zip(fields_to_include,fields_units):
            ## determine if you want to take the log of the field for Firefly
            log_flag = 'log(' in units 

            ## read the field array from the dataset
            this_field_array = region[ptype,field]

            ## fix the units string and prepend 'log' to the field for
            ##  the UI name
            if log_flag:
                units = units[len('log('):-1]
                field = 'log{}'.format(field)

            ## perform the unit conversion and take the log if 
            ##  necessary.
            this_field_array.convert_to_units(units)
            if log_flag:
                this_field_array = np.log10(this_field_array)

            ## add this array to the tracked arrays
            tracked_arrays +=[this_field_array]
            tracked_names = np.append(tracked_names,[field],axis=0)
        
        ## flag whether we want to filter and/or color by these fields
        ##  we'll assume yes for both cases, this can be changed after
        ##  the reader object is returned to the user.
        tracked_filter_flags = np.ones(len(tracked_names))
        tracked_colormap_flags = np.ones(len(tracked_names))

        ## create a firefly ParticleGroup for this particle type
        pg = ParticleGroup(
            UIname =  ptype,
            coordinates=region[ptype,'relative_particle_position'].convert_to_units(coordinate_units),
            tracked_arrays=tracked_arrays,
            tracked_names=tracked_names,
            tracked_filter_flags=tracked_filter_flags,
            tracked_colormap_flags=tracked_colormap_flags,
            decimation_factor=default_decimation_factor)
        
        ## bind this particle group to the firefly reader object
        reader.addParticleGroup(pg)

    return reader
