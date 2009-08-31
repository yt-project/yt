### A sample script to run the HaloProfiler.  ###
# See extensions/halo_profiler.html in the docs for full documentation.

import yt.extensions.HaloProfiler as HP

# Instantiate HaloProfiler for this dataset.
hp = HP.HaloProfiler("/Users/britton/EnzoRuns/cool_core_unreasonable/DD0242/DD0242")

# Add a filter to remove halos that have no profile points with overdensity 
# above 200, and with virial masses less than 1e14 solar masses.
# Also, return the virial mass and radius to be written out to a file.
hp.add_halo_filter(HP.VirialFilter,must_be_virialized=True,
                overdensity_field='ActualOverdensity',
                virial_overdensity=200,
                virial_filters=[['TotalMassMsun','>=','1e14']],
                virial_quantities=['TotalMassMsun','RadiusMpc'])

# Add profile fields.
hp.add_profile('CellVolume',weight_field=None,accumulation=True)
hp.add_profile('TotalMassMsun',weight_field=None,accumulation=True)
hp.add_profile('Density',weight_field=None,accumulation=False)
hp.add_profile('Temperature',weight_field='CellMassMsun',accumulation=False)
# Make profiles and output filtered halo list to FilteredQuantities.out.
hp.make_profiles(filename="FilteredQuantities.out")

# Add projection fields.
hp.add_projection('Density',weight_field=None)
hp.add_projection('Temperature',weight_field='Density')
hp.add_projection('Metallicity',weight_field='Density')
# Make projections for all three axes using the filtered halo list and 
# save data to hdf5 files.
hp.make_projections(save_cube=True,save_images=True,
                    halo_list='filtered',axes=[0,1,2])
