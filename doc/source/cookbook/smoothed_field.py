import yt

# Load a Gadget dataset following the demonstration notebook.
fname = 'GadgetDiskGalaxy/snapshot_200.hdf5'

unit_base = {'UnitLength_in_cm'         : 3.08568e+21,
             'UnitMass_in_g'            :   1.989e+43,
             'UnitVelocity_in_cm_per_s' :      100000}

bbox_lim = 1e5  # kpc

bbox = [[-bbox_lim, bbox_lim],
        [-bbox_lim, bbox_lim],
        [-bbox_lim, bbox_lim]]

ds = yt.load(fname, unit_base=unit_base, bounding_box=bbox)

# Create a derived field, the metal density.
def _metal_density(field, data):
    density = data['PartType0', 'Density']
    Z = data['PartType0', 'metallicity']
    return density * Z

# Add it to the dataset.
ds.add_field(('PartType0', 'metal_density'), function=_metal_density,
             units="g/cm**3", particle_type=True)


# Add the corresponding smoothed field to the dataset.
from yt.fields.particle_fields import add_volume_weighted_smoothed_field

add_volume_weighted_smoothed_field('PartType0', 'Coordinates', 'Masses',
                                   'SmoothingLength', 'Density',
                                   'metal_density', ds.field_info)

# Define the region where the disk galaxy is. (See the Gadget notebook for
# details. Here I make the box a little larger than needed to eliminate the
# margin effect.)
center = ds.arr([31996, 31474, 28970], "code_length")
box_size = ds.quan(250, "code_length")
left_edge = center - box_size/2*1.1
right_edge = center + box_size/2*1.1
box = ds.box(left_edge=left_edge, right_edge=right_edge)

# And make a projection plot!
yt.ProjectionPlot(ds, 'z',
                  ('deposit', 'PartType0_smoothed_metal_density'),
                  center=center, width=box_size, data_source=box).save()
