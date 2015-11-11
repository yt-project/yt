import yt
import numpy as np
from yt.visualization.volume_rendering.api import BoxSource, CoordinateVectorSource

# Load the dataset.
ds = yt.load("Enzo_64/DD0043/data0043")
sc = yt.create_scene(ds, ('gas','density'))
sc.get_source(0).transfer_function.grey_opacity=True

sc.annotate_domain(ds)
sc.render()
sc.save("%s_vr_domain.png" % ds)

sc.annotate_grids(ds)
sc.render()
sc.save("%s_vr_grids.png" % ds)

# Here we can draw the coordinate vectors on top of the image by processing
# it through the camera. Then save it out.
sc.annotate_axes()
sc.render()
sc.save("%s_vr_coords.png" % ds)
