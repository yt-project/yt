import yt
import numpy as np
from yt.visualization.volume_rendering.api import BoxSource, CoordinateVectorSource

# Load the dataset.
ds = yt.load("Enzo_64/DD0043/data0043")
im, sc = yt.volume_render(ds, ('gas','density'))
sc.get_source(0).transfer_function.grey_opacity=True

sc.annotate_domain(ds)
im = sc.render()
im.write_png("%s_vr_domain.png" % ds)

sc.annotate_grids(ds)
im = sc.render()
im.write_png("%s_vr_grids.png" % ds)

# Here we can draw the coordinate vectors on top of the image by processing
# it through the camera. Then save it out.
sc.annotate_axes()
im = sc.render()
im.write_png("%s_vr_coords.png" % ds)
