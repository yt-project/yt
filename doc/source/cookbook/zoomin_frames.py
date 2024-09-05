import numpy as np

import yt

# load data
fn = "IsolatedGalaxy/galaxy0030/galaxy0030"
ds = yt.load(fn)

# This is the number of frames to make -- below, you can see how this is used.
n_frames = 5

# This is the minimum size in smallest_dx of our last frame.
# Usually it should be set to something like 400, but for THIS
# dataset, we actually don't have that great of resolution.
min_dx = 40

frame_template = "frame_%05i"  # Template for frame filenames

p = yt.SlicePlot(ds, "z", ("gas", "density"))  # Add our slice, along z
p.annotate_contour(("gas", "temperature"))  # We'll contour in temperature

# What we do now is a bit fun.  "enumerate" returns a tuple for every item --
# the index of the item, and the item itself.  This saves us having to write
# something like "i = 0" and then inside the loop "i += 1" for ever loop.  The
# argument to enumerate is the 'logspace' function, which takes a minimum and a
# maximum and the number of items to generate.  It returns 10^power of each
# item it generates.

for i, v in enumerate(
    np.logspace(0, np.log10(ds.index.get_smallest_dx() * min_dx), n_frames)
):
    # We set our width as necessary for this frame
    p.set_width(v, "unitary")
    # save
    p.save(frame_template % (i))
