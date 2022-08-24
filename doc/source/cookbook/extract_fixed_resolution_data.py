# For this example we will use h5py to write to our output file.
import h5py

import yt

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

level = 2
dims = ds.domain_dimensions * ds.refine_by**level

# We construct an object that describes the data region and structure we want
# In this case, we want all data up to the maximum "level" of refinement
# across the entire simulation volume.  Higher levels than this will not
# contribute to our covering grid.
cube = ds.covering_grid(
    level,
    left_edge=[0.0, 0.0, 0.0],
    dims=dims,
    # And any fields to preload (this is optional!)
    fields=[("gas", "density")],
)

# Now we open our output file using h5py
# Note that we open with 'w' (write), which will overwrite existing files!
f = h5py.File("my_data.h5", mode="w")

# We create a dataset at the root, calling it "density"
f.create_dataset("/density", data=cube[("gas", "density")])

# We close our file
f.close()

# If we want to then access this datacube in the h5 file, we can now...
f = h5py.File("my_data.h5", mode="r")
print(f["density"][()])
