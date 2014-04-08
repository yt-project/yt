from yt.mods import *

# For this example we will use h5py to write to our output file.
import h5py

pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")

level = 2
dims = pf.domain_dimensions * pf.refine_by**level

# Now, we construct an object that describes the data region and structure we
# want
cube = pf.covering_grid(2, # The level we are willing to extract to; higher
                             # levels than this will not contribute to the data!
                          left_edge=[0.0, 0.0, 0.0], 
                          # And any fields to preload (this is optional!)
                          dims = dims,
                          fields=["density"])

# Now we open our output file using h5py
# Note that we open with 'w' which will overwrite existing files!
f = h5py.File("my_data.h5", "w") 

# We create a dataset at the root note, calling it density...
f.create_dataset("/density", data=cube["density"])

# We close our file
f.close()
