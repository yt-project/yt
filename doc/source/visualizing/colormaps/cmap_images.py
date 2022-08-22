import matplotlib as mpl

import yt

# Load the dataset.
ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

# Create projections using each colormap available.
p = yt.ProjectionPlot(ds, "z", "density", weight_field="density", width=0.4)

for cmap in mpl.colormaps:
    if cmap.startswith("idl"):
        continue
    p.set_cmap(field="density", cmap=cmap)
    p.annotate_title(cmap)
    p.save(f"Projection_{cmap.replace(' ', '_')}.png")
