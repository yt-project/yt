from yt.mods import *
import matplotlib.cm as cm

# Load the dataset.
pf = load(os.path.join(ytcfg.get("yt", "test_data_dir"), "IsolatedGalaxy/galaxy0030/galaxy0030"))

# Create projections using each colormap available.
p = ProjectionPlot(pf, "z", "density", weight_field = "density", width=0.4)

for cmap in cm.datad:
    if cmap.startswith("idl"):
        continue
    p.set_cmap(field="density", cmap=cmap)
    p.annotate_title(cmap)
    p.save('Projection_%s.png' % cmap.replace(' ', '_'))
