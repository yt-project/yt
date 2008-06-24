from yt.config import ytcfg

ytcfg["lagos","serialize"] = 'False'

from yt.mods import *

oso = lagos.OrionStaticOutput("tests/uniformCollapse_base64_maxLev4/plt0005")
pc = raven.PlotCollection(oso,center=3*[0.0])
ax = 0
pc.add_slice("z-velocity",ax)
for p in pc.plots:
    p.set_xlim(p.data['px'].min(), p.data['px'].max())
    p.set_ylim(p.data['py'].min(), p.data['py'].max())
    p.set_zlim(p.data['Density'].min(), p.data['Density'].max())
    p.add_callback(raven.GridBoundaryCallback())
    p.set_log_field(False)    
    v1 = "%s-velocity" % (lagos.axis_names[lagos.x_dict[ax]])
    v2 = "%s-velocity" % (lagos.axis_names[lagos.y_dict[ax]])
    print v1,v2
    p.add_callback(raven.QuiverCallback(v1,v2,32))


pc.save('vecTest')
