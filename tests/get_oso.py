from yt.config import ytcfg

ytcfg["lagos","serialize"] = 'False'

from yt.mods import *

oso = lagos.OrionStaticOutput("tests/uniformCollapse_base64_maxLev4/plt0005")
oso.conversion_factors["Density"] = 1.0
pc = raven.PlotCollection(oso,center=3*[0.0])
pc.add_projection("Density",0)
pc.add_projection("Density",1)
pc.add_projection("Density",2)
#pc.add_projection("Density",0)
#pc.add_projection("Density",1)
#pc.add_projection("Density",2)
pc.set_width((oso["DomainRightEdge"]-oso["DomainLeftEdge"])[0],'1')
for p in pc.plots:
    #p.center = oso["DomainRightEdge"]
    #p.data['px'] -= oso["DomainLeftEdge"][0]
    #p.data['py'] -= oso["DomainLeftEdge"][0]
    p.data["Density"] /= (oso["DomainRightEdge"]-oso["DomainLeftEdge"])[0]
    p.data["Density"] /= (oso["cm"])
    p.set_xlim(p.data['px'].min(), p.data['px'].max())
    p.set_ylim(p.data['py'].min(), p.data['py'].max())
    p.set_zlim(p.data['Density'].min(), p.data['Density'].max())
    x1,x2 = p._axes.get_xlim()
    y1,y2 = p._axes.get_ylim()
    #p.add_callback(raven.LinePlotCallback([(x1+x2)/2.0]*2, [y1,y2]))
    p.add_callback(raven.GridBoundaryCallback())
    p.set_log_field(False)
pc.save('hi')
#for p in pc.plots: p._axes.plot([0.0,0.0],p.ylim)

