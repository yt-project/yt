

from UniversalFields import *

add_field = add_orion_field

add_field("density", function=lambda a,b: None, take_log=True,
          validators = [ValidateDataField("density")],
          units=r"\rm{g}/\rm{cm}^3")

translation_dict = {"x-velocity": "xvel",
                    "y-velocity": "yvel",
                    "z-velocity": "zvel",
                    "Density": "density",
                    "Gas_Energy": "eden",
                    "Temperature": "temperature",
                    "x-momentum": "xmom",
                    "y-momentum": "ymom",
                    "z-momentum": "zmom"
                   }

add_field("Density", function=lambda a,b: b["density"])

#for f,v in translation_dict.items():
#    if v not in OrionFieldInfo:
#        add_field(v, function=lambda a,b: None, take_log=True,
#                  validators = [ValidateDataField(v)])
#    def func(a, b):
#        print "Returning %s (%s)" % (v, f)
#        return b[v]
#    print "Setting up translator from %s to %s" % (v, f)
#    add_field(f, function=func, take_log=True)
