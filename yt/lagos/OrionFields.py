

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

def _generate_translation(mine, theirs):
    add_field(theirs, function=lambda a, b: b[mine], take_log=True)

for f,v in translation_dict.items():
    if v not in OrionFieldInfo:
        add_field(v, function=lambda a,b: None, take_log=True,
                  validators = [ValidateDataField(v)])
    #print "Setting up translator from %s to %s" % (v, f)
    _generate_translation(v, f)
