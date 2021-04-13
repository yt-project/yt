vals = [
    "FluidQuantities",
    "Particles",
    "Parameters",
    "Units",
    "ReadOnDemand",
    "LoadRawData",
    "LevelOfSupport",
    "ContactPerson",
]


class CodeSupport:
    def __init__(self, **kwargs):
        self.support = {}
        for v in vals:
            self.support[v] = "N"
        for k, v in kwargs.items():
            if k in vals:
                self.support[k] = v


Y = "Y"
N = "N"

code_names = ["Enzo", "Orion", "FLASH", "RAMSES", "Chombo", "Gadget", "ART", "ZEUS"]


codes = dict(
    Enzo=CodeSupport(
        FluidQuantities=Y,
        Particles=Y,
        Parameters=Y,
        Units=Y,
        ReadOnDemand=Y,
        LoadRawData=Y,
        ContactPerson="Matt Turk",
        LevelOfSupport="Full",
    ),
    Orion=CodeSupport(
        FluidQuantities=Y,
        Particles=N,
        Parameters=Y,
        Units=Y,
        ReadOnDemand=Y,
        LoadRawData=Y,
        ContactPerson="Jeff Oishi",
        LevelOfSupport="Full",
    ),
    FLASH=CodeSupport(
        FluidQuantities=Y,
        Particles=N,
        Parameters=N,
        Units=Y,
        ReadOnDemand=Y,
        LoadRawData=Y,
        ContactPerson="John !ZuHone",
        LevelOfSupport="Partial",
    ),
    RAMSES=CodeSupport(
        FluidQuantities=Y,
        Particles=N,
        Parameters=N,
        Units=N,
        ReadOnDemand=Y,
        LoadRawData=Y,
        ContactPerson="Matt Turk",
        LevelOfSupport="Partial",
    ),
    Chombo=CodeSupport(
        FluidQuantities=Y,
        Particles=N,
        Parameters=N,
        Units=N,
        ReadOnDemand=Y,
        LoadRawData=Y,
        ContactPerson="Jeff Oishi",
        LevelOfSupport="Partial",
    ),
    Gadget=CodeSupport(
        FluidQuantities=N,
        Particles=Y,
        Parameters=Y,
        Units=Y,
        ReadOnDemand=N,
        LoadRawData=N,
        ContactPerson="Chris Moody",
        LevelOfSupport="Partial",
    ),
    ART=CodeSupport(
        FluidQuantities=N,
        Particles=N,
        Parameters=N,
        Units=N,
        ReadOnDemand=N,
        LoadRawData=N,
        ContactPerson="Matt Turk",
        LevelOfSupport="None",
    ),
    ZEUS=CodeSupport(
        FluidQuantities=N,
        Particles=N,
        Parameters=N,
        Units=N,
        ReadOnDemand=N,
        LoadRawData=N,
        ContactPerson="Matt Turk",
        LevelOfSupport="None",
    ),
)

print("|| . ||", end=" ")
for c in code_names:
    print(f"{c} || ", end=" ")
print()

for vn in vals:
    print(f"|| !{vn} ||", end=" ")
    for c in code_names:
        print(f"{codes[c].support[vn]} || ", end=" ")
    print()
