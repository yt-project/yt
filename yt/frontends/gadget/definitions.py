gadget_header_specs = {
    "default": (
        ("Npart", 6, "i"),
        ("Massarr", 6, "d"),
        ("Time", 1, "d"),
        ("Redshift", 1, "d"),
        ("FlagSfr", 1, "i"),
        ("FlagFeedback", 1, "i"),
        ("Nall", 6, "i"),
        ("FlagCooling", 1, "i"),
        ("NumFiles", 1, "i"),
        ("BoxSize", 1, "d"),
        ("Omega0", 1, "d"),
        ("OmegaLambda", 1, "d"),
        ("HubbleParam", 1, "d"),
        ("FlagAge", 1, "i"),
        ("FlagMetals", 1, "i"),
        ("NallHW", 6, "i"),
        ("unused", 16, "i"),
    ),
    "pad32": (("empty", 32, "c"),),
    "pad64": (("empty", 64, "c"),),
    "pad128": (("empty", 128, "c"),),
    "pad256": (("empty", 256, "c"),),
}

gadget_ptype_specs = {"default": ("Gas", "Halo", "Disk", "Bulge", "Stars", "Bndry")}

gadget_field_specs = {
    "default": (
        "Coordinates",
        "Velocities",
        "ParticleIDs",
        "Mass",
        ("InternalEnergy", "Gas"),
        ("Density", "Gas"),
        ("SmoothingLength", "Gas"),
    ),
    "agora_unlv": (
        "Coordinates",
        "Velocities",
        "ParticleIDs",
        "Mass",
        ("InternalEnergy", "Gas"),
        ("Density", "Gas"),
        ("Electron_Number_Density", "Gas"),
        ("HI_NumberDensity", "Gas"),
        ("SmoothingLength", "Gas"),
    ),
    "group0000": (
        "Coordinates",
        "Velocities",
        "ParticleIDs",
        "Mass",
        "Potential",
        ("Temperature", "Gas"),
        ("Density", "Gas"),
        ("ElectronNumberDensity", "Gas"),
        ("HI_NumberDensity", "Gas"),
        ("SmoothingLength", "Gas"),
        ("StarFormationRate", "Gas"),
        ("DelayTime", "Gas"),
        ("FourMetalFractions", ("Gas", "Stars")),
        ("MaxTemperature", ("Gas", "Stars")),
        ("NStarsSpawned", ("Gas", "Stars")),
        ("StellarAge", "Stars"),
    ),
    "magneticum_box2_hr": (
        "Coordinates",
        "Velocities",
        "ParticleIDs",
        "Mass",
        ("InternalEnergy", "Gas"),
        ("Density", "Gas"),
        ("SmoothingLength", "Gas"),
        ("ColdFraction", "Gas"),
        ("Temperature", "Gas"),
        ("StellarAge", ("Stars", "Bndry")),
        "Potential",
        ("InitialMass", "Stars"),
        ("ElevenMetalMasses", ("Gas", "Stars")),
        ("StarFormationRate", "Gas"),
        ("TrueMass", "Bndry"),
        ("AccretionRate", "Bndry"),
    ),
}

gadget_hdf5_ptypes = (
    "PartType0",
    "PartType1",
    "PartType2",
    "PartType3",
    "PartType4",
    "PartType5",
)

SNAP_FORMAT_2_OFFSET = 16

"""
Here we have a dictionary of possible element species defined in Gadget
datasets, keyed by the number of elements. In some cases, these are mass
fractions, in others, they are metals--the context for the dataset will
determine this. The "Ej" key is for the total mass of all elements that
are not explicitly listed.
"""
elem_names_opts = {
    4: ["C", "O", "Si", "Fe"],
    7: ["C", "N", "O", "Mg", "Si", "Fe", "Ej"],
    8: ["He", "C", "O", "Mg", "S", "Si", "Fe", "Ej"],
    11: ["He", "C", "Ca", "O", "N", "Ne", "Mg", "S", "Si", "Fe", "Ej"],
    15: [
        "He",
        "C",
        "Ca",
        "O",
        "N",
        "Ne",
        "Mg",
        "S",
        "Si",
        "Fe",
        "Na",
        "Al",
        "Ar",
        "Ni",
        "Ej",
    ],
}
