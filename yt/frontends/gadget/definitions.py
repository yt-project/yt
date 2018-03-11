"""
Gadget definitions




"""

#-----------------------------------------------------------------------------
# Copyright (c) 2014, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

gadget_header_specs = dict(
    default      = (('Npart', 6, 'i'),
                    ('Massarr', 6, 'd'),
                    ('Time', 1, 'd'),
                    ('Redshift', 1, 'd'),
                    ('FlagSfr', 1, 'i'),
                    ('FlagFeedback', 1, 'i'),
                    ('Nall', 6, 'i'),
                    ('FlagCooling', 1, 'i'),
                    ('NumFiles', 1, 'i'),
                    ('BoxSize', 1, 'd'),
                    ('Omega0', 1, 'd'),
                    ('OmegaLambda', 1, 'd'),
                    ('HubbleParam', 1, 'd'),
                    ('FlagAge', 1, 'i'),
                    ('FlagMetals', 1, 'i'),
                    ('NallHW', 6, 'i'),
                    ('unused', 16, 'i')),
    pad32       = (('empty',  32, 'c'),),
    pad64       = (('empty',  64, 'c'),),
    pad128      = (('empty', 128, 'c'),),
    pad256      = (('empty', 256, 'c'),),
)

gadget_ptype_specs = dict(
    default = ( "Gas",
                "Halo",
                "Disk",
                "Bulge",
                "Stars",
                "Bndry" )
)

gadget_field_specs = dict(
    default = ( "Coordinates",
                "Velocities",
                "ParticleIDs",
                "Mass",
                ("InternalEnergy", "Gas"),
                ("Density", "Gas"),
                ("SmoothingLength", "Gas"),
    ),
    agora_unlv = ( "Coordinates",
                   "Velocities",
                   "ParticleIDs",
                   "Mass",
                   ("InternalEnergy", "Gas"),
                   ("Density", "Gas"),
                   ("Electron_Number_Density", "Gas"),
                   ("HI_NumberDensity", "Gas"),
                   ("SmoothingLength", "Gas"),
    ),
    group0000 =  ( "Coordinates",
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
                   ("StellarAge", "Stars")
    ),
)

gadget_hdf5_ptypes  = (
    "PartType0",
    "PartType1",
    "PartType2",
    "PartType3",
    "PartType4",
    "PartType5"
)

SNAP_FORMAT_2_OFFSET = 16
