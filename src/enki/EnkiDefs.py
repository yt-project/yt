# We're going to be setting default vals here


MAX_SPHERE = 10             # Taken from Grid.h
MAX_DEPTH_OF_HIERARCHY = 50 # Taken from macros_and_parameters.h
DEFAULT_GHOST_ZONES    = 3  # Taken from macros_and_parameters.h


FieldTypes = {}

FieldTypes["Density"] = EnzoInterface.Density
FieldTypes["TotalEnergy"] = EnzoInterface.TotalEnergy
FieldTypes["GasEnergy"] = EnzoInterface.GasEnergy
FieldTypes["x-velocity"] = EnzoInterface.Velocity1
FieldTypes["y-velocity"] = EnzoInterface.Velocity2
FieldTypes["z-velocity"] = EnzoInterface.Velocity3
FieldTypes["colour"] = EnzoInterface.Metallicity
FieldTypes["Electron_Density"] = EnzoInterface.ElectronDensity
FieldTypes["HI_Density"] = EnzoInterface.HIDensity
FieldTypes["HII_Density"] = EnzoInterface.HIIDensity
FieldTypes["HeI_Density"] = EnzoInterface.HeIDensity
FieldTypes["HeII_Density"] = EnzoInterface.HeIIDensity
FieldTypes["HeIII_Density"] = EnzoInterface.HeIIIDensity
FieldTypes["HM_Density"] = EnzoInterface.HMDensity
FieldTypes["H2I_Density"] = EnzoInterface.H2IDensity
FieldTypes["H2II_Density"] = EnzoInterface.H2IIDensity
FieldTypes["DI_Density"] = EnzoInterface.DIDensity
FieldTypes["DII_Density"] = EnzoInterface.DIIDensity
FieldTypes["HDI_Density"] = EnzoInterface.HDIDensity
FieldTypes["Metal_Density"] = EnzoInterface.Metallicity
