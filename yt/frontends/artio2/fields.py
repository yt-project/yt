from yt.data_objects.field_info_container import \
    FieldInfoContainer, \
    FieldInfo, \
    NullFunc, \
    ValidateParameter, \
    ValidateDataField, \
    ValidateProperty, \
    ValidateSpatial, \
    ValidateGridType
import yt.data_objects.universal_fields
import numpy as np

KnownArtio2Fields = FieldInfoContainer()
add_artio2_field = KnownArtio2Fields.add_field

Artio2FieldInfo = FieldInfoContainer.create_with_fallback(FieldInfo)
add_field = Artio2FieldInfo.add_field

# maybe?
#translation_dict = {
#    "Density":("HVAR_GAS_DENSITY","Density"
#    "GasEnergy":"HVAR_GAS_ENERGY",
#    "Pressure":"HVAR_PRESSURE",
#    "Gamma":"HVAR_GAMMA",
#    "InternalEnergy":"HVAR_INTERNAL_ENERGY",
#    "XMomentumDensity":"HVAR_MOMENTUM_X",
#    "YMomentumDensity":"HVAR_MOMENTUM_Y",
#    "ZMomentumDensity":"HVAR_MOMENTUM_Z",
#    "ElectronInternalEnergy":"HVAR_ELECTRON_INTERNAL_ENERGY",
#    "HI_Density":"RT_HVAR_HI",
#    "HII_Density":"RT_HVAR_HII",
#    "HeI Density":"RT_HVAR_HeI",
#    "HeII_Density":"RT_HVAR_HeII",
#    "HeIII Density":"RT_HVAR_HeIII",
#    "H2_Density":"RT_HVAR_H2",
#    "MetalDensitySNII":"HVAR_METAL_DENSITY_II",
#    "MetalDensitySNIa":"HVAR_METAL_DENSITY_Ia",
#    "BlastwaveTime":"HVAR_BLASTWAVE_TIME",
#    "GravPotential":"VAR_POTENTIAL",
#    "GravPotentialOld":"VAR_POTENTIAL_HYDRO" }
#
#for yt_name,artio_name in translation_dict.items() :
#    add_field(

add_artio2_field( "Density", function=NullFunc, take_log=True,
            validators = [ValidateDataField("Density")])
