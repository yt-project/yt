from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer

m_units = "Msun/h"
p_units = "kpccm/h"
r_units = "kpccm/h"
v_units = "km/s"


class AHFHalosFieldInfo(FieldInfoContainer):

    # See http://popia.ft.uam.es/AHF/files/AHF.pdf
    # and search for '*.AHF_halos'.
    known_particle_fields: KnownFieldsT = (
        ("ID", ("", ["particle_identifier"], None)),
        ("hostHalo", ("", [], None)),
        ("numSubStruct", ("", [], None)),
        ("Mvir", (m_units, ["particle_mass"], "Virial Mass")),
        ("npart", ("", [], None)),
        ("Xc", (p_units, ["particle_position_x"], None)),
        ("Yc", (p_units, ["particle_position_y"], None)),
        ("Zc", (p_units, ["particle_position_z"], None)),
        ("VXc", (v_units, ["particle_velocity_x"], None)),
        ("VYc", (v_units, ["particle_velocity_y"], None)),
        ("VZc", (v_units, ["particle_velocity_z"], None)),
        ("Rvir", (r_units, ["virial_radius"], "Virial Radius")),
        ("Rmax", (r_units, [], None)),
        ("r2", (r_units, [], None)),
        ("mbp_offset", (r_units, [], None)),
        ("com_offset", (r_units, [], None)),
        ("Vmax", (v_units, [], None)),
        ("v_sec", (v_units, [], None)),
        ("sigV", (v_units, [], None)),
        ("lambda", ("", [], None)),
        ("lambdaE", ("", [], None)),
        ("Lx", ("", [], None)),
        ("Ly", ("", [], None)),
        ("Lz", ("", [], None)),
        ("b", ("", [], None)),
        ("c", ("", [], None)),
        ("Eax", ("", [], None)),
        ("Eay", ("", [], None)),
        ("Eaz", ("", [], None)),
        ("Ebx", ("", [], None)),
        ("Eby", ("", [], None)),
        ("Ebz", ("", [], None)),
        ("Ecx", ("", [], None)),
        ("Ecy", ("", [], None)),
        ("Ecz", ("", [], None)),
        ("ovdens", ("", [], None)),
        ("nbins", ("", [], None)),
        ("fMhires", ("", [], None)),
        ("Ekin", ("Msun/h*(km/s)**2", [], None)),
        ("Epot", ("Msun/h*(km/s)**2", [], None)),
        ("SurfP", ("Msun/h*(km/s)**2", [], None)),
        ("Phi0", ("(km/s)**2", [], None)),
        ("cNFW", ("", [], None)),
    )
