import numpy as np

from yt.frontends.gizmo.fields import GizmoFieldInfo, metal_elements
from yt.testing import fake_particle_ds


class FakeData(dict):
    def __init__(self, ds, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ds = ds


def test_gizmo_nuclei_mass_density_fields():
    ds = fake_particle_ds()
    field_list = [
        ("PartType0", "Metallicity_00"),
        ("PartType0", "density"),
    ] + [("PartType0", f"{species}_metallicity") for species in metal_elements]

    fi = GizmoFieldInfo(ds, field_list)
    fi.setup_particle_fields("PartType0")

    for species in metal_elements:
        fkey = ("PartType0", f"{species}_nuclei_mass_density")
        assert fkey in fi

    data = FakeData(
        ds,
        {
            ("PartType0", "density"): np.array([1.0, 2.0]),
            ("PartType0", "He_metallicity"): np.array([0.25, 0.25]),
        },
    )

    he_field = fi["PartType0", "He_nuclei_mass_density"]
    res = he_field(data)
    np.testing.assert_allclose(res, [0.25, 0.5])
