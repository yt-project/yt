from functools import partial

ax = "txyz"


def metric(mu: int, nu: int):
    if (mu, nu) == (0, 0):
        return -1
    elif mu == nu:
        return 1
    else:
        return 0


def setup_stress_energy_ideal(registry, ftype="gas"):

    pc = registry.ds.units.physical_constants
    inv_c2 = 1.0 / (pc.clight * pc.clight)

    def _T(field, data, mu: int, nu: int):
        Umu = data[ftype, f"four_velocity_{ax[mu]}"]
        Unu = data[ftype, f"four_velocity_{ax[nu]}"]
        p = data[ftype, "pressure"]
        e = data[ftype, "thermal_energy_density"]
        rho = data[ftype, "density"]
        return (rho + (e + p) * inv_c2) * Umu * Unu + metric(mu, nu) * p

    for mu in range(5):
        for nu in range(5):
            registry.add_field(
                (ftype, f"T^{mu}{nu}"),
                sampling_type="local",
                function=partial(_T, mu=mu, nu=nu),
                units=registry.ds.unit_system["pressure"],
            )
