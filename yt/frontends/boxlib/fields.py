import re
import string

import numpy as np

from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer
from yt.units import YTQuantity
from yt.utilities.physical_constants import amu_cgs, boltzmann_constant_cgs, c

rho_units = "code_mass / code_length**3"
mom_units = "code_mass / (code_time * code_length**2)"
eden_units = "code_mass / (code_time**2 * code_length)"  # erg / cm^3

spec_finder = re.compile(r".*\((\D*)(\d*)\).*")


def _thermal_energy_density(field, data):
    # What we've got here is UEINT:
    # u here is velocity
    # E is energy density from the file
    #   rho e = rho E - rho * u * u / 2
    ke = (
        0.5
        * (
            data[("gas", "momentum_density_x")] ** 2
            + data[("gas", "momentum_density_y")] ** 2
            + data[("gas", "momentum_density_z")] ** 2
        )
        / data[("gas", "density")]
    )
    return data[("boxlib", "eden")] - ke


def _specific_thermal_energy(field, data):
    # This is little e, so we take thermal_energy_density and divide by density
    return data[("gas", "thermal_energy_density")] / data[("gas", "density")]


def _temperature(field, data):
    mu = data.ds.parameters["mu"]
    gamma = data.ds.parameters["gamma"]
    tr = data[("gas", "thermal_energy_density")] / data[("gas", "density")]
    tr *= mu * amu_cgs / boltzmann_constant_cgs
    tr *= gamma - 1.0
    return tr


class WarpXFieldInfo(FieldInfoContainer):

    known_other_fields: KnownFieldsT = (
        ("Bx", ("T", ["magnetic_field_x", "B_x"], None)),
        ("By", ("T", ["magnetic_field_y", "B_y"], None)),
        ("Bz", ("T", ["magnetic_field_z", "B_z"], None)),
        ("Ex", ("V/m", ["electric_field_x", "E_x"], None)),
        ("Ey", ("V/m", ["electric_field_y", "E_y"], None)),
        ("Ez", ("V/m", ["electric_field_z", "E_z"], None)),
        ("jx", ("A", ["current_x", "Jx", "J_x"], None)),
        ("jy", ("A", ["current_y", "Jy", "J_y"], None)),
        ("jz", ("A", ["current_z", "Jz", "J_z"], None)),
    )

    known_particle_fields: KnownFieldsT = (
        ("particle_weight", ("", ["particle_weighting"], None)),
        ("particle_position_x", ("m", [], None)),
        ("particle_position_y", ("m", [], None)),
        ("particle_position_z", ("m", [], None)),
        ("particle_velocity_x", ("m/s", [], None)),
        ("particle_velocity_y", ("m/s", [], None)),
        ("particle_velocity_z", ("m/s", [], None)),
        ("particle_momentum_x", ("kg*m/s", [], None)),
        ("particle_momentum_y", ("kg*m/s", [], None)),
        ("particle_momentum_z", ("kg*m/s", [], None)),
    )

    extra_union_fields = (
        ("kg", "particle_mass"),
        ("C", "particle_charge"),
        ("", "particle_ones"),
    )

    def __init__(self, ds, field_list):
        super().__init__(ds, field_list)

        # setup nodal flag information
        for field in ds.index.raw_fields:
            finfo = self.__getitem__(("raw", field))
            finfo.nodal_flag = ds.nodal_flags[field]

    def setup_fluid_fields(self):
        for field in self.known_other_fields:
            fname = field[0]
            self.alias(("mesh", fname), ("boxlib", fname))

    def setup_fluid_aliases(self):
        super().setup_fluid_aliases("mesh")

    def setup_particle_fields(self, ptype):
        def get_mass(field, data):
            species_mass = data.ds.index.parameters[ptype + "_mass"]
            return data[(ptype, "particle_weight")] * YTQuantity(species_mass, "kg")

        self.add_field(
            (ptype, "particle_mass"),
            sampling_type="particle",
            function=get_mass,
            units="kg",
        )

        def get_charge(field, data):
            species_charge = data.ds.index.parameters[ptype + "_charge"]
            return data[(ptype, "particle_weight")] * YTQuantity(species_charge, "C")

        self.add_field(
            (ptype, "particle_charge"),
            sampling_type="particle",
            function=get_charge,
            units="C",
        )

        def get_energy(field, data):
            p2 = (
                data[ptype, "particle_momentum_x"] ** 2
                + data[ptype, "particle_momentum_y"] ** 2
                + data[ptype, "particle_momentum_z"] ** 2
            )
            return np.sqrt(p2 * c**2 + data[ptype, "particle_mass"] ** 2 * c**4)

        self.add_field(
            (ptype, "particle_energy"),
            sampling_type="particle",
            function=get_energy,
            units="J",
        )

        def get_velocity_x(field, data):
            return (
                c**2
                * data[ptype, "particle_momentum_x"]
                / data[ptype, "particle_energy"]
            )

        def get_velocity_y(field, data):
            return (
                c**2
                * data[ptype, "particle_momentum_y"]
                / data[ptype, "particle_energy"]
            )

        def get_velocity_z(field, data):
            return (
                c**2
                * data[ptype, "particle_momentum_z"]
                / data[ptype, "particle_energy"]
            )

        self.add_field(
            (ptype, "particle_velocity_x"),
            sampling_type="particle",
            function=get_velocity_x,
            units="m/s",
        )

        self.add_field(
            (ptype, "particle_velocity_y"),
            sampling_type="particle",
            function=get_velocity_y,
            units="m/s",
        )

        self.add_field(
            (ptype, "particle_velocity_z"),
            sampling_type="particle",
            function=get_velocity_z,
            units="m/s",
        )

        super().setup_particle_fields(ptype)


class NyxFieldInfo(FieldInfoContainer):
    known_particle_fields: KnownFieldsT = (
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
    )


class BoxlibFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        ("density", (rho_units, ["density"], None)),
        ("eden", (eden_units, ["total_energy_density"], None)),
        ("xmom", (mom_units, ["momentum_density_x"], None)),
        ("ymom", (mom_units, ["momentum_density_y"], None)),
        ("zmom", (mom_units, ["momentum_density_z"], None)),
        ("temperature", ("K", ["temperature"], None)),
        ("Temp", ("K", ["temperature"], None)),
        ("x_velocity", ("cm/s", ["velocity_x"], None)),
        ("y_velocity", ("cm/s", ["velocity_y"], None)),
        ("z_velocity", ("cm/s", ["velocity_z"], None)),
        ("xvel", ("cm/s", ["velocity_x"], None)),
        ("yvel", ("cm/s", ["velocity_y"], None)),
        ("zvel", ("cm/s", ["velocity_z"], None)),
    )

    known_particle_fields: KnownFieldsT = (
        ("particle_mass", ("code_mass", [], None)),
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_momentum_x", ("code_mass*code_length/code_time", [], None)),
        ("particle_momentum_y", ("code_mass*code_length/code_time", [], None)),
        ("particle_momentum_z", ("code_mass*code_length/code_time", [], None)),
        # Note that these are *internal* agmomen
        ("particle_angmomen_x", ("code_length**2/code_time", [], None)),
        ("particle_angmomen_y", ("code_length**2/code_time", [], None)),
        ("particle_angmomen_z", ("code_length**2/code_time", [], None)),
        ("particle_id", ("", ["particle_index"], None)),
        ("particle_mdot", ("code_mass/code_time", [], None)),
        # "mlast",
        # "r",
        # "mdeut",
        # "n",
        # "burnstate",
        # "luminosity",
    )

    def setup_particle_fields(self, ptype):
        def _get_vel(axis):
            def velocity(field, data):
                return (
                    data[(ptype, f"particle_momentum_{axis}")]
                    / data[(ptype, "particle_mass")]
                )

            return velocity

        for ax in "xyz":
            self.add_field(
                (ptype, f"particle_velocity_{ax}"),
                sampling_type="particle",
                function=_get_vel(ax),
                units="code_length/code_time",
            )

        super().setup_particle_fields(ptype)

    def setup_fluid_fields(self):
        unit_system = self.ds.unit_system
        # Now, let's figure out what fields are included.
        if any(f[1] == "xmom" for f in self.field_list):
            self.setup_momentum_to_velocity()
        elif any(f[1] == "xvel" for f in self.field_list):
            self.setup_velocity_to_momentum()
        self.add_field(
            ("gas", "specific_thermal_energy"),
            sampling_type="cell",
            function=_specific_thermal_energy,
            units=unit_system["specific_energy"],
        )
        self.add_field(
            ("gas", "thermal_energy_density"),
            sampling_type="cell",
            function=_thermal_energy_density,
            units=unit_system["pressure"],
        )
        if ("gas", "temperature") not in self.field_aliases:
            self.add_field(
                ("gas", "temperature"),
                sampling_type="cell",
                function=_temperature,
                units=unit_system["temperature"],
            )

    def setup_momentum_to_velocity(self):
        def _get_vel(axis):
            def velocity(field, data):
                return data[("boxlib", f"{axis}mom")] / data[("boxlib", "density")]

            return velocity

        for ax in "xyz":
            self.add_field(
                ("gas", f"velocity_{ax}"),
                sampling_type="cell",
                function=_get_vel(ax),
                units=self.ds.unit_system["velocity"],
            )

    def setup_velocity_to_momentum(self):
        def _get_mom(axis):
            def momentum(field, data):
                return data[("boxlib", f"{axis}vel")] * data[("boxlib", "density")]

            return momentum

        for ax in "xyz":
            self.add_field(
                ("gas", f"momentum_density_{ax}"),
                sampling_type="cell",
                function=_get_mom(ax),
                units=mom_units,
            )


class CastroFieldInfo(FieldInfoContainer):

    known_other_fields: KnownFieldsT = (
        ("density", ("g/cm**3", ["density"], r"\rho")),
        ("xmom", ("g/(cm**2 * s)", ["momentum_density_x"], r"\rho u")),
        ("ymom", ("g/(cm**2 * s)", ["momentum_density_y"], r"\rho v")),
        ("zmom", ("g/(cm**2 * s)", ["momentum_density_z"], r"\rho w")),
        # velocity components are not always present
        ("x_velocity", ("cm/s", ["velocity_x"], r"u")),
        ("y_velocity", ("cm/s", ["velocity_y"], r"v")),
        ("z_velocity", ("cm/s", ["velocity_z"], r"w")),
        ("rho_E", ("erg/cm**3", ["total_energy_density"], r"\rho E")),
        # internal energy density (not just thermal)
        ("rho_e", ("erg/cm**3", [], r"\rho e")),
        ("Temp", ("K", ["temperature"], r"T")),
        ("grav_x", ("cm/s**2", [], r"\mathbf{g} \cdot \mathbf{e}_x")),
        ("grav_y", ("cm/s**2", [], r"\mathbf{g} \cdot \mathbf{e}_y")),
        ("grav_z", ("cm/s**2", [], r"\mathbf{g} \cdot \mathbf{e}_z")),
        ("pressure", ("dyne/cm**2", [], r"p")),
        (
            "kineng",
            ("erg/cm**3", ["kinetic_energy_density"], r"\frac{1}{2}\rho|\mathbf{U}|^2"),
        ),
        ("soundspeed", ("cm/s", ["sound_speed"], "Sound Speed")),
        ("Machnumber", ("", ["mach_number"], "Mach Number")),
        ("entropy", ("erg/(g*K)", ["entropy"], r"s")),
        ("magvort", ("1/s", ["vorticity_magnitude"], r"|\nabla \times \mathbf{U}|")),
        ("divu", ("1/s", ["velocity_divergence"], r"\nabla \cdot \mathbf{U}")),
        ("eint_E", ("erg/g", [], r"e(E,U)")),
        ("eint_e", ("erg/g", [], r"e")),
        ("magvel", ("cm/s", ["velocity_magnitude"], r"|\mathbf{U}|")),
        ("radvel", ("cm/s", ["radial_velocity"], r"\mathbf{U} \cdot \mathbf{e}_r")),
        ("magmom", ("g*cm/s", ["momentum_magnitude"], r"\rho |\mathbf{U}|")),
        ("maggrav", ("cm/s**2", [], r"|\mathbf{g}|")),
        ("phiGrav", ("erg/g", [], r"\Phi")),
        ("enuc", ("erg/(g*s)", [], r"\dot{e}_{\rm{nuc}}")),
        ("rho_enuc", ("erg/(cm**3*s)", [], r"\rho \dot{e}_{\rm{nuc}}")),
        ("angular_momentum_x", ("g/(cm*s)", [], r"\ell_x")),
        ("angular_momentum_y", ("g/(cm*s)", [], r"\ell_y")),
        ("angular_momentum_z", ("g/(cm*s)", [], r"\ell_z")),
        ("phiRot", ("erg/g", [], r"\Phi_{\rm{rot}}")),
        ("rot_x", ("cm/s**2", [], r"\mathbf{f}_{\rm{rot}} \cdot \mathbf{e}_x")),
        ("rot_y", ("cm/s**2", [], r"\mathbf{f}_{\rm{rot}} \cdot \mathbf{e}_y")),
        ("rot_z", ("cm/s**2", [], r"\mathbf{f}_{\rm{rot}} \cdot \mathbf{e}_z")),
    )

    known_particle_fields: KnownFieldsT = (
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
    )

    def setup_fluid_fields(self):
        # add X's
        for _, field in self.ds.field_list:
            if field.startswith("X("):
                # We have a fraction
                nice_name, tex_label = _nice_species_name(field)
                self.alias(
                    ("gas", f"{nice_name}_fraction"), ("boxlib", field), units=""
                )
                func = _create_density_func(("gas", f"{nice_name}_fraction"))
                self.add_field(
                    name=("gas", f"{nice_name}_density"),
                    sampling_type="cell",
                    function=func,
                    units=self.ds.unit_system["density"],
                )
                # We know this will either have one letter, or two.
                if field[3] in string.ascii_letters:
                    element, weight = field[2:4], field[4:-1]
                else:
                    element, weight = field[2:3], field[3:-1]  # NOQA

                # Here we can, later, add number density
                # right now element and weight inferred above are unused


class MaestroFieldInfo(FieldInfoContainer):

    known_other_fields: KnownFieldsT = (
        ("density", ("g/cm**3", ["density"], None)),
        ("x_vel", ("cm/s", ["velocity_x"], r"\tilde{u}")),
        ("y_vel", ("cm/s", ["velocity_y"], r"\tilde{v}")),
        ("z_vel", ("cm/s", ["velocity_z"], r"\tilde{w}")),
        (
            "magvel",
            (
                "cm/s",
                ["velocity_magnitude"],
                r"|\tilde{\mathbf{U}} + w_0 \mathbf{e}_r|",
            ),
        ),
        (
            "radial_velocity",
            ("cm/s", ["radial_velocity"], r"\mathbf{U}\cdot \mathbf{e}_r"),
        ),
        ("circum_velocity", ("cm/s", ["tangential_velocity"], r"U - U\cdot e_r")),
        ("tfromp", ("K", [], "T(\\rho,p,X)")),
        ("tfromh", ("K", [], "T(\\rho,h,X)")),
        ("Machnumber", ("", ["mach_number"], "M")),
        ("S", ("1/s", [], None)),
        ("ad_excess", ("", [], r"\nabla - \nabla_\mathrm{ad}")),
        ("deltaT", ("", [], "[T(\\rho,h,X) - T(\\rho,p,X)]/T(\\rho,h,X)")),
        ("deltagamma", ("", [], r"\Gamma_1 - \overline{\Gamma_1}")),
        ("deltap", ("", [], "[p(\\rho,h,X) - p_0] / p_0")),
        ("divw0", ("1/s", [], r"\nabla \cdot \mathbf{w}_0")),
        # Specific entropy
        ("entropy", ("erg/(g*K)", ["entropy"], "s")),
        ("entropypert", ("", [], r"[s - \overline{s}] / \overline{s}")),
        ("enucdot", ("erg/(g*s)", [], r"\dot{\epsilon}_{nuc}")),
        ("Hext", ("erg/(g*s)", [], "H_{ext}")),
        # Perturbational pressure grad
        ("gpi_x", ("dyne/cm**3", [], r"\left(\nabla\pi\right)_x")),
        ("gpi_y", ("dyne/cm**3", [], r"\left(\nabla\pi\right)_y")),
        ("gpi_z", ("dyne/cm**3", [], r"\left(\nabla\pi\right)_z")),
        ("h", ("erg/g", [], "h")),
        ("h0", ("erg/g", [], "h_0")),
        # Momentum cannot be computed because we need to include base and
        # full state.
        ("momentum", ("g*cm/s", ["momentum_magnitude"], r"\rho |\mathbf{U}|")),
        ("p0", ("erg/cm**3", [], "p_0")),
        ("p0pluspi", ("erg/cm**3", [], r"p_0 + \pi")),
        ("pi", ("erg/cm**3", [], r"\pi")),
        ("pioverp0", ("", [], r"\pi/p_0")),
        # Base state density
        ("rho0", ("g/cm**3", [], "\\rho_0")),
        ("rhoh", ("erg/cm**3", ["enthalpy_density"], "(\\rho h)")),
        # Base state enthalpy density
        ("rhoh0", ("erg/cm**3", [], "(\\rho h)_0")),
        ("rhohpert", ("erg/cm**3", [], "(\\rho h)^\\prime")),
        ("rhopert", ("g/cm**3", [], "\\rho^\\prime")),
        ("soundspeed", ("cm/s", ["sound_speed"], None)),
        ("sponge", ("", [], None)),
        ("tpert", ("K", [], r"T - \overline{T}")),
        # Again, base state -- so we can't compute ourselves.
        ("vort", ("1/s", ["vorticity_magnitude"], r"|\nabla\times\tilde{U}|")),
        # Base state
        ("w0_x", ("cm/s", [], "(w_0)_x")),
        ("w0_y", ("cm/s", [], "(w_0)_y")),
        ("w0_z", ("cm/s", [], "(w_0)_z")),
    )

    def setup_fluid_fields(self):
        unit_system = self.ds.unit_system
        # pick the correct temperature field
        if self.ds.parameters["use_tfromp"]:
            self.alias(
                ("gas", "temperature"),
                ("boxlib", "tfromp"),
                units=unit_system["temperature"],
            )
        else:
            self.alias(
                ("gas", "temperature"),
                ("boxlib", "tfromh"),
                units=unit_system["temperature"],
            )

        # Add X's and omegadots, units of 1/s
        for _, field in self.ds.field_list:
            if field.startswith("X("):
                # We have a mass fraction
                nice_name, tex_label = _nice_species_name(field)
                # Overwrite field to use nicer tex_label display_name
                self.add_output_field(
                    ("boxlib", field),
                    sampling_type="cell",
                    units="",
                    display_name=tex_label,
                )
                self.alias(
                    ("gas", f"{nice_name}_fraction"), ("boxlib", field), units=""
                )
                func = _create_density_func(("gas", f"{nice_name}_fraction"))
                self.add_field(
                    name=("gas", f"{nice_name}_density"),
                    sampling_type="cell",
                    function=func,
                    units=unit_system["density"],
                    display_name=r"\rho %s" % tex_label,
                )

                # Most of the time our species will be of the form
                # element name + atomic weight (e.g. C12), but
                # sometimes we make up descriptive names (e.g. ash)
                if any(char.isdigit() for char in field):
                    # We know this will either have one letter, or two.
                    if field[3] in string.ascii_letters:
                        element, weight = field[2:4], field[4:-1]
                    else:
                        element, weight = field[2:3], field[3:-1]  # NOQA
                    weight = int(weight)

                # Here we can, later, add number density using 'element' and
                # 'weight' inferred above

            elif field.startswith("omegadot("):
                nice_name, tex_label = _nice_species_name(field)
                display_name = r"\dot{\omega}\left[%s\right]" % tex_label
                # Overwrite field to use nicer tex_label'ed display_name
                self.add_output_field(
                    ("boxlib", field),
                    sampling_type="cell",
                    units=unit_system["frequency"],
                    display_name=display_name,
                )
                self.alias(
                    ("gas", f"{nice_name}_creation_rate"),
                    ("boxlib", field),
                    units=unit_system["frequency"],
                )


def _nice_species_name(field):
    spec_match = spec_finder.search(field)
    nice_name = "".join(spec_match.groups())
    # if the species field is a descriptive name, then the match
    # on the integer will be blank
    # modify the tex string in this case to remove spurious tex spacing
    lab = r"X\left(^{%s}%s\right)"
    if spec_match.groups()[-1] == "":
        lab = r"X\left(%s%s\right)"
    tex_label = lab % spec_match.groups()[::-1]
    return nice_name, tex_label


def _create_density_func(field_name):
    def _func(field, data):
        return data[field_name] * data["gas", "density"]

    return _func
