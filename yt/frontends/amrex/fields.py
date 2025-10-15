import re
from typing import TypeAlias

import numpy as np

from yt._typing import KnownFieldsT
from yt.fields.field_info_container import FieldInfoContainer
from yt.units import YTQuantity
from yt.utilities.physical_constants import amu_cgs, boltzmann_constant_cgs, c

rho_units = "code_mass / code_length**3"
mom_units = "code_mass / (code_time * code_length**2)"
eden_units = "code_mass / (code_time**2 * code_length)"  # erg / cm^3


def _thermal_energy_density(data):
    # What we've got here is UEINT:
    # u here is velocity
    # E is energy density from the file
    #   rho e = rho E - rho * u * u / 2
    ke = (
        0.5
        * (
            data["gas", "momentum_density_x"] ** 2
            + data["gas", "momentum_density_y"] ** 2
            + data["gas", "momentum_density_z"] ** 2
        )
        / data["gas", "density"]
    )
    return data["boxlib", "eden"] - ke


def _specific_thermal_energy(data):
    # This is little e, so we take thermal_energy_density and divide by density
    return data["gas", "thermal_energy_density"] / data["gas", "density"]


def _temperature(data):
    mu = data.ds.parameters["mu"]
    gamma = data.ds.parameters["gamma"]
    tr = data["gas", "thermal_energy_density"] / data["gas", "density"]
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
        def get_mass(data):
            species_mass = data.ds.index.parameters[ptype + "_mass"]
            return data[ptype, "particle_weight"] * YTQuantity(species_mass, "kg")

        self.add_field(
            (ptype, "particle_mass"),
            sampling_type="particle",
            function=get_mass,
            units="kg",
        )

        def get_charge(data):
            species_charge = data.ds.index.parameters[ptype + "_charge"]
            return data[ptype, "particle_weight"] * YTQuantity(species_charge, "C")

        self.add_field(
            (ptype, "particle_charge"),
            sampling_type="particle",
            function=get_charge,
            units="C",
        )

        def get_energy(data):
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

        def get_velocity_x(data):
            return (
                c**2
                * data[ptype, "particle_momentum_x"]
                / data[ptype, "particle_energy"]
            )

        def get_velocity_y(data):
            return (
                c**2
                * data[ptype, "particle_momentum_y"]
                / data[ptype, "particle_energy"]
            )

        def get_velocity_z(data):
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
            def velocity(data):
                return (
                    data[ptype, f"particle_momentum_{axis}"]
                    / data[ptype, "particle_mass"]
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
            def velocity(data):
                return data["boxlib", f"{axis}mom"] / data["boxlib", "density"]

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
            def momentum(data):
                return data["boxlib", f"{axis}vel"] * data["boxlib", "density"]

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
        ("MachNumber", ("", ["mach_number"], "Mach Number")),
        ("abar", ("", [], r"$\bar{A}$")),
        ("Ye", ("", [], r"$Y_e$")),
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
                sub = Substance(field)
                # Overwrite field to use nicer tex_label display_name
                self.add_output_field(
                    ("boxlib", field),
                    sampling_type="cell",
                    units="",
                    display_name=rf"X\left({sub.to_tex()}\right)",
                )
                self.alias(("gas", f"{sub}_fraction"), ("boxlib", field), units="")
                func = _create_density_func(("gas", f"{sub}_fraction"))
                self.add_field(
                    name=("gas", f"{sub}_density"),
                    sampling_type="cell",
                    function=func,
                    units=self.ds.unit_system["density"],
                    display_name=rf"\rho {sub.to_tex()}",
                )


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
        tfromp = False
        if "use_tfromp" in self.ds.parameters:
            # original MAESTRO (F90) code
            tfromp = self.ds.parameters["use_tfromp"]
        elif "maestro.use_tfromp" in self.ds.parameters:
            # new MAESTROeX (C++) code
            tfromp = self.ds.parameters["maestro.use_tfromp"]

        if tfromp:
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
                sub = Substance(field)
                # Overwrite field to use nicer tex_label display_name
                self.add_output_field(
                    ("boxlib", field),
                    sampling_type="cell",
                    units="",
                    display_name=rf"X\left({sub.to_tex()}\right)",
                )
                self.alias(("gas", f"{sub}_fraction"), ("boxlib", field), units="")
                func = _create_density_func(("gas", f"{sub}_fraction"))
                self.add_field(
                    name=("gas", f"{sub}_density"),
                    sampling_type="cell",
                    function=func,
                    units=unit_system["density"],
                    display_name=rf"\rho {sub.to_tex()}",
                )

            elif field.startswith("omegadot("):
                sub = Substance(field)
                display_name = rf"\dot{{\omega}}\left[{sub.to_tex()}\right]"
                # Overwrite field to use nicer tex_label'ed display_name
                self.add_output_field(
                    ("boxlib", field),
                    sampling_type="cell",
                    units=unit_system["frequency"],
                    display_name=display_name,
                )
                self.alias(
                    ("gas", f"{sub}_creation_rate"),
                    ("boxlib", field),
                    units=unit_system["frequency"],
                )


class QuokkaFieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        ("gasDensity", ("code_mass / code_length**3", ["density"], r"\rho")),
        (
            "gasEnergy",
            (
                "code_mass / (code_length * code_time**2)",
                ["total_energy_density"],
                r"\rho E",
            ),
        ),
        (
            "gasInternalEnergy",
            (
                "code_mass / (code_length * code_time**2)",
                ["internal_energy_density"],
                r"\rho e",
            ),
        ),
        (
            "x-GasMomentum",
            (
                "code_mass / (code_length**2 * code_time)",
                ["momentum_density_x"],
                r"\rho u",
            ),
        ),
        (
            "y-GasMomentum",
            (
                "code_mass / (code_length**2 * code_time)",
                ["momentum_density_y"],
                r"\rho v",
            ),
        ),
        (
            "z-GasMomentum",
            (
                "code_mass / (code_length**2 * code_time)",
                ["momentum_density_z"],
                r"\rho w",
            ),
        ),
        # Temperature field is not present in early Quokka datasets
        ("gasTemperature", ("K", ["temperature"], r"T")),
        # Scalar fields are not always present
        ("scalar_0", ("", ["scalar_0"], "Scalar 0")),
        ("scalar_1", ("", ["scalar_1"], "Scalar 1")),
        ("scalar_2", ("", ["scalar_2"], "Scalar 2")),
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
    )

    def setup_fluid_fields(self):
        # Define momentum-to-velocity conversion
        def _get_cell_velocity(axis):
            def velocity(field, data):
                # Divide momentum by density for cell-centered velocity
                return (
                    data["boxlib", f"{axis}-GasMomentum"] / data["boxlib", "gasDensity"]
                )

            return velocity

        # Add cell-centered velocity fields dynamically for each axis
        for axis in "xyz":
            field_name = f"velocity_{axis}"
            momentum_name = f"{axis}-GasMomentum"

            if ("boxlib", momentum_name) in self.field_list:
                self.add_field(
                    ("gas", field_name),
                    sampling_type="cell",
                    function=_get_cell_velocity(axis),
                    units=self.ds.unit_system["length"] / self.ds.unit_system["time"],
                    display_name=f"v_{axis} (cell-centered)",
                )

        # Add face-centered velocities dynamically for each axis
        for axis in "xyz":
            face_velocity_name = f"velocity_face_{axis}"

            if ("boxlib", f"{axis}-velocity") in self.field_list:
                self.add_field(
                    ("gas", face_velocity_name),
                    sampling_type="cell",
                    function=lambda field, data, axis=axis: data[
                        "boxlib", f"{axis}-velocity"
                    ]
                    * self.ds.unit_system["length"]
                    / self.ds.unit_system["time"],
                    units=self.ds.unit_system["length"] / self.ds.unit_system["time"],
                    display_name=f"v_{axis} (face-centered)",
                )

        # Call the Bfields and radiation fields setup
        self.setup_Bfields()  # magnetic fields are still placeholder here
        self.setup_radiation_fields()

    def setup_Bfields(self):
        """
        Dynamically add magnetic fields based on presence of Bfield fields in ds.parameters['fields']
        """
        # Check if any field name contains 'Bfield'
        if not any("Bfield" in field for field in self.ds.parameters.get("fields", [])):
            return

        for axis in "xyz":
            boxlib_bfield = f"{axis}-BField"

            if ("boxlib", boxlib_bfield) in self.field_list:
                self.add_field(
                    ("mag", f"{axis}-field"),
                    sampling_type="cell",
                    function=lambda field, data, axis=axis: data[
                        "boxlib", f"{axis}-BField"
                    ]
                    * self.ds.unit_system["magnetic_field_strength"],
                    units=self.ds.unit_system["magnetic_field_strength"],
                    display_name=f"B_{axis} (magnetic field)",
                )

    def setup_radiation_fields(self):
        # Dynamically add radiation fields
        num_groups = self.ds.parameters.get("radiation_field_groups", 0)
        for group in range(num_groups):
            # Add energy density
            energy_field = f"radEnergy-Group{group}"
            if ("boxlib", energy_field) in self.field_list:
                self.add_field(
                    ("rad", f"energy_density_{group}"),
                    sampling_type="cell",
                    function=lambda _, data, ef=energy_field: data["boxlib", ef]
                    * self.ds.unit_system["energy"]
                    / self.ds.unit_system["length"] ** 3,
                    units=self.ds.unit_system["energy"]
                    / self.ds.unit_system["length"] ** 3,
                    display_name=f"Radiation Energy Density Group {group}",
                )

            # Add flux density for each axis
            for axis in "xyz":
                flux_field = f"{axis}-RadFlux-Group{group}"
                if ("boxlib", flux_field) in self.field_list:
                    self.add_field(
                        ("rad", f"flux_density_{axis}_{group}"),
                        sampling_type="cell",
                        function=lambda field, data, flux_field=flux_field: data[
                            "boxlib", flux_field
                        ]
                        * self.ds.unit_system["energy"]
                        / self.ds.unit_system["length"] ** 2
                        / self.ds.unit_system["time"],
                        units=self.ds.unit_system["energy"]
                        / self.ds.unit_system["length"] ** 2
                        / self.ds.unit_system["time"],
                        display_name=f"Radiation Flux Density {axis.upper()} Group {group}",
                    )

        # Dynamically set up custom particle fields (convert `real_comp` to physics) for all particle types
        for ptype in self.ds.particle_types:
            self.setup_custom_particle_fields(ptype)

    def setup_custom_particle_fields(self, ptype):
        """
        Dynamically set up custom particle fields (real_comp) for the given particle type,
        interpreting dimensional expressions and applying appropriate unit conversions.

        Parameters:
        -----------
        ptype : str
            The particle type (e.g., 'Rad_particles', 'CIC_particles', etc.) for which fields will be dynamically added.
        """
        particle_info = self.ds.parameters.get("particle_info", {}).get(ptype, {})
        field_map = particle_info.get("fields", [])
        unit_map = particle_info.get("units", {})

        def interpret_units(dim_expr):
            """
            Interpret a dimensional expression like 'M^a L^b T^c Θ^d' and return the corresponding unit factor.

            Parameters:
            -----------
            dim_expr : str
                Dimensional expression (e.g., 'M^1 L^2 T^-3' for luminosity).

            Returns:
            --------
            conversion_factor : unyt_quantity
                The unit conversion factor constructed from the dimensional expression.
            """
            if not dim_expr or dim_expr == "dimensionless":
                return 1.0  # No conversion needed for dimensionless fields

            # Parse the dimensional expression
            dimensions = {"M": 0, "L": 0, "T": 0, "Θ": 0}  # Default to zero exponents
            for term in dim_expr.split():
                base, exponent = term.split("^")
                dimensions[base] = int(exponent)

            # Construct the conversion factor using the unit system. Use code_mass, code_length, code_time format for units
            conversion_factor = (
                f"code_mass**{dimensions['M']} * "
                f"code_length**{dimensions['L']} * "
                f"code_time**{dimensions['T']} * "
                f"code_temperature**{dimensions['Θ']}"
            )
            # Convert the string expression to actual unit quantity
            conversion_factor = self.ds.quan(1.0, conversion_factor)
            return conversion_factor

        for idx, field_name in enumerate(field_map):
            # Construct the `real_comp` field name
            real_comp_field = f"particle_real_comp{idx}"

            # Check if the `real_comp` field exists in the dataset
            if (ptype, real_comp_field) in self.field_list:
                # Retrieve the dimensional expression for the new field
                dim_expr = unit_map.get(field_name, "dimensionless")

                # Interpret the dimensional expression to get the conversion factor
                conversion_factor = interpret_units(dim_expr)

                # Add the field with the appropriate units
                self.add_field(
                    (ptype, field_name),
                    sampling_type="particle",
                    function=lambda field,
                    data,
                    real_comp_field=real_comp_field,
                    conv=conversion_factor: (data[ptype, real_comp_field] * conv),
                    units=conversion_factor.units
                    if hasattr(conversion_factor, "units")
                    else "",
                    display_name=field_name.replace("_", " ").capitalize(),
                )


substance_expr_re = re.compile(r"\(([a-zA-Z][a-zA-Z0-9]*)\)")
substance_elements_re = re.compile(r"(?P<element>[a-zA-Z]+)(?P<digits>\d*)")
SubstanceSpec: TypeAlias = list[tuple[str, int]]


class Substance:
    def __init__(self, data: str) -> None:
        if (m := substance_expr_re.search(data)) is None:
            raise ValueError(f"{data!r} doesn't match expected regular expression")
        sub_str = m.group()
        constituents = substance_elements_re.findall(sub_str)

        # 0 is used as a sentinel value to mark descriptive names
        default_value = 1 if len(constituents) > 1 else 0
        self._spec: SubstanceSpec = [
            (name, int(count or default_value)) for (name, count) in constituents
        ]

    def get_spec(self) -> SubstanceSpec:
        return self._spec.copy()

    def is_isotope(self) -> bool:
        return len(self._spec) == 1 and self._spec[0][1] > 0

    def is_molecule(self) -> bool:
        return len(self._spec) != 1

    def is_descriptive_name(self) -> bool:
        return len(self._spec) == 1 and self._spec[0][1] == 0

    def __str__(self) -> str:
        return "".join(
            f"{element}{count if count > 1 else ''}" for element, count in self._spec
        )

    def _to_tex_isotope(self) -> str:
        element, count = self._spec[0]
        return rf"^{{{count}}}{element}"

    def _to_tex_molecule(self) -> str:
        return "".join(
            rf"{element}_{{{count if count > 1 else ''}}}"
            for element, count in self._spec
        )

    def _to_tex_descriptive(self) -> str:
        return str(self)

    def to_tex(self) -> str:
        if self.is_isotope():
            return self._to_tex_isotope()
        elif self.is_molecule():
            return self._to_tex_molecule()
        elif self.is_descriptive_name():
            return self._to_tex_descriptive()
        else:
            # should only be reachable in case of a regular expression defect
            raise RuntimeError


def _create_density_func(field_name):
    def _func(data):
        return data[field_name] * data["gas", "density"]

    return _func
