from yt.fields.field_info_container import FieldInfoContainer

# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.


class CM1FieldInfo(FieldInfoContainer):
    known_other_fields = (
        # Each entry here is of the form
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
        ("uinterp", ("m/s", ["x velocity"], "uinterp")),
        ("vinterp", ("m/s", ["y velocity"], "vinterp")),
        ("winterp", ("m/s", ["z velocity"], "winterp")),
        ("u", ("m/s", ["x velocity"], "u")),
        ("v", ("m/s", ["y velocity"], "v")),
        ("w", ("m/s", ["z velocity"], "w")),
        ("hwin_sr", ("m/s", ["storm relative horizontal wind speed"], "hwind_sr")),
        ("windmag_sr", ("m/s", ["storm relative 3D wind speed"], "windmag_sr")),
        ("hwin_gr", ("m/s", ["ground relative horizontal wind speed"], "hwind_gr")),
        ("thpert", ("K", ["potential temperature perturbation"], r"$\theta'$")),
        ("th", ("K", ["potential temperature"], r"$\theta$")),
        ("thrhopert", ("K", ["density potential temperature perturbation"], r"$\theta_{\rho}'$")),
        ("prespert", ("hPa", ["presure perturbation"], "p'")),
        ("pipert", ("#", ["nondimensional pressure perturbation"], r"$\pi'$")),
        ("rhopert", ("kg/m**3", ["density perturbation"], r"\rho'")),
        ("dbz", ("dBZ", ["simulated reflectivity"], "dbZ")),
        ("qvpert", ("g/kg", ["water vapor mixing ratio perturbation"], "qv'")),
        ("qc", ("g/kg", ["cloud liquid water mixing ratio"], "qc")),
        ("qr", ("g/kg", ["rain mixing ratio"], "qr")),
        ("qi", ("g/kg", ["cloud ice mixing ratio"], "qi")),
        ("qs", ("g/kg", ["snow mixing ratio"], "qs")),
        ("qg", ("g/kg", ["graupel or hail mixing ratio"], "qg")),
        ("qcloud", ("g/kg", ["sum of cloud water and cloud ice mixing ratios"], "qc+qi")),
        ("qprecip", ("g/kg", ["sum of rain, graupel, and snow mixing ratios"], "qr+qs+qg")),
        ("nci", ("1/cm**3", ["cloud ice number concentration"], "nci")),
        ("ncr", ("1/cm**3", ["rain number concentration"], "ncr")),
        ("ncs", ("1/cm**3", ["snow number concentration"], "ncs")),
        ("ncg", ("1/cm**3", ["hail number concentration"], "ncg")),
        ("qhl", ("g/kg", ["hail mixing ratio (NSSL)"], "qhl")),
        ("ccn", ("1/kg", ["cloud condensation nuclei number concentration (NSSL)"], "ccn")),
        ("ccw", ("1/kg", ["cloud droplet number concentration (NSSL)"], "ccw")),
        ("crw", ("1/kg", ["rain number concentration (NSSL)"], "crw")),
        ("cci", ("1/kg", ["cloud ice number concentration (NSSL)"], "cci")),
        ("csw", ("1/kg", ["snow number concentration (NSSL)"], "csw")),
        ("chw", ("1/kg", ["graupel number concentration (NSSL)"], "chw")),
        ("chl", ("1/kg", ["hail number concentration (NSSL)"], "chl")),
        ("vhw", ("m**3/kg", ["graupel volume (NSSL)"], "vhw")),
        ("vhl", ("m**3/kg", ["hail volume (NSSL)"], "vhl")),
        ("xvort", ("1/s", ["x vorticity"], r"$\xi$")),
        ("xvort_tilt", ("1/s**2", ["x vorticity tilting"], r"$\xi$"+" tilting")),
        ("xvort_stretch", ("1/s**2", ["x vorticity stretching"], r"$\xi$"+" stretching")),
        ("xvort_baro", ("1/s**2", ["x vorticity baroclinic"], r"$\xi$"+" baroclinic")),
        ("xvort_solenoid", ("1/s**2", ["x vorticity solenoid"], r"$\xi$"+" solenoid")),
        ("yvort", ("1/s", ["y vorticity"], r"$\eta")),
        ("yvort_tilt", ("1/s**2", ["y vorticity tilting"], r"$\eta$"+" tilting")),
        ("yvort_stretch", ("1/s**2", ["y vorticity stretching"], r"$\eta$"+" stretching")),
        ("yvort_baro", ("1/s**2", ["y vorticity baroclinic"], r"$\eta$"+" baroclinic")),
        ("yvort_solenoid", ("1/s**2", ["y vorticity solenoid"], r"$\eta$"+" solenoid")),
        ("zvort", ("1/s", ["z vorticity"], r"$\zeta")),
        ("zvort_tilt", ("1/s**2", ["z vorticity tilting"], r"$\zeta$"+" tilting")),
        ("zvort_stretch", ("1/s**2", ["z vorticity stretching"], r"$\zeta$"+" stretching")),
        ("zvort_solenoid", ("1/s**2", ["z vorticity solenoid"], r"$\zeta$"+" solenoid")),
        ("hvort", ("1/s", ["horizontal vorticity magnitude"], r"$\omega_h$")),
        ("vortmag", ("1/s", ["vorticity magnitude"], r"$\omega$")),
        ("streamvort", ("1/s", ["streamwise vorticity"], "streamwise_vorticity")),
        ("khh", ("m**2/s", ["horizontal subgrid eddy diffusivity"], "khh")),
        ("khv", ("m**2/s", ["vertical subgrid eddy diffusivity"], "khv")),
        ("kmh", ("m**2/s", ["horizontal subgrid eddy viscosity"], "kmh")),
        ("kmv", ("m**2/s", ["vertical subgrid eddy viscosity"], "kmv")),
        ("pgrad_u", ("m/s**2", ["x pressure gradient acceleration"], "pgrad_u")),
        ("pgrad_v", ("m/s**2", ["y pressure gradient acceleration"], "pgrad_v")),
        ("pgrad_w", ("m/s**2", ["z pressure gradient acceleration"], "pgrad_w")),
        ("buoyancy", ("m/s**2", ["z buoyancy acceleration"], "buoyancy")),
    )

    known_particle_fields = (
        # Identical form to above
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
    )

    def setup_fluid_fields(self):
        # Here we do anything that might need info about the dataset.
        # You can use self.alias, self.add_output_field (for on-disk fields)
        # and self.add_field (for derived fields).
        pass

    def setup_particle_fields(self, ptype):
        super().setup_particle_fields(ptype)
        # This will get called for every particle type.
