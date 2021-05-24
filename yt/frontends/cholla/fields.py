from yt.fields.field_info_container import FieldInfoContainer

# Copied from Athena frontend
pres_units = "code_pressure"
erg_units = "code_mass * (code_length/code_time)**2"
rho_units = "code_mass / code_length**3"
mom_units = "code_mass / code_length**2 / code_time"

def velocity_field(comp):
    def _velocity(field, data):
        return data["cholla", f"momentum_{comp}"] / data["cholla", "density"]

    return _velocity


# We need to specify which fields we might have in our dataset.  The field info
# container subclass here will define which fields it knows about.  There are
# optionally methods on it that get called which can be subclassed.


class ChollaFieldInfo(FieldInfoContainer):
    known_other_fields = (
        # Each entry here is of the form
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
        ("density", (rho_units, ["density"], None)),
        ("momentum_x", (mom_units, ["momentum_x"], None)),
        ("momentum_y", (mom_units, ["momentum_y"], None)),
        ("momentum_z", (mom_units, ["momentum_z"], None)),
        ("Energy", ("code_pressure", ["total_energy_density"], None)),
    )

    known_particle_fields = (
        # Identical form to above
        # ( "name", ("units", ["fields", "to", "alias"], # "display_name")),
    )

    def __init__(self, ds, field_list):
        super().__init__(ds, field_list)
        # If you want, you can check self.field_list

    # In Cholla, conservative variables are written out.
    # By default, yt concerns itself with primitive variables. The following
    # field definitions allow for conversions to primitive variables.


    def setup_fluid_fields(self):
        # Here we do anything that might need info about the dataset.
        # You can use self.alias, self.add_output_field (for on-disk fields)
        # and self.add_field (for derived fields).

        unit_system = self.ds.unit_system

        # Add velocity fields
        for comp in "xyz":
              self.add_field(
                  ("gas", f"velocity_{comp}"),
                  sampling_type="cell",
                  function=velocity_field(comp),
                  units=unit_system["velocity"],
              )

        def kinetic_energy_density(data):
            return (
                     0.5*data["cholla", "density"]
                     * (data["gas", "velocity_x"] * data["gas", "velocity_x"]
                      + data["gas", "velocity_y"] * data["gas", "velocity_y"]
                      + data["gas", "velocity_z"] * data["gas", "velocity_z"]
                     )
            )

        # Add pressure field
        if ("cholla", "GasEnergy") in self.field_list:
           self.add_output_field(
               ("cholla", "GasEnergy"), sampling_type="cell", units=pres_units
           )
           self.alias(
               ("gas", "thermal_energy"),
               ("cholla", "GasEnergy"),
               units=unit_system["pressure"],
           )

           def _pressure(field, data):
               return (
                  (data.ds.gamma - 1.0)
                  * data["cholla", "GasEnergy"]
               )
        else:
          def _pressure(field, data):
              return (
                  (data.ds.gamma - 1.0)
                  * (data["cholla", "Energy"] - 
                     kinetic_energy_density(data)
                    )
                  )

        self.add_field(
            ("gas", "pressure"),
            sampling_type="cell",
            function=_pressure,
            units=unit_system["pressure"],
        )
        
        def _specific_total_energy(field, data):
            return data["cholla", "Energy"] / data["cholla", "density"]

        self.add_field(
            ("gas", "specific_total_energy"),
            sampling_type="cell",
            function=_specific_total_energy,
            units=unit_system["specific_energy"],
        )


        # Add temperature field
        def _temperature(field, data):
            return (
                data.ds.mu
                * data["gas", "pressure"]
                / data["gas", "density"]
                * mh
                / kboltz
            )

        self.add_field(
            ("gas", "temperature"),
            sampling_type="cell",
            function=_temperature,
            units=unit_system["temperature"],
        )
        

    def setup_particle_fields(self, ptype):
        super().setup_particle_fields(ptype)
        # This will get called for every particle type.
