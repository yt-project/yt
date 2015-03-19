import numpy as np
import yt

class SmallEnzoSuite:
    dsname = "IsolatedGalaxy/galaxy0030/galaxy0030"
    def setup(self):
        self.ds = yt.load(self.dsname)

    def time_all_particles(self):
        dd = self.ds.all_data()
        dd["all", "particle_velocity_x"]
        dd["all", "particle_velocity_y"]
        dd["all", "particle_velocity_z"]

    def time_all_particles_derived(self):
        dd = self.ds.all_data()
        dd["all", "particle_velocity_magnitude"]

    def time_gas_read(self):
        dd = self.ds.all_data()
        dd["gas", "density"]

    def time_gas_derived(self):
        dd = self.ds.all_data()
        dd["gas", "velocity_magnitude"]

    def time_project_unweight(self):
        proj = self.ds.proj("density", 0)

    def time_project_weight(self):
        proj = self.ds.proj("density", 0, "density")

    def time_ghostzones(self):
        dd = self.ds.all_data()
        dd["velocity_divergence"]

    def time_particle_quantities(self):
        dd = self.ds.all_data()
        dd.quantities.extrema("particle_mass")
        dd.quantities.extrema("particle_velocity_magnitude")
        dd.quantities.extrema(["particle_velocity_%s" % ax for ax in 'xyz'])

    def time_gas_quantites(self):
        dd = self.ds.all_data()
        dd.quantities.extrema("density")
        dd.quantities.extrema(["velocity_x", "velocity_y", "velocity_z"])
