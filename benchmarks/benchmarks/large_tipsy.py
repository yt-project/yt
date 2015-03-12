import numpy as np
import yt
from yt.utilities.answer_testing.framework import data_dir_load

class PKDGravTipsySuite:
    dsname = "halo1e11_run1.00400/halo1e11_run1.00400"
    def setup(self):
        cosmology_parameters = dict(current_redshift = 0.0,
                                    omega_lambda = 0.728,
                                    omega_matter = 0.272,
                                    hubble_constant = 0.702)
        kwargs = dict(field_dtypes = {"Coordinates": "d"},
                      cosmology_parameters = cosmology_parameters,
                      unit_base = {'length': (1.0/60.0, "Mpccm/h")},
                      n_ref = 64)
        self.ds = data_dir_load(self.dsname, yt.TipsyDataset, (), kwargs)

    def time_all_particles(self):
        dd = self.ds.all_data()
        dd["all", "particle_velocity_x"]
        dd["all", "particle_velocity_y"]
        dd["all", "particle_velocity_z"]

    def time_all_particles_derived(self):
        dd = self.ds.all_data()
        dd["all", "particle_velocity_magnitude"]

    def time_all_sph_kernel(self):
        dd = self.ds.all_data()
        dd["gas", "density"]

    def time_project_unweight(self):
        proj = self.ds.proj("density", 0)

    def time_project_weight(self):
        proj = self.ds.proj("density", 0, "density")

    def time_particle_quantities(self):
        dd = self.ds.all_data()
        dd.quantities.extrema("particle_mass")
        dd.quantities.extrema("particle_velocity_magnitude")
        dd.quantities.extrema(["particle_velocity_%s" % ax for ax in 'xyz'])

class GasolineTipsySuite(PKDGravTipsySuite):
    dsname = "agora_1e11.00400/agora_1e11.00400"
    def setup(self):
        cosmology_parameters = dict(current_redshift = 0.0,
                                    omega_lambda = 0.728,
                                    omega_matter = 0.272,
                                    hubble_constant = 0.702)
        kwargs = dict(cosmology_parameters = cosmology_parameters,
                      unit_base = {'length': (1.0/60.0, "Mpccm/h")},
                      n_ref = 64)
        self.ds = data_dir_load(self.dsname, yt.TipsyDataset, (), kwargs)
