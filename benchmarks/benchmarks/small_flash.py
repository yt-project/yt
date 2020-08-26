import numpy as np
import yt


class SmallFlashSuite:
    dsname = "GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_0690"

    def setup(self):
        self.ds = yt.load(self.dsname)

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
        dd = self.ds.sphere(self.ds.domain_center, self.ds.domain_width[0] * 0.25)
        dd["velocity_divergence"]

    def time_gas_quantities(self):
        dd = self.ds.all_data()
        dd.quantities.extrema("density")
        dd.quantities.extrema(["velocity_x", "velocity_y", "velocity_z"])
