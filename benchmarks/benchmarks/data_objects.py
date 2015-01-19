import numpy as np
import yt

if yt.__version__.startswith('2'):
    from yt.mods import load, ColorTransferFunction
# else:
#     from yt.visualization.volume_rendering.old_camera import Camera

DSNAME = "HiresIsolatedGalaxy/DD0044/DD0044"

class Suite:
    def setup(self):
        if yt.__version__.startswith('3'):
            self.ds = yt.load(DSNAME)
            self.ad = self.ds.all_data()
            self.field_name = "density"
        else:
            self.ds = load(DSNAME)
            self.ad = self.ds.h.all_data()
            self.field_name = "Density"
        # Warmup hdd
        self.ad[self.field_name]
        if yt.__version__.startswith('3'):
            mi, ma = self.ad.quantities['Extrema'](self.field_name)
            self.tf = yt.ColorTransferFunction((np.log10(mi)+1, np.log10(ma)))
        else:
            mi, ma = self.ad.quantities['Extrema'](self.field_name)[0]
            self.tf = ColorTransferFunction((np.log10(mi)+1, np.log10(ma)))
        self.tf.add_layers(5, w=0.02, colormap="spectral")
        self.c = [0.5, 0.5, 0.5]
        self.L = [0.5, 0.2, 0.7]
        self.W = 1.0
        self.Npixels = 512

    if yt.__version__.startswith('3'):
        def time_load_all_data(self):
            self.ds.all_data()
    else:
        def time_load_all_data(self):
            self.ds.h.all_data()

    def time_extrema_quantities(self):
        self.ad.quantities['Extrema'](self.field_name)

    if yt.__version__.startswith('3'):
        def time_alldata_projection(self):
            self.ds.proj(self.field_name, 0)
    else:
        def time_alldata_projection(self):
           self.ds.h.proj(0, self.field_name)

    if yt.__version__.startswith('3'):
        def time_slice(self):
            slc = self.ds.slice(0, 0.5)
            slc[self.field_name]
    else:
        def time_slice(self):
            slc = self.ds.h.slice(0, 0.5, self.field_name)
            slc[self.field_name]

#    if yt.__version__.startswith('3'):
#        def command(self):
#            cam = Camera(self.c, self.L, self.W, self.Npixels, self.tf, ds=self.ds)
#            cam.snapshot("%s_volume_rendered.png" % self.ds, clip_ratio=8.0)
#    else:
#        def command(self):
#            cam = self.ds.h.camera(self.c, self.L, self.W, self.Npixels, self.tf)
#            cam.snapshot("%s_volume_rendered.png" % self.ds, clip_ratio=8.0)
