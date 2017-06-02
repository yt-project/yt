import yt
from matplotlib.animation import FuncAnimation
from matplotlib import rc_context

ts = yt.load('GasSloshingLowRes/sloshing_low_res_hdf5_plt_cnt_*')

plot = yt.SlicePlot(ts[0], 'z', 'density')
plot.set_zlim('density', 8e-29, 3e-26)

fig = plot.plots['density'].figure

# animate must accept an integer frame number. We use the frame number
# to identify which dataset in the time series we want to load
def animate(i):
    ds = ts[i]
    plot._switch_ds(ds)

animation = FuncAnimation(fig, animate, frames=len(ts))

# Override matplotlib's defaults to get a nicer looking font
with rc_context({'mathtext.fontset': 'stix'}):
    animation.save('animation.mp4')
