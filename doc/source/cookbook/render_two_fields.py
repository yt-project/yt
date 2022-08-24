import yt
from yt.visualization.volume_rendering.api import Scene, create_volume_source

filePath = "Sedov_3d/sedov_hdf5_chk_0003"
ds = yt.load(filePath)
ds.force_periodicity()

sc = Scene()

# set up camera
cam = sc.add_camera(ds, lens_type="perspective")
cam.resolution = [400, 400]

cam.position = ds.arr([1, 1, 1], "cm")
cam.switch_orientation()

# add rendering of density field
dens = create_volume_source(ds, field=("flash", "dens"))
dens.use_ghost_zones = True
sc.add_source(dens)
sc.save("density.png", sigma_clip=6)

# add rendering of x-velocity field
vel = create_volume_source(ds, field=("flash", "velx"))
vel.use_ghost_zones = True
sc.add_source(vel)
sc.save("density_any_velocity.png", sigma_clip=6)
