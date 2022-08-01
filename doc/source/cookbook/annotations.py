import yt

ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
p = yt.ProjectionPlot(ds, "z", ("gas", "density"))
p.annotate_sphere([0.54, 0.72], radius=(1, "Mpc"), coord_system="axis", text="Halo #7")
p.annotate_sphere(
    [0.65, 0.38, 0.3],
    radius=(1.5, "Mpc"),
    coord_system="data",
    circle_args={"color": "green", "linewidth": 4, "linestyle": "dashed"},
)
p.annotate_arrow([0.87, 0.59, 0.2], coord_system="data", color="red")
p.annotate_text([10, 20], "Some halos", coord_system="plot")
p.annotate_marker([0.45, 0.1, 0.4], coord_system="data", color="yellow", s=500)
p.annotate_line([0.2, 0.4], [0.3, 0.9], coord_system="axis")
p.annotate_timestamp(redshift=True)
p.annotate_scale()
p.save()
