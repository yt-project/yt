import yt

ts = yt.load("enzo_tiny_cosmology/DD000?/DD000?")
for ds in ts:
    p = yt.ProjectionPlot(ds, "z", ("gas", "density"))
    p.annotate_timestamp(corner="upper_left", redshift=True, draw_inset_box=True)
    p.annotate_scale(corner="upper_right")
    p.save()
