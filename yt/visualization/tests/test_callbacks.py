"""
Tests for callbacks



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import os, tempfile, shutil
from yt.testing import \
    fake_amr_pf
from test_plotwindow import assert_fname
from yt.visualization.api import \
    SlicePlot, ProjectionPlot, OffAxisSlicePlot, OffAxisProjectionPlot
import contextlib

# These are a very simple set of tests that verify that each callback is or is
# not working.  They only check that it functions without an error; they do not
# check that it is providing correct results.

# These are the callbacks still to test:
#
#  X velocity
#  X magnetic_field
#  X quiver
#  X contour
#    grids
#    streamlines
#    axis_label
#    units
#    line
#    image_line
#    cquiver
#    clumps
#    arrow
#    point
#    marker
#    sphere
#    hop_circles
#    hop_particles
#    coord_axes
#    text
#    particles
#    title
#    flash_ray_data
#    timestamp
#    material_boundary

@contextlib.contextmanager
def _cleanup_fname():
    tmpdir = tempfile.mkdtemp()
    yield tmpdir
    shutil.rmtree(tmpdir)

def test_velocity_callback():
    with _cleanup_fname() as prefix:
        pf = fake_amr_pf(fields =
            ("Density", "x-velocity", "y-velocity", "z-velocity"))
        for ax in 'xyz':
            p = ProjectionPlot(pf, ax, "Density")
            p.annotate_velocity()
            yield assert_fname, p.save(prefix)[0]
            p = ProjectionPlot(pf, ax, "Density", weight_field="Density")
            p.annotate_velocity()
            yield assert_fname, p.save(prefix)[0]
            p = SlicePlot(pf, ax, "Density")
            p.annotate_velocity()
            yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(pf, "x", "Density")
        p.annotate_velocity(factor=8, scale=0.5, scale_units="inches",
                            normalize = True)
        p.save()

def test_magnetic_callback():
    with _cleanup_fname() as prefix:
        pf = fake_amr_pf(fields = ("Density", "Bx", "By", "Bz"))
        for ax in 'xyz':
            p = ProjectionPlot(pf, ax, "Density")
            p.annotate_magnetic_field()
            yield assert_fname, p.save(prefix)[0]
            p = ProjectionPlot(pf, ax, "Density", weight_field="Density")
            p.annotate_magnetic_field()
            yield assert_fname, p.save(prefix)[0]
            p = SlicePlot(pf, ax, "Density")
            p.annotate_magnetic_field()
            yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(pf, "x", "Density")
        p.annotate_magnetic_field(factor=8, scale=0.5,
            scale_units="inches", normalize = True)
        p.save()

def test_quiver_callback():
    with _cleanup_fname() as prefix:
        pf = fake_amr_pf(fields = ("Density", "Vx", "Vy"))
        for ax in 'xyz':
            p = ProjectionPlot(pf, ax, "Density")
            p.annotate_quiver("Vx", "Vy")
            yield assert_fname, p.save(prefix)[0]
            p = ProjectionPlot(pf, ax, "Density", weight_field="Density")
            p.annotate_quiver("Vx", "Vy")
            yield assert_fname, p.save(prefix)[0]
            p = SlicePlot(pf, ax, "Density")
            p.annotate_quiver("Vx", "Vy")
            yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(pf, "x", "Density")
        p.annotate_quiver("Vx", "Vy", factor=8, scale=0.5,
            scale_units="inches", normalize = True,
            bv_x = 0.5, bv_y = 0.5)
        p.save()

def test_contour_callback():
    with _cleanup_fname() as prefix:
        pf = fake_amr_pf(fields = ("Density", "Temperature"))
        for ax in 'xyz':
            p = ProjectionPlot(pf, ax, "Density")
            p.annotate_contour("Temperature")
            yield assert_fname, p.save(prefix)[0]
            p = ProjectionPlot(pf, ax, "Density", weight_field="Density")
            p.annotate_contour("Temperature")
            yield assert_fname, p.save(prefix)[0]
            p = SlicePlot(pf, ax, "Density")
            p.annotate_contour("Temperature")
            yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(pf, "x", "Density")
        p.annotate_contour("Temperature", ncont=10, factor=8,
            take_log=False, clim=(0.4, 0.6),
            plot_args={'lw':2.0}, label=True,
            label_args={'text-size':'x-large'})
        p.save()

        p = SlicePlot(pf, "x", "Density")
        s2 = pf.h.slice(0.2, 0)
        p.annotate_contour("Temperature", ncont=10, factor=8,
            take_log=False, clim=(0.4, 0.6),
            plot_args={'lw':2.0}, label=True,
            label_args={'text-size':'x-large'},
            data_source=s2)
        p.save()
