"""
Tests for callbacks



"""
from __future__ import absolute_import

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
import os, tempfile, shutil
from yt.testing import \
    fake_amr_ds
import yt.units as u
from .test_plotwindow import assert_fname
from yt.visualization.api import \
    SlicePlot, ProjectionPlot, OffAxisSlicePlot
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
#  X grids
#    streamlines
#    units
#  X line
#    cquiver
#    clumps
#  X arrow
#  X marker
#  X sphere
#    hop_circles
#    hop_particles
#    coord_axes
#  X text
#    particles
#    title
#    flash_ray_data
#  X timestamp
#  X scale
#    material_boundary
#  X ray

@contextlib.contextmanager
def _cleanup_fname():
    tmpdir = tempfile.mkdtemp()
    yield tmpdir
    shutil.rmtree(tmpdir)

def test_timestamp_callback():
    with _cleanup_fname() as prefix:
        ax = 'z'
        vector = [1.0,1.0,1.0]
        ds = fake_amr_ds(fields = ("density",))
        p = ProjectionPlot(ds, ax, "density")
        p.annotate_timestamp()
        yield assert_fname, p.save(prefix)[0]
        p = SlicePlot(ds, ax, "density")
        p.annotate_timestamp()
        yield assert_fname, p.save(prefix)[0]
        p = OffAxisSlicePlot(ds, vector, "density")
        p.annotate_timestamp()
        yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", "density")
        p.annotate_timestamp(corner='lower_right', redshift=True, 
                             draw_inset_box=True)
        p.save(prefix)

def test_scale_callback():
    with _cleanup_fname() as prefix:
        ax = 'z'
        vector = [1.0,1.0,1.0]
        ds = fake_amr_ds(fields = ("density",))
        p = ProjectionPlot(ds, ax, "density")
        p.annotate_scale()
        yield assert_fname, p.save(prefix)[0]
        p = SlicePlot(ds, ax, "density")
        p.annotate_scale()
        yield assert_fname, p.save(prefix)[0]
        p = OffAxisSlicePlot(ds, vector, "density")
        p.annotate_scale()
        yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", "density")
        p.annotate_scale(corner='upper_right', coeff=10., unit='kpc')
        p.save(prefix)

def test_line_callback():
    with _cleanup_fname() as prefix:
        ax = 'z'
        vector = [1.0,1.0,1.0]
        ds = fake_amr_ds(fields = ("density",))
        p = ProjectionPlot(ds, ax, "density")
        p.annotate_line([0.1,0.1,0.1],[0.5,0.5,0.5])
        yield assert_fname, p.save(prefix)[0]
        p = SlicePlot(ds, ax, "density")
        p.annotate_line([0.1,0.1,0.1],[0.5,0.5,0.5])
        yield assert_fname, p.save(prefix)[0]
        p = OffAxisSlicePlot(ds, vector, "density")
        p.annotate_line([0.1,0.1,0.1],[0.5,0.5,0.5])
        yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", "density")
        p.annotate_line([0.1,0.1],[0.5,0.5], coord_system='axis', 
                        plot_args={'color':'red'})
        p.save(prefix)

def test_ray_callback():
    with _cleanup_fname() as prefix:
        ax = 'z'
        vector = [1.0,1.0,1.0]
        ds = fake_amr_ds(fields = ("density",))
        ray = ds.ray((0.1, 0.2, 0.3), (1.6, 1.8, 1.5))
        oray = ds.ortho_ray(0, (0.3, 0.4))
        p = ProjectionPlot(ds, ax, "density")
        p.annotate_ray(oray)
        p.annotate_ray(ray)
        yield assert_fname, p.save(prefix)[0]
        p = SlicePlot(ds, ax, "density")
        p.annotate_ray(oray)
        p.annotate_ray(ray)
        yield assert_fname, p.save(prefix)[0]
        p = OffAxisSlicePlot(ds, vector, "density")
        p.annotate_ray(oray)
        p.annotate_ray(ray)
        yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", "density")
        p.annotate_ray(oray)
        p.annotate_ray(ray, plot_args={'color':'red'})
        p.save(prefix)

def test_arrow_callback():
    with _cleanup_fname() as prefix:
        ax = 'z'
        vector = [1.0,1.0,1.0]
        ds = fake_amr_ds(fields = ("density",))
        p = ProjectionPlot(ds, ax, "density")
        p.annotate_arrow([0.5,0.5,0.5])
        yield assert_fname, p.save(prefix)[0]
        p = SlicePlot(ds, ax, "density")
        p.annotate_arrow([0.5,0.5,0.5])
        yield assert_fname, p.save(prefix)[0]
        p = OffAxisSlicePlot(ds, vector, "density")
        p.annotate_arrow([0.5,0.5,0.5])
        yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", "density")
        p.annotate_arrow([0.5,0.5], coord_system='axis', length=0.05)
        p.save(prefix)

def test_marker_callback():
    with _cleanup_fname() as prefix:
        ax = 'z'
        vector = [1.0,1.0,1.0]
        ds = fake_amr_ds(fields = ("density",))
        p = ProjectionPlot(ds, ax, "density")
        p.annotate_marker([0.5,0.5,0.5])
        yield assert_fname, p.save(prefix)[0]
        p = SlicePlot(ds, ax, "density")
        p.annotate_marker([0.5,0.5,0.5])
        yield assert_fname, p.save(prefix)[0]
        p = OffAxisSlicePlot(ds, vector, "density")
        p.annotate_marker([0.5,0.5,0.5])
        yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", "density")
        p.annotate_marker([0.5,0.5], coord_system='axis', marker='*')
        p.save(prefix)

def test_sphere_callback():
    with _cleanup_fname() as prefix:
        ax = 'z'
        vector = [1.0,1.0,1.0]
        ds = fake_amr_ds(fields = ("density",))
        p = ProjectionPlot(ds, ax, "density")
        p.annotate_sphere([0.5,0.5,0.5], 0.1)
        yield assert_fname, p.save(prefix)[0]
        p = SlicePlot(ds, ax, "density")
        p.annotate_sphere([0.5,0.5,0.5], 0.1)
        yield assert_fname, p.save(prefix)[0]
        p = OffAxisSlicePlot(ds, vector, "density")
        p.annotate_sphere([0.5,0.5,0.5], 0.1)
        yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", "density")
        p.annotate_sphere([0.5,0.5], 0.1, coord_system='axis', text='blah')
        p.save(prefix)

def test_text_callback():
    with _cleanup_fname() as prefix:
        ax = 'z'
        vector = [1.0,1.0,1.0]
        ds = fake_amr_ds(fields = ("density",))
        p = ProjectionPlot(ds, ax, "density")
        p.annotate_text([0.5,0.5,0.5], 'dinosaurs!')
        yield assert_fname, p.save(prefix)[0]
        p = SlicePlot(ds, ax, "density")
        p.annotate_text([0.5,0.5,0.5], 'dinosaurs!')
        yield assert_fname, p.save(prefix)[0]
        p = OffAxisSlicePlot(ds, vector, "density")
        p.annotate_text([0.5,0.5,0.5], 'dinosaurs!')
        yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", "density")
        p.annotate_text([0.5,0.5], 'dinosaurs!', coord_system='axis',
                        text_args={'color':'red'})
        p.save(prefix)

def test_velocity_callback():
    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields =
            ("density", "velocity_x", "velocity_y", "velocity_z"))
        for ax in 'xyz':
            p = ProjectionPlot(ds, ax, "density", weight_field="density")
            p.annotate_velocity()
            yield assert_fname, p.save(prefix)[0]
            p = SlicePlot(ds, ax, "density")
            p.annotate_velocity()
            yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", "density")
        p.annotate_velocity(factor=8, scale=0.5, scale_units="inches",
                            normalize = True)
        p.save(prefix)

def test_magnetic_callback():
    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields = ("density", "magnetic_field_x",
          "magnetic_field_y", "magnetic_field_z"))
        for ax in 'xyz':
            p = ProjectionPlot(ds, ax, "density", weight_field="density")
            p.annotate_magnetic_field()
            yield assert_fname, p.save(prefix)[0]
            p = SlicePlot(ds, ax, "density")
            p.annotate_magnetic_field()
            yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", "density")
        p.annotate_magnetic_field(factor=8, scale=0.5,
            scale_units="inches", normalize = True)
        p.save(prefix)

def test_quiver_callback():
    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields =
            ("density", "velocity_x", "velocity_y", "velocity_z"))
        for ax in 'xyz':
            p = ProjectionPlot(ds, ax, "density")
            p.annotate_quiver("velocity_x", "velocity_y")
            yield assert_fname, p.save(prefix)[0]
            p = ProjectionPlot(ds, ax, "density", weight_field="density")
            p.annotate_quiver("velocity_x", "velocity_y")
            yield assert_fname, p.save(prefix)[0]
            p = SlicePlot(ds, ax, "density")
            p.annotate_quiver("velocity_x", "velocity_y")
            yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", "density")
        p.annotate_quiver("velocity_x", "velocity_y", factor=8, scale=0.5,
            scale_units="inches", normalize = True,
            bv_x = 0.5 * u.cm / u.s,
            bv_y = 0.5 * u.cm / u.s)
        p.save(prefix)

def test_contour_callback():
    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields = ("density", "temperature"))
        for ax in 'xyz':
            p = ProjectionPlot(ds, ax, "density")
            p.annotate_contour("temperature")
            yield assert_fname, p.save(prefix)[0]
            p = ProjectionPlot(ds, ax, "density", weight_field="density")
            p.annotate_contour("temperature")
            yield assert_fname, p.save(prefix)[0]
            p = SlicePlot(ds, ax, "density")
            p.annotate_contour("temperature") # BREAKS WITH ndarray
            yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", "density")
        p.annotate_contour("temperature", ncont=10, factor=8,
            take_log=False, clim=(0.4, 0.6),
            plot_args={'lw':2.0}, label=True,
            text_args={'text-size':'x-large'})
        p.save(prefix)

        p = SlicePlot(ds, "x", "density")
        s2 = ds.slice(0, 0.2)
        p.annotate_contour("temperature", ncont=10, factor=8,
            take_log=False, clim=(0.4, 0.6),
            plot_args={'lw':2.0}, label=True,
            text_args={'text-size':'x-large'},
            data_source=s2)
        p.save(prefix)

def test_grids_callback():
    with _cleanup_fname() as prefix:
        ds = fake_amr_ds(fields = ("density",))
        for ax in 'xyz':
            p = ProjectionPlot(ds, ax, "density")
            p.annotate_grids()
            yield assert_fname, p.save(prefix)[0]
            p = ProjectionPlot(ds, ax, "density", weight_field="density")
            p.annotate_grids()
            yield assert_fname, p.save(prefix)[0]
            p = SlicePlot(ds, ax, "density")
            p.annotate_grids()
            yield assert_fname, p.save(prefix)[0]
        # Now we'll check a few additional minor things
        p = SlicePlot(ds, "x", "density")
        p.annotate_grids(alpha=0.7, min_pix=10, min_pix_ids=30,
            draw_ids=True, periodic=False, min_level=2,
            max_level=3, cmap="gist_stern")
        p.save(prefix)

