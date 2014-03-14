.. function:: annotate_arrow(self, pos, code_size, plot_args=None):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.ArrowCallback`.)

   This adds an arrow pointing at *pos* with size
   *code_size* in code units.  *plot_args* is a dict fed to
   matplotlib with arrow properties.

.. python-script::

   from yt.mods import *
   pf = load("HiresIsolatedGalaxy/DD0044/DD0044")
   slc = SlicePlot(pf, 'z', 'density', width=(10,'kpc'), center='max')
   slc.annotate_arrow((0.53, 0.53, 0.53), 1/pf['kpc'])
   slc.save()

-------------

.. function:: annotate_clumps(self, clumps, plot_args=None):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.ClumpContourCallback`.)

   Take a list of *clumps* and plot them as a set of
   contours.

.. python-script::

   from yt.mods import *
   from yt.analysis_modules.level_sets.api import *

   pf = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   data_source = pf.disk([0.5, 0.5, 0.5], [0., 0., 1.],
                           8./pf.units['kpc'], 1./pf.units['kpc'])

   c_min = 10**na.floor(na.log10(data_source['density']).min()  )
   c_max = 10**na.floor(na.log10(data_source['density']).max()+1)

   function = 'self.data[\'Density\'].size > 20'
   master_clump = Clump(data_source, None, 'density', function=function)
   find_clumps(master_clump, c_min, c_max, 2.0)
   leaf_clumps = get_lowest_clumps(master_clump)

   prj = ProjectionPlot(pf, 2, 'density', center='c', width=(20,'kpc'))
   prj.annotate_clumps(leaf_clumps)
   prj.save('clumps')

-------------

.. function:: annotate_contour(self, field, ncont=5, factor=4, take_log=False, clim=None, plot_args=None):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.ContourCallback`.)

   Add contours in *field* to the plot.  *ncont* governs the
   number of contours generated, *factor* governs the number
   of points used in the interpolation, *take_log* governs
   how it is contoured and *clim* gives the (upper, lower)
   limits for contouring.

.. python-script::
   
   from yt.mods import *
   pf = load("Enzo_64/DD0043/data0043")
   s = SlicePlot(pf, "x", ["density"], center="max")
   s.annotate_contour("temperature")
   s.save()

-------------

.. function:: annotate_cquiver(self, field_x, field_y, factor):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.CuttingQuiverCallback`.)

   Get a quiver plot on top of a cutting plane, using
   *field_x* and *field_y*, skipping every *factor*
   datapoint in the discretization.

.. python-script::

   from yt.mods import *
   pf = load("Enzo_64/DD0043/data0043")
   s = OffAxisSlicePlot(pf, [1,1,0], ["density"], center="c")
   s.annotate_cquiver('CuttingPlaneVelocityX', 'CuttingPlaneVelocityY', 10)
   s.zoom(1.5)
   s.save()

-------------

.. function:: annotate_grids(self, alpha=1.0, min_pix=1, annotate=False, periodic=True):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.GridBoundaryCallback`.)

   Adds grid boundaries to a plot, optionally with
   *alpha*-blending. Cuttoff for display is at *min_pix*
   wide. *annotate* puts the grid id in the corner of the
   grid.  (Not so great in projections...)

.. python-script::

   from yt.mods import *
   pf = load("HiresIsolatedGalaxy/DD0044/DD0044")
   slc = SlicePlot(pf, 'z', 'density', width=(10,'kpc'), center='max')
   slc.annotate_grids()
   slc.save()

-------------

.. function:: annotate_hop_circles(self, hop_output, max_number=None, annotate=False, min_size=20, max_size=10000000, font_size=8, print_halo_size=False, print_halo_mass=False, width=None):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.HopCircleCallback`.)

   Accepts a :class:`yt.HopList` *hop_output* and plots up
   to *max_number* (None for unlimited) halos as circles.

.. python-script::

   from yt.mods import *
   pf = load("Enzo_64/DD0043/data0043")
   halos = HaloFinder(pf)
   p = ProjectionPlot(pf, "z", "density")
   p.annotate_hop_circles(halos)
   p.save()

-------------

.. function:: annotate_hop_particles(self, hop_output, max_number, p_size=1.0, min_size=20, alpha=0.2):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.HopParticleCallback`.)

   Adds particle positions for the members of each halo as
   identified by HOP. Along *axis* up to *max_number* groups
   in *hop_output* that are larger than *min_size* are
   plotted with *p_size* pixels per particle;  *alpha*
   determines the opacity of each particle.

.. python-script::

   from yt.mods import *
   pf = load("Enzo_64/DD0043/data0043")
   halos = HaloFinder(pf)
   p = ProjectionPlot(pf, "x", "density", center='m', width=(10, 'Mpc'))
   p.annotate_hop_particles(halos, max_number=100, p_size=5.0)
   p.save()

-------------

.. function:: annotate_image_line(self, p1, p2, data_coords=False, plot_args=None):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.ImageLineCallback`.)

   Plot from *p1* to *p2* (normalized image plane coordinates) with
   *plot_args* fed into the plot.

.. python-script::

   from yt.mods import *
   pf = load("HiresIsolatedGalaxy/DD0044/DD0044")
   p = ProjectionPlot(pf, 'z', 'density', center='m', width=(10, 'kpc'))
   p.annotate_image_line((0.3, 0.4), (0.8, 0.9), plot_args={'linewidth':5})
   p.save()

-------------

.. function:: annotate_line(self, x, y, plot_args=None):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.LinePlotCallback`.)

   Over plot *x* and *y* (in code units) with *plot_args* fed into the plot.

.. python-script::

   from yt.mods import *
   pf = load("HiresIsolatedGalaxy/DD0044/DD0044")
   p = ProjectionPlot(pf, 'z', 'density', center='m', width=(10, 'kpc'))
   p.annotate_line([-6, -4, -2, 0, 2, 4, 6], [3.6, 1.6, 0.4, 0, 0.4, 1.6, 3.6], plot_args={'linewidth':5})
   p.save()

-------------

.. function:: annotate_magnetic_field(self, factor=16, scale=None, scale_units=None, normalize=False):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.MagFieldCallback`.)

   Adds a 'quiver' plot of magnetic field to the plot,
   skipping all but every *factor* datapoint. *scale* is the
   data units per arrow length unit using *scale_units* (see
   matplotlib.axes.Axes.quiver for more info). if
   *normalize* is True, the magnetic fields will be scaled
   by their local (in-plane) length, allowing morphological
   features to be more clearly seen for fields with
   substantial variation in field strength.

.. python-script::

   from yt.mods import *
   pf = load("MHDSloshing/virgo_low_res.0054.vtk",
             parameters={"TimeUnits":3.1557e13, "LengthUnits":3.0856e24,
                         "DensityUnits":6.770424595218825e-27})
   p = ProjectionPlot(pf, 'z', 'density', center='c', width=(300, 'kpc'))
   p.annotate_magnetic_field()
   p.save()

-------------

.. function:: annotate_marker(self, pos, marker='x', plot_args=None):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.MarkerAnnotateCallback`.)

   Adds text *marker* at *pos* in code coordinates.
   *plot_args* is a dict that will be forwarded to the plot
   command.

.. python-script::

   from yt.mods import *
   pf = load("HiresIsolatedGalaxy/DD0044/DD0044")
   s = SlicePlot(pf, 'z', 'density', center='m', width=(10, 'kpc'))
   s.annotate_marker([0.53, 0.53, 0.53], plot_args={'s':10000})
   s.save()   

-------------

.. function:: annotate_particles(self, width, p_size=1.0, col='k', marker='o', stride=1.0, ptype=None, stars_only=False, dm_only=False, minimum_mass=None):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.ParticleCallback`.)

   Adds particle positions, based on a thick slab along
   *axis* with a *width* along the line of sight.  *p_size*
   controls the number of pixels per particle, and *col*
   governs the color.  *ptype* will restrict plotted
   particles to only those that are of a given type.
   *minimum_mass* will require that the particles be of a
   given mass, calculated via ParticleMassMsun, to be
   plotted.

.. python-script::

   from yt.mods import *
   pf = load("Enzo_64/DD0043/data0043")
   p = ProjectionPlot(pf, "x", "density", center='m', width=(10, 'Mpc'))
   p.annotate_particles(10/pf['Mpc'])
   p.save()

-------------

.. function:: annotate_point(self, pos, text, text_args=None):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.PointAnnotateCallback`.)

   This adds *text* at position *pos*, where *pos* is in
   code-space. *text_args* is a dict fed to the text
   placement code.

.. python-script::

   from yt.mods import *
   pf = load("HiresIsolatedGalaxy/DD0044/DD0044")
   p = ProjectionPlot(pf, 'z', 'density', center='m', width=(10, 'kpc'))
   p.annotate_point([0.53, 0.526, 0.53], "What's going on here?", text_args={'size':'xx-large', 'color':'w'})
   p.save()

-------------

.. function:: annotate_quiver(self, field_x, field_y, factor, scale=None, scale_units=None, normalize=False):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.QuiverCallback`.)

   Adds a 'quiver' plot to any plot, using the *field_x* and
   *field_y* from the associated data, skipping every
   *factor* datapoints *scale* is the data units per arrow
   length unit using *scale_units*  (see
   matplotlib.axes.Axes.quiver for more info)

.. python-script::

   from yt.mods import *
   pf = load("HiresIsolatedGalaxy/DD0044/DD0044")
   p = ProjectionPlot(pf, 'z', 'density', center=[0.53, 0.53, 0.53], 
                      weight_field='density', width=(20, 'kpc'))
   p.annotate_quiver('velocity_x', 'velocity_y', 16)
   p.save()

-------------

.. function:: annotate_sphere(self, center, radius, circle_args=None, text=None, text_args=None):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.SphereCallback`.)

   A sphere centered at *center* in code units with radius
   *radius* in code units will be created, with optional
   *circle_args*, *text*, and *text_args*.

.. python-script::

   from yt.mods import *
   pf = load("HiresIsolatedGalaxy/DD0044/DD0044")
   p = ProjectionPlot(pf, 'z', 'density', center=[0.53, 0.53, 0.53], width=(20, 'kpc'))
   p.annotate_sphere([0.53, 0.53, 0.53], 2/pf['kpc'], {'fill':True})
   p.save()

-------------

.. function:: annotate_streamlines(self, field_x, field_y, factor=6.0, nx=16, ny=16, xstart=(0, 1), ystart=(0, 1), nsample=256, start_at_xedge=False, start_at_yedge=False, plot_args=None):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.StreamlineCallback`.)

   Add streamlines to any plot, using the *field_x* and
   *field_y* from the associated data, using *nx* and *ny*
   starting points that are bounded by *xstart* and
   *ystart*.  To begin streamlines from the left edge of the
   plot, set *start_at_xedge* to True; for the bottom edge,
   use *start_at_yedge*.  A line with the qmean vector
   magnitude will cover 1.0/*factor* of the image.

.. python-script::

   from yt.mods import *
   pf = load("HiresIsolatedGalaxy/DD0044/DD0044")
   s = SlicePlot(pf, 'z', 'density', center=[0.53, 0.53, 0.53], width=(20, 'kpc'))
   s.annotate_streamlines('velocity_x', 'velocity_y')
   s.save()

-------------

.. function:: annotate_text(self, pos, text, data_coords=False, text_args=None):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.TextLabelCallback`.)

   Accepts a position in (0..1, 0..1) of the image, some
   text and optionally some text arguments. If data_coords
   is True, position will be in code units instead of image
   coordinates.

.. python-script::

   from yt.mods import *
   pf = load("HiresIsolatedGalaxy/DD0044/DD0044")
   s = SlicePlot(pf, 'z', 'density', center='m', width=(10, 'kpc'))
   s.annotate_text((0.53, 0.53), 'Sample text', text_args={'size':'xx-large', 'color':'w'})
   s.save()

-------------

.. function:: annotate_title(self, title='Plot'):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.TitleCallback`.)

   Accepts a *title* and adds it to the plot

.. python-script::

   from yt.mods import *
   pf = load("HiresIsolatedGalaxy/DD0044/DD0044")
   p = ProjectionPlot(pf, 'z', 'density', center=[0.53, 0.53, 0.53], width=(20, 'kpc'))
   p.annotate_title('Density plot')
   p.save()

-------------

.. function:: annotate_velocity(self, factor=16, scale=None, scale_units=None, normalize=False):

   (This is a proxy for :class:`~yt.visualization.plot_modifications.VelocityCallback`.)

   Adds a 'quiver' plot of velocity to the plot, skipping
   all but every *factor* datapoint. *scale* is the data
   units per arrow length unit using *scale_units* (see
   matplotlib.axes.Axes.quiver for more info). if
   *normalize* is True, the velocity fields will be scaled
   by their local (in-plane) length, allowing morphological
   features to be more clearly seen for fields with
   substantial variation in field strength (normalize is not
   implemented and thus ignored for Cutting Planes).

.. python-script::

   from yt.mods import *
   pf = load("HiresIsolatedGalaxy/DD0044/DD0044")
   p = SlicePlot(pf, 'z', 'density', center='m', width=(10, 'kpc'))
   p.annotate_velocity()
   p.save()
