yt Concepts and History
=======================

.. sectionauthor:: J. S. Oishi <jsoishi@gmail.com>

History
-------

yt was originally written by Matthew Turk in the course of his graduate
studies.  However, it is now a community-developed project with contributions
from many people, the hospitality of several institutions, and benefiting from
many different grants.  With this community-driven approach and contributions
from many external developers, it has evolved from a simple data-reader and
exporter into what a fully-featured toolkit for analysis and visualization of
adaptive mesh refinement data.

yt was designed to be a completely Free (as in beer *and* as in freedom --
"free and libre" as the saying goes) user-extensible framework for analyzing
and visualizing astrophysical data, currently working with several different
codes, including the "flagship" codes Enzo, Orion, Nyx and FLASH.  It relies on
no proprietary software -- although it can be and has been extended to
interface with proprietary software and libraries -- and has been designed from
the ground up to enable users to be as immersed in the data as they desire.

yt is currently being developed by a team of developers.  We used to include
the entire list here, but it got too long, so now the best place to find it is
to take a look at the current `CREDITS
<http://hg.yt-project.org/yt/src/yt/CREDITS>`_ file.  All development is
conducted in the open, accessible at http://yt-project.org/ .


What functionality does yt offer?
---------------------------------

yt has evolved substantially over the time of its development.  Here is a
**non-comprehensive** list of features:

* Data Objects

  * Arbitrary data objects (Spheres, cylinders, rectangular prisms, arbitrary index selection)
  * Covering grids (smoothed and raw) for automatic ghost-zone generation
  * Identification of topologically-connected sets in arbitrary fields
  * Projections, orthogonal slices, oblique slices
  * Axially-aligned rays
  * Memory-conserving 1-, 2- and 3-D distribution functions of *arbitrary* fields and objects
  * Halo-finding (HOP) algorithm with full particle information and sphere access
  * Nearly **all** operations can be conducted in parallel

* Data access

  * Arbitrary field definition
  * Derived quantities (average values, spin parameter, bulk velocity, etc)
  * Custom C-written HDF5 backend for packed and unpacked AMR, NumPy-based HDF4 backend
  * Flexible units system with CGS by default
  * Per-user field and quantity plugins

* Visualization

  * Mathtext TeX-like text formatting
  * Slices, projections, oblique slices
  * Profiles and phase diagrams
  * Linked zooms, colormaps, and saving across multiple plots
  * Contours, vector plots, annotated boxes, grid boundary plot overlays.
  * Parallel volume rendering with arbitrary emission and absorption
    coefficients
  * Off-axis projections and line integrals

* Analysis Modules

  * `Parallel Halo Finding <http://adsabs.harvard.edu/abs/2010ApJS..191...43S>`_
  * Light cone generation
  * Halo merger tree creation
  * Level set extraction and clump identification
  * Halo Mass Function
  * SED generation
  * Arbitrary global or sub-region two-point functions

* Command-line tools

  * Time-series movies
  * HOP Halo Finding
  * Quick slice and projection images
  * Mapserver
  * Bug reports
  * IPython frontend

* Access to components

  * Monetary cost: **FREE**.
  * Source code availability: **FULL**.
  * Portability: **YES**.


The Goals and Design of yt
--------------------------

There are many tools available for analysis and visualization of AMR
data; there are many just for ``enzo``. So why yt? Along the road
to answering that question, we shall take a somewhat philosophical
scenic route. For the more pragmatically minded, the answer is simple:
what yt does not yet do, you can make it do so. This is not as
glib as it may seem: it is in fact the main philosophical tennant that
underlies yt. In this section, it is not our goal to show you just
how much yt already does. Instead, we will discuss how it is that
yt does anything at all. In doing so, we hope to give you a sense
of whether or not yt will align with your science goals.

At its core, yt is not a set of scripts to visualize AMR data, nor
is it a set of low-level routines that return a homo- or even
heterogeneous set of gridded data to your favorite scientific
programming language--though yt incorporates both of these things,
if your favorite scientific language is python. Instead, yt
provides a series of objects, some common AMR code structures (such as
hierarchies and levels in a nested mesh) and some physical (a
cylinder, cube, or sphere somewhere in the problem domain), that allow
you to process AMR data in order to get at the fundamental underlying
physics. 


yt evolved naturally out of three design goals, though when Matt
was busy writing it, he never really thought about them.  Over
time, it became clear that they are real and furthermore that they
are important to understanding how to use yt.  These three goals
are directed analysis, repeatability, and data exploration. 

Directed Analysis: Answering a Question
+++++++++++++++++++++++++++++++++++++++

One of the main motivators for yt is to make it possible to sit
down with a definite question about an AMR dataset and code up a
script that will provide an answer to that question. Indeed much of its
object-oriented nature can be viewed as a way perform operations on a
data object. Given that AMR simulations are usually used to track some
kind of structure formation, be it shocks, stars, or galaxies, the
data object may not be the entire domain, but some region within it
that is interesting. This data object in hand, yt makes it easy
(if possible: some tasks yt can merely make *possible*) to
manipulate that data in such a way to answer a question related to
your research.

Repeatability
+++++++++++++

In any scientific analysis, being able to repeat the set of steps that
prepared an answer or physical quantity is essential.  To that end,
much of the usage of yt is focused around running scripts,
describing actions and plots programmatically.  Being able to write a
script or conducting a set of commands that will reproduce identical
results is fundamentally important, and yt will attempt to make
that easy.  It's for this reason that the interactive features of
yt are not always as advanced as they might otherwise be. We are
actively working on integrating the SAGE notebook system into yt,
which our preliminary tests suggest is a nice compromise between
interactivity and repeatability. 

Exploration
+++++++++++

However, it is the serendipitous nature of science that often finding
the right question is not obvious at first. This is certainly true for
astrophysical simulation, especially so for simulations of structure
formation. What are we looking for, and how will we know when we find
it? 

Quite often, the best way forward is to explore the simulation data as
freely as possible.  Without the ability for spot-examination,
serendipitous discovery or general wandering, the code would be simply
a pipeline, rather than a general tool. The flexible extensibility of
yt, that is, the ability to create new derived quantities easily,
as well as the ability to extract and display data regions in a
variety of ways allows for this exploration.

.. _how-yt-thinks-about-data:

Object Methodology
------------------

yt follows a strong object-oriented methodology.  There is no real
global state of yt; all state is contained within objects that
encapsulate an AMR code object or physical region.

Physical Objects vs Code Objects
++++++++++++++++++++++++++++++++

The best way to think about doing things with yt is to think first
of objects. The AMR code puts a number of objects on disk, and yt
has a matching set of objects to mimic these closely as possible. Your
code runs (hopefully) a simulacrum of the physical universe, and thus
in order to make sense of the output data, yt provides a set of
objects meant to mimic the kinds of physical regions and processes you
are interested in. For example, in a simulation of star formation out
of some larger structure (the cosmic dark matter web, a turbulent
molecular cloud), you might be interested in a sphere one parsec in
radius around the point of maximum density. In a simulation of an
accretion disk, you might want a cylindrical region of 1000 AU in
radius and 10 AU in height with its axial vector aligned with the net
angular momentum vector, which may be arbitrary with respect to the
simulation cardinal axes. These are physical objects, and yt has a
set of these too. Finally, you may wish to reduce the data to produce
some essential data that represent a specific process. These
reductions are also objects, and they are included in yt as well.

Somewhat separate from this, but in the same spirit, are plots. In
yt, plots are also objects that one can create, manipulate, and
save. In the case of plots, however, you tell yt what you want to
see, and it can fetch data from the appropriate source. 

In list form,

   Code Objects
     These are things that are on the disk that the AMR code knows about --
     things like grids, data dumps, the grid index and so on.
   Physical Objects
     These are objects like spheres, rectangular prisms, slices, and
     so on. These are collections of related data arranged by physical
     properties, and they are not necessarily associated with a single
     code object.
   Reduced Objects
     These are objects created by taking a set of data and reducing it
     into a smaller format, suitable for a specific purpose.
     Histograms, 1-D profiles, and averages are all members of this
     category.
   Plots
     Plots are somewhat different than other objects, as they are
     neither physical nor code. Instead, the plotting interface
     accepts information about what you want to see, then goes and
     fetches what is necessary--from code, physical, and reduced
     objects as necessary.

.. _intro_to_projections:

Flexible Projections: an Example of Reusable Data Reduction
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

AMR data is best applied when the dynamic range in a quantity of
interest (length or mass scales, typically) is large, but the volume
filling factor of such interesting regions is small. In astronomy,
virtually all available observations are projections on the sky, with
little radial information about the observed structure. In order to
compare with these observations, *projections* are an extremely useful
data reduction for simulations. It is often useful to project to a
given resolution, which may be as high as the highest subdomain in the
AMR data set. However, projecting in a given direction through the
full AMR data set can be quite costly in computing time. yt's
project tool saves an *adaptive* projection when it completes this
costly step, allowing you to make 2D images at whatever resolution you
like with very modest computational resources. This idea, that of
saving as much information as you need (and no more) to make the data
reduction flexible for reuse is another core idea behind yt. You
should not have to spend computational resources and precious time to
replot a projection from a 1000x1000 image to a 2000x2000 image. As a
side note, in this specific case, because the 2D data product yt
produces is "smart", it never needs to use an array in memory as large
as the full effective AMR resolution (which could be very large, and
nearly devoid of unique information).


Derived Fields and Derived Quantities
-------------------------------------

While the primary attraction of yt is the large set of basic code,
physical, reduced, and plot objects already developed, at its core is the fact
that any of the objects can be used as starting points for creating fields and
quantities of your own devices. Derived quantities and derived fields are the
physical objects yt creates from the "primitive" variables the simulation
code stores. These may or may not be the so-called primitive variables of fluid
dynamics (density, velocity, energy): they are whatever your simulation code
writes to disk. 

Derived quantities are those data products derived from these variables such
that the total amount of returned data is *less* than the number of cells.
Derived fields, on the other hand, return a field with *equal* numbers of cells
and the same geometry as the primitive variables from which it was derived. For
example, yt could compute the gravitational potential at every point in
space reconstructed from the density field.

yt already includes a large number of both :ref:`derived fields <field-list>` 
and :ref:`derived quantities <derived-quantities>`, but its real power is 
that it is easy to create your own. See :ref:`creating-derived-fields` for 
detailed instructions on creating derived fields. 
