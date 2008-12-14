Introduction
============

History
-------

My name is Matthew Turk, and I am the primary author of yt.  I designed and
implemented it during the course of my graduate studies working with Prof. Tom
Abel at Stanford University, under the auspices of the Department of Energy
through the SLAC National Accelerator Center and, briefly, at Los Alamos
National Lab.  It has evolved from a simple data-reader and exporter into what
I believe is a fully-featured toolkit for analysis and visualization.

yt was designed to be a completely Free (as in beer *and* as in freedom)
user-extensible framework for analyzing and visualizing adaptive mesh
refinement data.  It relies on no proprietary software -- although it can be
and has been extended to interface with proprietary software and libraries --
and has been designed from the ground up to enable users to be as immersed in
the data as they desire.

Originally, yt was written as a frontend to HippoDraw, an extensible and
comprehensive plotting package built at SLAC by Paul Kunz.  Over time
it has been generalized to rely on commodity Python packages (mostly), its
dependencies have been reduced, and thus its installation made significantly easier.

In 2008, yt was extended by Jeff Oishi to include support for the Orion code.
Additionally, it was selected for inclusion in the 1.5 release of the
Enzo code.

What yt is and is not
---------------------

In some sense, yt is also designed to be rather utilitarian.  By virtue of the
fact that it has been written in an interpreted language, it can be somewhat
slower than purely C-based analysis codes, although I believe that to be
mitigated by a cleverness of algorithms and a substantially improved
development time for the user.  Several of the most computatioanlly intensive
problems have been written in C, or rely exclusively on C-based numerical
libraries.

The primary goal has been, and will continue to be, to present an interface to
the user that enables selection and analysis of arbitrary subsets of data.

What functionality does yt offer?
---------------------------------

yt has evolved substantially over the time of its development.  Here is a
non-comprehensive list of features:

* Data Objects

  * Arbitrary data objects (Spheres, cylinders, rectangular prisms, arbitrary index selection)
  * Covering grids (smoothed and raw) for automatic ghost-zone generation
  * Identification of topologically-connected sets in arbitrary fields
  * Projections, orthogonal slices, oblique slices
  * Axially-aligned rays
  * Memory-conserving 1-, 2- and 3-D profiles of arbitrary fields and objects.
  * Halo-finding (HOP) algorithm with full particle information and sphere access
  * `PyTables <http://www.pytables.org/>`_, the data-storage backend for
    marshalling of data.

* Data access

  * Arbitrary field definition
  * Derived quantities (average values, spin parameter, bulk velocity, etc)
  * Custom C- written HDF5 backend for packed and unpacked AMR, NumPy-based HDF4 backend
  * CGS units used everywhere
  * Per-user field and quantity plugins

* Plotting

  * Mathtext TeX-like text formatting
  * Slices, projections, oblique slices
  * Profiles and phase diagrams
  * Linked zooms, colormaps, and saving across multiple plots
  * Contours, vector plots, annotated boxes, grid boundary plot overlays.
  * Simple 3D plotting of phase plots and volume-rendered boxes via hooks into the S2PLOT library

* GUI

  * Linked zooming via slider
  * Interactive re-centering
  * Length scales in human-readable coordinates
  * Drawing of circles for generation of data objects and phase plots
  * Image saving
  * Arbitrary plots within the GUI namespace
  * Full interpreter access to data objects
  * Macros and other scripts able to be run from within the namespace

* Command-line tools

  * Zooming movies
  * Time-series movies
  * Particle movies
  * Radial profile tool

* Access to components

  * Monetary cost: **FREE**.
  * Source code availability: **FULL**.
  * Portability: **YES**.

How do I cite yt?
-----------------

If you use some of the advanced features of yt and would like to cite it in
a publication, you should feel free to cite the 
`Proceedings Paper <http://conference.scipy.org/proceedings/SciPy2008/paper_11>`_ 
with the following BibTeX entry: ::

   @InProceedings{SciPyProceedings_46,
     author =       {Matthew Turk},
     title =        {Analysis and Visualization of Multi-Scale Astrophysical
   Simulations Using Python and NumPy},
     booktitle =   {Proceedings of the 7th Python in Science Conference},
     pages =     {46 - 50},
     address = {Pasadena, CA USA},
     year =      {2008},
     editor =    {Ga\"el Varoquaux and Travis Vaught and Jarrod Millman},
   }
