Introduction
============

History
-------

My name is Matthew Turk, and I am the primary author of yt.  I designed and
implemented it during the course of my graduate studies working with Prof. Tom
Abel at Stanford University, under the auspices of the Department of Energy
through the Stanford Linear Accelerator Center and Los Alamos National Lab.  It
has evolved from a simple data-reader and exporter into what I believe is a
fully-featured toolkit.

yt was designed to be a completely Free (as in beer *and* as in freedom)
user-extensible framework for analyzing and visualizing adaptive mesh
refinement data.  It relies on no proprietary software -- although it can be
and has been extended to interface with proprietary software and libraries --
and has been designed from the ground up to enable users to be as immersed in
the data as they desire.

Originally, yt was written as a frontend to HippoDraw, an extensible and
comprehensive plotting package built at SLAC by Paul Kunz.  Over time, however,
it has been generalized to rely on mostly commodity Python packages, and its
dependencies reduced and ease of installation increased.

What yt is and is not
---------------------

In some sense, yt is also designed to be rather utilitarian.  By virtue of the
fact that it has been written in an interpreted language, it can be somewhat
slower than purely C-based analysis codes, although I believe that to be
mitigated by a cleverness of algorithms and a substantially improved
development time for the user in the case of a desire to expand the
functionality.

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
 * Memory-conserving 1-, 2- and 3-D profiles of arbitrary fields.
 * Halo-finding (HOP) algorithm with full particle information and sphere access

* Data access

 * Arbitrary field definition
 * Derived quantities (average values, spin parameter, bulk velocity, etc)
 * Fast-HDF5 backend for packed and unpacked AMR, NumPy-based HDF4 backend
 * CGS units used everywhere
 * Per-user field and quantity plugins

* Plotting

 * Mathtext TeX-like text formatting
 * Slices, projections, oblique slices
 * Profiles and phase diagrams
 * Linked zooms, colormaps, and saving
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

