# The yt Project

[![Users' Mailing List](https://img.shields.io/badge/Users-List-lightgrey.svg)](http://lists.spacepope.org/listinfo.cgi/yt-users-spacepope.org/)
[![Devel Mailing List](https://img.shields.io/badge/Devel-List-lightgrey.svg)](http://lists.spacepope.org/listinfo.cgi/yt-dev-spacepope.org/)
[![Build Status](https://img.shields.io/travis/yt-project/yt.svg?branch=master)](https://travis-ci.org/yt-project/yt)
[![Latest Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](http://yt-project.org/docs/dev/)
[![Data Hub](https://img.shields.io/badge/data-hub-orange.svg)](https://hub.yt/)
                
[![yt logo](doc/source/_static/yt_logo.png "yt logo")](http://yt-project.org/)

yt is an open-source, permissively-licensed python package for analyzing and
visualizing volumetric data.

yt supports structured, variable-resolution meshes, unstructured meshes, and
discrete or sampled data such as particles. Focused on driving
physically-meaningful inquiry, yt has been applied in domains such as
astrophysics, seismology, nuclear engineering, molecular dynamics, and
oceanography. Composed of a friendly community of users and developers, we want
to make it easy to use and develop — we'd love it if you got involved!

We've written a [method
paper](http://adsabs.harvard.edu/abs/2011ApJS..192....9T) you may be interested
in; if you use yt in the preparation of a publication, please consider citing
it.

## Installation

If you're using conda with [conda-forge](http://conda-forge.github.io/), you
can install the most recent stable version by running:

```
conda install yt
```

Or, if you want the latest nightly build, you can manually install from our
repository:

```
conda install -c http://use.yt/with_conda yt
```

To get set up with a development version, you can clone this repository and
install like this:

```
git clone https://github.com/yt-project/yt
cd yt
python setup.py develop
```

## Getting Started

yt is designed to provide meaningful analysis of data.  We have some Quickstart
example notebooks in the repository:

 * [Introduction](doc/source/quickstart/1\)_Introduction.ipynb)
 * [Data Inspection](doc/source/quickstart/2\)_Data_Inspection.ipynb)
 * [Simple Visualization](doc/source/quickstart/3\)_Simple_Visualization.ipynb)
 * [Data Objects and Time Series](doc/source/quickstart/4\)_Data_Objects_and_Time_Series.ipynb)
 * [Derived Fields and Profiles](doc/source/quickstart/5\)_Derived_Fields_and_Profiles.ipynb)
 * [Volume Rendering](doc/source/quickstart/6\)_Volume_Rendering.ipynb)

If you'd like to try these online, you can visit our [yt Hub](https://hub.yt/)
and run a notebook next to some of our example data.

## Contributing

We love contributions!  yt is open source, built on open source, and we'd love
to have you hang out in our community.

We have developed some [guidelines](CONTRIBUTING.rst) for contributing to yt.

## Resources

We have some community and documentation resources available.

 * Our latest documentation is always at http://yt-project.org/docs/dev/ and it
   includes recipes, tutorials, and API documentation
 * The [discussion mailing
   list](http://lists.spacepope.org/listinfo.cgi/yt-users-spacepope.org/)
   should be your first stop for general questions
 * The [development mailing
   list](http://lists.spacepope.org/listinfo.cgi/yt-dev-spacepope.org/) is
   better suited for more development issues
 * You can also join us on Slack at yt-project.slack.com
