# The yt Project

[![Users' Mailing List](https://img.shields.io/badge/Users-List-lightgrey.svg)](https://mail.python.org/mm3/archives/list/yt-users@python.org//)
[![Devel Mailing List](https://img.shields.io/badge/Devel-List-lightgrey.svg)](https://mail.python.org/mm3/archives/list/yt-dev@python.org//)
[![Build Status](https://img.shields.io/travis/yt-project/yt.svg?branch=master)](https://travis-ci.org/yt-project/yt)
[![Latest Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](http://yt-project.org/docs/dev/)
[![Data Hub](https://img.shields.io/badge/data-hub-orange.svg)](https://hub.yt/)
[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](http://numfocus.org)
                
<a href="http://yt-project.org"><img src="doc/source/_static/yt_logo.png" width="300"></a>

yt is an open-source, permissively-licensed python package for analyzing and
visualizing volumetric data.

yt supports structured, variable-resolution meshes, unstructured meshes, and
discrete or sampled data such as particles. Focused on driving
physically-meaningful inquiry, yt has been applied in domains such as
astrophysics, seismology, nuclear engineering, molecular dynamics, and
oceanography. Composed of a friendly community of users and developers, we want
to make it easy to use and develop - we'd love it if you got involved!

We've written a [method
paper](http://adsabs.harvard.edu/abs/2011ApJS..192....9T) you may be interested
in; if you use yt in the preparation of a publication, please consider citing
it.

## Code of Conduct

yt abides by a code of conduct partially modified from the PSF code of conduct,
and is found [in our contributing
guide](http://yt-project.org/docs/dev/developing/developing.html#yt-community-code-of-conduct).

## Installation

If you're using conda with [conda-forge](http://conda-forge.github.io/), you
can install the most recent stable version by running:

```
conda install -c conda-forge yt
```

or by doing:

```
pip install yt
```

If you want the latest nightly build, you can manually install from our
repository:

```
conda install -c http://use.yt/with_conda yt
```

To get set up with a development version, you can clone this repository and
install like this:

```
git clone https://github.com/yt-project/yt yt-git
cd yt-git
pip install -e .
```

To set up yt in a virtualenv (and there are [many good
reasons](https://packaging.python.org/installing/#creating-virtual-environments)
to do so!) you can follow this prescription:

```
# Assuming you have cd'd into yt-git
# It is conventional to create virtualenvs at ~/.virtualenv/
$ mkdir -p ~/.virtualenv
# Assuming your version of Python 3 is 3.4 or higher,
# create a virtualenv named yt
$ python3 -m venv ~/.virtualenv/yt
# Activate it
$ source ~/.virtualenv/yt/bin/activate
# Make sure you run the latest version of pip
$ pip install --upgrade pip
$ pip install -e .
# Output installed packages
$ pip freeze
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

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

(This disclaimer was originally written by
[Adrienne Lowe](https://github.com/adriennefriend) for a
[PyCon talk](https://www.youtube.com/watch?v=6Uj746j9Heo), and was adapted by yt
based on its use in the README file for the
[MetPy project](https://github.com/Unidata/MetPy))

## Resources

We have some community and documentation resources available.

 * Our latest documentation is always at http://yt-project.org/docs/dev/ and it
   includes recipes, tutorials, and API documentation
 * The [discussion mailing
   list](https://mail.python.org/mm3/archives/list/yt-users@python.org//)
   should be your first stop for general questions
 * The [development mailing
   list](https://mail.python.org/mm3/archives/list/yt-dev@python.org//) is
   better suited for more development issues
 * You can also join us on Slack at yt-project.slack.com ([request an
   invite](http://yt-project.org/slack.html))
