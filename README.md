# The yt Project

[![PyPI](https://img.shields.io/pypi/v/yt)](https://pypi.org/project/yt)
[![Supported Python Versions](https://img.shields.io/pypi/pyversions/yt)](https://pypi.org/project/yt/)
[![Latest Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](http://yt-project.org/docs/dev/)
[![Users' Mailing List](https://img.shields.io/badge/Users-List-lightgrey.svg)](https://mail.python.org/archives/list/yt-users@python.org//)
[![Devel Mailing List](https://img.shields.io/badge/Devel-List-lightgrey.svg)](https://mail.python.org/archives/list/yt-dev@python.org//)
[![Data Hub](https://img.shields.io/badge/data-hub-orange.svg)](https://hub.yt/)
[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](http://numfocus.org)
[![Sponsor our Project](https://img.shields.io/badge/donate-to%20yt-blueviolet)](https://numfocus.salsalabs.org/donate-to-yt/index.html)

<!--- Tests and style --->
![Build and Test](https://github.com/yt-project/yt/workflows/Build%20and%20Test/badge.svg?branch=main)
[![CI (bleeding edge)](https://github.com/yt-project/yt/actions/workflows/bleeding-edge.yaml/badge.svg)](https://github.com/yt-project/yt/actions/workflows/bleeding-edge.yaml)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/yt-project/yt/main.svg)](https://results.pre-commit.ci/latest/github/yt-project/yt/main)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)
<!--- [![codecov](https://codecov.io/gh/yt-project/yt/branch/main/graph/badge.svg)](https://codecov.io/gh/yt-project/yt) --->

<a href="http://yt-project.org"><img src="https://raw.githubusercontent.com/yt-project/yt/main/doc/source/_static/yt_logo.png" width="300"></a>

yt is an open-source, permissively-licensed Python library for analyzing and
visualizing volumetric data.

yt supports structured, variable-resolution meshes, unstructured meshes, and
discrete or sampled data such as particles. Focused on driving
physically-meaningful inquiry, yt has been applied in domains such as
astrophysics, seismology, nuclear engineering, molecular dynamics, and
oceanography. Composed of a friendly community of users and developers, we want
to make it easy to use and develop - we'd love it if you got involved!

We've written a [method
paper](https://ui.adsabs.harvard.edu/abs/2011ApJS..192....9T) you may be interested
in; if you use yt in the preparation of a publication, please consider citing
it.

## Code of Conduct

yt abides by a code of conduct partially modified from the PSF code of conduct,
and is found [in our contributing
guide](http://yt-project.org/docs/dev/developing/developing.html#yt-community-code-of-conduct).

## Installation

You can install the most recent stable version of yt either with conda from
[conda-forge](https://conda-forge.org/):

```shell
conda install -c conda-forge yt
```

or with pip:

```shell
python -m pip install yt
```

More information on the various ways to install yt, and in particular to install from source,
can be found on [the project's website](https://yt-project.org/docs/dev/installing.html).

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
   list](https://mail.python.org/archives/list/yt-users@python.org//)
   should be your first stop for general questions
 * The [development mailing
   list](https://mail.python.org/archives/list/yt-dev@python.org//) is
   better suited for more development issues
 * You can also join us on Slack at yt-project.slack.com ([request an
   invite](https://yt-project.org/slack.html))

Is your code compatible with yt ? Great ! Please consider giving us a shoutout as a shiny badge in your README

- markdown
```markdown
[![yt-project](https://img.shields.io/static/v1?label="works%20with"&message="yt"&color="blueviolet")](https://yt-project.org)
```
- rst
```reStructuredText
|yt-project|

.. |yt-project| image:: https://img.shields.io/static/v1?label="works%20with"&message="yt"&color="blueviolet"
   :target: https://yt-project.org
```

## Powered by NumFOCUS

yt is a fiscally sponsored project of [NumFOCUS](https://numfocus.org/).
If you're interested in
supporting the active maintenance and development of this project, consider
[donating to the project](https://numfocus.salsalabs.org/donate-to-yt/index.html).
