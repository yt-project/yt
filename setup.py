from distutils.core import setup, Extension
from numarray.numarrayext import NumarrayExtension
import sys, time

if not hasattr(sys, 'version_info') or sys.version_info < (2,2,0,'alpha',0):
    raise SysError

setup(name = "yt",
    version = time.strftime("%y%m%d"),
    description = "A set of classes for manipulating Enzo data",
    url = "http://www.stanford.edu/~mturk/raven.html",
    author="Matthew Turk",
    author_email="mturk@stanford.edu",
    package_dir={"yt":"src"},
    packages = ["yt.enki", "yt.fido", "yt.lagos", "yt", "yt.raven", "yt.deliverator"],
    license="GPL-2",
    #ext_modules=[NumarrayExtension("lagos.RavenCombine",['lagos/point_combine.c'],\
        #include_dirs=["./"],
        #library_dirs=["./"],
        #libraries=['m'])]
    )
