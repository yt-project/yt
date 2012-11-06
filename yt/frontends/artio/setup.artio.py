from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

sourcefiles = ['artio_caller.pyx','artio_headers/artio_grid.c','artio_headers/artio_endian.c','artio_headers/artio_mpi.c', 'artio_headers/artio_parameter.c','artio_headers/artio_particle.c','artio_headers/artio_posix.c']

ext_modules = [Extension("artio_caller", 
                          sourcefiles
                          )]

setup(
  name = 'artio header app',
  cmdclass = {'build_ext': build_ext},
#  ext_modules = ext_modules
  ext_modules = [Extension("modartio",sourcefiles)]
)
