import os
import os.path
import glob
import sys
import time
import subprocess
import shutil
import glob
import distribute_setup
distribute_setup.use_setuptools()

from distutils.command.build_py import build_py
from numpy.distutils.command import build as np_build
from numpy.distutils.misc_util import appendpath
from numpy.distutils import log
from distutils import version

from distutils.core import Command


class BuildForthon(Command):

    """Command for building Forthon modules"""

    description = "Build Forthon modules"
    user_options = tuple()

    def initialize_options(self):

        """init options"""

        pass

    def finalize_options(self):

        """finalize options"""

        pass

    def run(self):

        """runner"""

        cwd = os.getcwd()
        os.chdir(os.path.join(cwd, 'yt/utilities/kdtree'))
        cmd = ["Forthon", "-F", "gfortran", "--compile_first", "fKD_source",
               "--no2underscores", "--fopt", "'-O3'", "fKD",
               "fKD_source.f90"]
        subprocess.check_call(cmd, shell=False)
        shutil.move(glob.glob('build/lib*/fKDpy.so')[0], os.getcwd())
        os.chdir(cwd)

REASON_FILES = []
REASON_DIRS = [
    "",
    "resources",
    "resources/ux",
    "resources/images",
    "resources/css",
    "resources/css/images",
    "app",
    "app/store",
    "app/store/widgets",
    "app/view",
    "app/view/widgets",
    "app/model",
    "app/controller",
    "app/controller/widgets",
    "app/templates",
]

for subdir in REASON_DIRS:
    dir_name = "yt/gui/reason/html/%s/" % (subdir)
    files = []
    for ext in ["js", "html", "css", "png", "ico", "gif"]:
        files += glob.glob("%s/*.%s" % (dir_name, ext))
    REASON_FILES.append((dir_name, files))

# Verify that we have Cython installed
try:
    import Cython
    if version.LooseVersion(Cython.__version__) < version.LooseVersion('0.16'):
        needs_cython = True
    else:
        needs_cython = False
except ImportError as e:
    needs_cython = True

if needs_cython:
    print "Cython is a build-time requirement for the source tree of yt."
    print "Please either install yt from a provided, release tarball,"
    print "or install Cython (version 0.16 or higher)."
    print "You may be able to accomplish this by typing:"
    print "     pip install -U Cython"
    sys.exit(1)

######
# This next bit comes from Matthew Brett, to get Cython working with NumPy
# distutils.  I added a bit to get C++ Cython working.
from os.path import join as pjoin, dirname
from distutils.dep_util import newer_group
from distutils.errors import DistutilsError


def generate_a_pyrex_source(self, base, ext_name, source, extension):
    ''' Monkey patch for numpy build_src.build_src method

    Uses Cython instead of Pyrex.

    Assumes Cython is present
    '''
    if self.inplace:
        target_dir = dirname(base)
    else:
        target_dir = appendpath(self.build_src, dirname(base))
    if extension.language == "c++":
        cplus = True
        file_ext = ".cpp"
    else:
        cplus = False
        file_ext = ".c"
    target_file = pjoin(target_dir, ext_name + file_ext)
    depends = [source] + extension.depends
    if self.force or newer_group(depends, target_file, 'newer'):
        import Cython.Compiler.Main
        log.info("cythonc:> %s" % (target_file))
        self.mkpath(target_dir)
        options = Cython.Compiler.Main.CompilationOptions(
            defaults=Cython.Compiler.Main.default_options,
            include_path=extension.include_dirs,
            language=extension.language, cplus=cplus,
            output_file=target_file)
        cython_result = Cython.Compiler.Main.compile(source,
                                                     options=options)
        if cython_result.num_errors != 0:
            raise DistutilsError("%d errors while compiling %r with Cython"
                                 % (cython_result.num_errors, source))
    return target_file


from numpy.distutils.command import build_src
build_src.build_src.generate_a_pyrex_source = generate_a_pyrex_source
# End snippet
######

import setuptools

VERSION = "2.5dev"

if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')


def get_mercurial_changeset_id(target_dir):
    """adapted from a script by Jason F. Harris, published at

    http://jasonfharris.com/blog/2010/05/versioning-your-application-with-the-mercurial-changeset-hash/

    """
    import subprocess
    import re
    get_changeset = subprocess.Popen('hg identify -b -i',
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     shell=True)

    if (get_changeset.stderr.read() != ""):
        print "Error in obtaining current changeset of the Mercurial repository"
        changeset = None

    changeset = get_changeset.stdout.read().strip()
    if (not re.search("^[0-9a-f]{12}", changeset)):
        print "Current changeset of the Mercurial repository is malformed"
        changeset = None

    return changeset


class my_build(np_build.build):
    def run(self):
        self.run_command("build_forthon")
        np_build.build.run(self)


class my_build_py(build_py):
    def run(self):
        # honor the --dry-run flag
        if not self.dry_run:
            target_dir = os.path.join(self.build_lib, 'yt')
            src_dir = os.getcwd()
            changeset = get_mercurial_changeset_id(src_dir)
            self.mkpath(target_dir)
            with open(os.path.join(target_dir, '__hg_version__.py'), 'w') as fobj:
                fobj.write("hg_version = '%s'\n" % changeset)

            self.run_command("build_forthon")
            build_py.run(self)


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.make_config_py()
    # config.make_svn_version_py()
    config.add_subpackage('yt', 'yt')
    config.add_scripts("scripts/*")

    return config


def setup_package():

    from numpy.distutils.core import setup

    setup(
        name="yt",
        version=VERSION,
        description="An analysis and visualization toolkit for Astrophysical "
                    + "simulations, focusing on Adaptive Mesh Refinement data "
                      "from Enzo, Orion, FLASH, and others.",
        classifiers=["Development Status :: 5 - Production/Stable",
                     "Environment :: Console",
                     "Intended Audience :: Science/Research",
                     "License :: OSI Approved :: GNU General Public License (GPL)",
                     "Operating System :: MacOS :: MacOS X",
                     "Operating System :: POSIX :: AIX",
                     "Operating System :: POSIX :: Linux",
                     "Programming Language :: C",
                     "Programming Language :: Python",
                     "Topic :: Scientific/Engineering :: Astronomy",
                     "Topic :: Scientific/Engineering :: Physics",
                     "Topic :: Scientific/Engineering :: Visualization"],
        keywords='astronomy astrophysics visualization ' +
        'amr adaptivemeshrefinement',
        entry_points={'console_scripts': [
        'yt = yt.utilities.command_line:run_main',
        ],
            'nose.plugins.0.10': [
                'answer-testing = yt.utilities.answer_testing.framework:AnswerTesting'
            ]
        },
        author="Matthew J. Turk",
        author_email="matthewturk@gmail.com",
        url="http://yt-project.org/",
        license="GPL-3",
        configuration=configuration,
        zip_safe=False,
        data_files=REASON_FILES,
        cmdclass={'build_py': my_build_py, 'build_forthon': BuildForthon,
                  'build': my_build},
    )
    return

if __name__ == '__main__':
    setup_package()
