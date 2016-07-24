"""
A means of running standalone commands with a shared set of options.



"""
from __future__ import print_function
from __future__ import absolute_import

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import argparse
import base64
import getpass
import numpy as np
import os
import sys
import subprocess
import tempfile
import json
import pprint

from yt.config import ytcfg
ytcfg["yt","__command_line"] = "True"
from yt.startup_tasks import parser, subparsers
from yt.funcs import \
    ensure_list, \
    get_hg_version, \
    get_yt_version, \
    mylog, \
    ensure_dir_exists, \
    update_hg, \
    enable_plugins
from yt.extern.six import add_metaclass, string_types
from yt.extern.six.moves import urllib, input
from yt.convenience import load
from yt.visualization.plot_window import \
    SlicePlot, \
    ProjectionPlot
from yt.utilities.metadata import get_metadata
from yt.utilities.exceptions import \
    YTOutputNotIdentified, YTFieldNotParseable

# loading field plugins for backward compatibility, since this module
# used to do "from yt.mods import *"
if ytcfg.getboolean("yt","loadfieldplugins"):
    enable_plugins()

_default_colormap = ytcfg.get("yt", "default_colormap")

def _fix_ds(arg):
    if os.path.isdir("%s" % arg) and \
        os.path.exists("%s/%s" % (arg,arg)):
        ds = load("%s/%s" % (arg,arg))
    elif os.path.isdir("%s.dir" % arg) and \
        os.path.exists("%s.dir/%s" % (arg,arg)):
        ds = load("%s.dir/%s" % (arg,arg))
    elif arg.endswith(".index"):
        ds = load(arg[:-10])
    else:
        ds = load(arg)
    return ds

def _add_arg(sc, arg):
    if isinstance(arg, string_types):
        arg = _common_options[arg].copy()
    argc = dict(arg.items())
    argnames = []
    if "short" in argc: argnames.append(argc.pop('short'))
    if "longname" in argc: argnames.append(argc.pop('longname'))
    sc.add_argument(*argnames, **argc)

def _print_failed_source_update(reinstall=False):
    print()
    print("The yt package is not installed from a mercurial repository,")
    print("so you must update this installation manually.")
    if 'Continuum Analytics' in sys.version or 'Anaconda' in sys.version:
        # see http://stackoverflow.com/a/21318941/1382869 for why we need
        # to check both Continuum *and* Anaconda
        print()
        print("Since it looks like you are using a python installation")
        print("that is managed by conda, you may want to do:")
        print()
        print("    $ conda update yt")
        print()
        print("to update your yt installation.")
        if reinstall is True:
            print()
            print("To update all of your packages, you can do:")
            print()
            print("    $ conda update --all")

def _print_installation_information(path):
    import yt
    print()
    print("yt module located at:")
    print("    %s" % (path))
    if "YT_DEST" in os.environ:
        spath = os.path.join(
            os.environ["YT_DEST"], "src", "yt-supplemental")
        if os.path.isdir(spath):
            print("The supplemental repositories are located at:")
            print("    %s" % (spath))
    print()
    print("The current version of yt is:")
    print()
    print("---")
    print("Version = %s" % yt.__version__)
    vstring = get_hg_version(path)
    if vstring is not None:
        print("Changeset = %s" % vstring.strip().decode("utf-8"))
    print("---")
    return vstring

class YTCommandSubtype(type):
    def __init__(cls, name, b, d):
        type.__init__(cls, name, b, d)
        if cls.name is not None:
            names = ensure_list(cls.name)
            for name in names:
                sc = subparsers.add_parser(name,
                    description = cls.description,
                    help = cls.description)
                sc.set_defaults(func=cls.run)
                for arg in cls.args:
                    _add_arg(sc, arg)

@add_metaclass(YTCommandSubtype)
class YTCommand(object):
    args = ()
    name = None
    description = ""
    aliases = ()
    ndatasets = 1

    @classmethod
    def run(cls, args):
        self = cls()
        # Check for some things we know; for instance, comma separated
        # field names should be parsed as tuples.
        if getattr(args, 'field', None) is not None and ',' in args.field:
            if args.field.count(",") > 1:
                raise YTFieldNotParseable(args.field)
            args.field = tuple(_.strip() for _ in args.field.split(","))
        if getattr(args, 'weight', None) is not None and ',' in args.weight:
            if args.weight.count(",") > 1:
                raise YTFieldNotParseable(args.weight)
            args.weight = tuple(_.strip() for _ in args.weight.split(","))
        # Some commands need to be run repeatedly on datasets
        # In fact, this is the rule and the opposite is the exception
        # BUT, we only want to parse the arguments once.
        if cls.ndatasets > 1:
            self(args)
        else:
            ds_args = getattr(args, "ds", [])
            if len(ds_args) > 1:
                datasets = args.ds
                for ds in datasets:
                    args.ds = ds
                    self(args)
            elif len(ds_args) == 0:
                datasets = []
                self(args)
            else:
                args.ds = getattr(args, 'ds', [None])[0]
                self(args)

class GetParameterFiles(argparse.Action):
    def __call__(self, parser, namespace, values, option_string = None):
        if len(values) == 1:
            datasets = values
        elif len(values) == 2 and namespace.basename is not None:
            datasets = ["%s%04i" % (namespace.basename, r)
                   for r in range(int(values[0]), int(values[1]), namespace.skip) ]
        else:
            datasets = values
        namespace.ds = [_fix_ds(ds) for ds in datasets]

_common_options = dict(
    all     = dict(longname="--all", dest="reinstall",
                   default=False, action="store_true",
                   help="Reinstall the full yt stack in the current location."),
    ds      = dict(short="ds", action=GetParameterFiles,
                   nargs="+", help="datasets to run on"),
    ods     = dict(action=GetParameterFiles, dest="ds",
                   nargs="*", help="(Optional) datasets to run on"),
    axis    = dict(short="-a", longname="--axis",
                   action="store", type=int,
                   dest="axis", default=4,
                   help="Axis (4 for all three)"),
    log     = dict(short="-l", longname="--log",
                   action="store_true",
                   dest="takelog", default=True,
                   help="Use logarithmic scale for image"),
    linear  = dict(longname="--linear",
                   action="store_false",
                   dest="takelog",
                   help="Use linear scale for image"),
    text    = dict(short="-t", longname="--text",
                   action="store", type=str,
                   dest="text", default=None,
                   help="Textual annotation"),
    field   = dict(short="-f", longname="--field",
                   action="store", type=str,
                   dest="field", default="density",
                   help=("Field to color by, "
                         "use a comma to separate field tuple values")),
    weight  = dict(short="-g", longname="--weight",
                   action="store", type=str,
                   dest="weight", default=None,
                   help=("Field to weight projections with, "
                         "use a comma to separate field tuple values")),
    cmap    = dict(longname="--colormap",
                   action="store", type=str,
                   dest="cmap", default=_default_colormap,
                   help="Colormap name"),
    zlim    = dict(short="-z", longname="--zlim",
                   action="store", type=float,
                   dest="zlim", default=None,
                   nargs=2,
                   help="Color limits (min, max)"),
    dex     = dict(longname="--dex",
                   action="store", type=float,
                   dest="dex", default=None,
                   nargs=1,
                   help="Number of dex above min to display"),
    width   = dict(short="-w", longname="--width",
                   action="store", type=float,
                   dest="width", default=None,
                   help="Width in specified units"),
    unit    = dict(short="-u", longname="--unit",
                   action="store", type=str,
                   dest="unit", default='1',
                   help="Desired units"),
    center  = dict(short="-c", longname="--center",
                   action="store", type=float,
                   dest="center", default=None,
                   nargs=3,
                   help="Center, space separated (-1 -1 -1 for max)"),
    max     = dict(short="-m", longname="--max",
                   action="store_true",
                   dest="max",default=False,
                   help="Center the plot on the density maximum"),
    bn      = dict(short="-b", longname="--basename",
                   action="store", type=str,
                   dest="basename", default=None,
                   help="Basename of datasets"),
    output  = dict(short="-o", longname="--output",
                   action="store", type=str,
                   dest="output", default="frames/",
                   help="Folder in which to place output images"),
    outputfn= dict(short="-o", longname="--output",
                   action="store", type=str,
                   dest="output", default=None,
                   help="File in which to place output"),
    skip    = dict(short="-s", longname="--skip",
                   action="store", type=int,
                   dest="skip", default=1,
                   help="Skip factor for outputs"),
    proj    = dict(short="-p", longname="--projection",
                   action="store_true",
                   dest="projection", default=False,
                   help="Use a projection rather than a slice"),
    maxw    = dict(longname="--max-width",
                   action="store", type=float,
                   dest="max_width", default=1.0,
                   help="Maximum width in code units"),
    minw    = dict(longname="--min-width",
                   action="store", type=float,
                   dest="min_width", default=50,
                   help="Minimum width in units of smallest dx (default: 50)"),
    nframes = dict(short="-n", longname="--nframes",
                   action="store", type=int,
                   dest="nframes", default=100,
                   help="Number of frames to generate"),
    slabw   = dict(longname="--slab-width",
                   action="store", type=float,
                   dest="slab_width", default=1.0,
                   help="Slab width in specified units"),
    slabu   = dict(short="-g", longname="--slab-unit",
                   action="store", type=str,
                   dest="slab_unit", default='1',
                   help="Desired units for the slab"),
    ptype   = dict(longname="--particle-type",
                   action="store", type=int,
                   dest="ptype", default=2,
                   help="Particle type to select"),
    agecut  = dict(longname="--age-cut",
                   action="store", type=float,
                   dest="age_filter", default=None,
                   nargs=2,
                   help="Bounds for the field to select"),
    uboxes  = dict(longname="--unit-boxes",
                   action="store_true",
                   dest="unit_boxes",
                   help="Display heldsul unit boxes"),
    thresh  = dict(longname="--threshold",
                   action="store", type=float,
                   dest="threshold", default=None,
                   help="Density threshold"),
    dm_only = dict(longname="--all-particles",
                   action="store_false",
                   dest="dm_only", default=True,
                   help="Use all particles"),
    grids   = dict(longname="--show-grids",
                   action="store_true",
                   dest="grids", default=False,
                   help="Show the grid boundaries"),
    time    = dict(longname="--time",
                   action="store_true",
                   dest="time", default=False,
                   help="Print time in years on image"),
    contours    = dict(longname="--contours",
                   action="store",type=int,
                   dest="contours", default=None,
                   help="Number of Contours for Rendering"),
    contour_width  = dict(longname="--contour_width",
                   action="store",type=float,
                   dest="contour_width", default=None,
                   help="Width of gaussians used for rendering."),
    enhance   = dict(longname="--enhance",
                   action="store_true",
                   dest="enhance", default=False,
                   help="Enhance!"),
    valrange  = dict(short="-r", longname="--range",
                   action="store", type=float,
                   dest="valrange", default=None,
                   nargs=2,
                   help="Range, space separated"),
    up  = dict(longname="--up",
                   action="store", type=float,
                   dest="up", default=None,
                   nargs=3,
                   help="Up, space separated"),
    viewpoint  = dict(longname="--viewpoint",
                   action="store", type=float,
                   dest="viewpoint", default=[1., 1., 1.],
                   nargs=3,
                   help="Viewpoint, space separated"),
    pixels    = dict(longname="--pixels",
                   action="store",type=int,
                   dest="pixels", default=None,
                   help="Number of Pixels for Rendering"),
    halos   = dict(longname="--halos",
                   action="store", type=str,
                   dest="halos",default="multiple",
                   help="Run halo profiler on a 'single' halo or 'multiple' halos."),
    halo_radius = dict(longname="--halo_radius",
                       action="store", type=float,
                       dest="halo_radius",default=0.1,
                       help="Constant radius for profiling halos if using hop output files with no radius entry. Default: 0.1."),
    halo_radius_units = dict(longname="--halo_radius_units",
                             action="store", type=str,
                             dest="halo_radius_units",default="1",
                             help="Units for radius used with --halo_radius flag. Default: '1' (code units)."),
    halo_hop_style = dict(longname="--halo_hop_style",
                          action="store", type=str,
                          dest="halo_hop_style",default="new",
                          help="Style of hop output file.  'new' for yt_hop files and 'old' for enzo_hop files."),
    halo_dataset = dict(longname="--halo_dataset",
                               action="store", type=str,
                               dest="halo_dataset",default=None,
                               help="HaloProfiler dataset."),
    make_profiles = dict(longname="--make_profiles",
                         action="store_true", default=False,
                         help="Make profiles with halo profiler."),
    make_projections = dict(longname="--make_projections",
                            action="store_true", default=False,
                            help="Make projections with halo profiler.")

    )

def _get_yt_stack_date():
    if "YT_DEST" not in os.environ:
        print("Could not determine when yt stack was last updated.")
        return
    date_file = os.path.join(os.environ["YT_DEST"], ".yt_update")
    if not os.path.exists(date_file):
        print("Could not determine when yt stack was last updated.")
        return
    print("".join(open(date_file, 'r').readlines()))
    print("To update all dependencies, run \"yt update --all\".")

def _update_yt_stack(path):
    "Rerun the install script to updated all dependencies."

    if "YT_DEST" not in os.environ:
        print()
        print("This yt installation does not appear to be managed by the")
        print("source-based install script, but 'update --all' was specified.")
        print("You will need to update your dependencies manually.")
        return

    install_script = os.path.join(path, "doc/install_script.sh")
    if not os.path.exists(install_script):
        print()
        print("Install script not found!")
        print("The install script should be here: %s," % install_script)
        print("but it was not.")
        return

    print()
    print("We will now attempt to update the yt stack located at:")
    print("    %s." % os.environ["YT_DEST"])
    print()
    print("[hit enter to continue or Ctrl-C to stop]")
    try:
        input()
    except:
        sys.exit(0)
    os.environ["REINST_YT"] = "1"
    ret = subprocess.call(["bash", install_script])
    print()
    if ret:
        print("The install script seems to have failed.")
        print("Check the output above.")
    else:
        print("The yt stack has been updated successfully.")
        print("Now get back to work!")

# This code snippet is modified from Georg Brandl
def bb_apicall(endpoint, data, use_pass = True):
    uri = 'https://api.bitbucket.org/1.0/%s/' % endpoint
    # since bitbucket doesn't return the required WWW-Authenticate header when
    # making a request without Authorization, we cannot use the standard urllib2
    # auth handlers; we have to add the requisite header from the start
    if data is not None:
        data = urllib.parse.urlencode(data)
    req = urllib.request.Request(uri, data)
    if use_pass:
        username = input("Bitbucket Username? ")
        password = getpass.getpass()
        upw = '%s:%s' % (username, password)
        req.add_header('Authorization', 'Basic %s' % base64.b64encode(upw).strip())
    return urllib.request.urlopen(req).read()

class YTBugreportCmd(YTCommand):
    name = "bugreport"
    description = \
        """
        Report a bug in yt

        """

    def __call__(self, args):
        print("===============================================================")
        print()
        print("Hi there!  Welcome to the yt bugreport taker.")
        print()
        print("===============================================================")
        print("At any time in advance of the upload of the bug, you should feel free")
        print("to ctrl-C out and submit the bug report manually by going here:")
        print("   http://bitbucket.org/yt_analysis/yt/issues/new")
        print()
        print("Also, in order to submit a bug through this interface, you")
        print("need a Bitbucket account. If you don't have one, exit this ")
        print("bugreport now and run the 'yt bootstrap_dev' command to create one.")
        print()
        print("Have you checked the existing bug reports to make")
        print("sure your bug has not already been recorded by someone else?")
        print("   http://bitbucket.org/yt_analysis/yt/issues?status=new&status=open")
        print()
        print("Finally, are you sure that your bug is, in fact, a bug? It might")
        print("simply be a misunderstanding that could be cleared up by")
        print("visiting the yt irc channel or getting advice on the email list:")
        print("   http://yt-project.org/irc.html")
        print("   http://lists.spacepope.org/listinfo.cgi/yt-users-spacepope.org")
        print()
        summary = input("Press <enter> if you remain firm in your conviction to continue.")
        print()
        print()
        print("Okay, sorry about that. How about a nice, pithy ( < 12 words )")
        print("summary of the bug?  (e.g. 'Particle overlay problem with parallel ")
        print("projections')")
        print()
        try:
            current_version = get_yt_version()
        except:
            current_version = "Unavailable"
        summary = input("Summary? ")
        bugtype = "bug"
        data = dict(title = summary, type=bugtype)
        print()
        print("Okay, now let's get a bit more information.")
        print()
        print("Remember that if you want to submit a traceback, you can run")
        print("any script with --paste or --detailed-paste to submit it to")
        print("the pastebin and then include the link in this bugreport.")
        if "EDITOR" in os.environ:
            print()
            print("Press enter to spawn your editor, %s" % os.environ["EDITOR"])
            input()
            tf = tempfile.NamedTemporaryFile(delete=False)
            fn = tf.name
            tf.close()
            subprocess.call("$EDITOR %s" % fn, shell = True)
            content = open(fn).read()
            try:
                os.unlink(fn)
            except:
                pass
        else:
            print()
            print("Couldn't find an $EDITOR variable.  So, let's just take")
            print("take input here.  Type up your summary until you're ready")
            print("to be done, and to signal you're done, type --- by itself")
            print("on a line to signal your completion.")
            print()
            print("(okay, type now)")
            print()
            lines = []
            while 1:
                line = input()
                if line.strip() == "---": break
                lines.append(line)
            content = "\n".join(lines)
        content = "Reporting Version: %s\n\n%s" % (current_version, content)
        endpoint = "repositories/yt_analysis/yt/issues"
        data['content'] = content
        print()
        print("===============================================================")
        print()
        print("Okay, we're going to submit with this:")
        print()
        print("Summary: %s" % (data['title']))
        print()
        print("---")
        print(content)
        print("---")
        print()
        print("===============================================================")
        print()
        print("Is that okay?  If not, hit ctrl-c.  Otherwise, enter means")
        print("'submit'.  Next we'll ask for your Bitbucket Username.")
        print("If you don't have one, run the 'yt bootstrap_dev' command.")
        print()
        input()
        retval = bb_apicall(endpoint, data, use_pass=True)
        import json
        retval = json.loads(retval)
        url = "http://bitbucket.org/yt_analysis/yt/issue/%s" % retval['local_id']
        print()
        print("===============================================================")
        print()
        print("Thanks for your bug report!  Together we'll make yt totally bug free!")
        print("You can view bug report here:")
        print("   %s" % url)
        print()
        print("Keep in touch!")
        print()


class YTHubRegisterCmd(YTCommand):
    name = "hub_register"
    description = \
        """
        Register a user on the Hub: http://hub.yt-project.org/
        """
    def __call__(self, args):
        # We need these pieces of information:
        #   1. Name
        #   2. Email
        #   3. Username
        #   4. Password (and password2)
        #   5. (optional) URL
        #   6. "Secret" key to make it epsilon harder for spammers
        if ytcfg.get("yt","hub_api_key") != "":
            print("You seem to already have an API key for the hub in")
            print("~/.yt/config .  Delete this if you want to force a")
            print("new user registration.")
        print("Awesome!  Let's start by registering a new user for you.")
        print("Here's the URL, for reference: http://hub.yt-project.org/ ")
        print()
        print("As always, bail out with Ctrl-C at any time.")
        print()
        print("What username would you like to go by?")
        print()
        username = input("Username? ")
        if len(username) == 0: sys.exit(1)
        print()
        print("To start out, what's your name?")
        print()
        name = input("Name? ")
        if len(name) == 0: sys.exit(1)
        print()
        print("And your email address?")
        print()
        email = input("Email? ")
        if len(email) == 0: sys.exit(1)
        print()
        print("Please choose a password:")
        print()
        while 1:
            password1 = getpass.getpass("Password? ")
            password2 = getpass.getpass("Confirm? ")
            if len(password1) == 0: continue
            if password1 == password2: break
            print("Sorry, they didn't match!  Let's try again.")
            print()
        print()
        print("Would you like a URL displayed for your user?")
        print("Leave blank if no.")
        print()
        url = input("URL? ")
        print()
        print("Okay, press enter to register.  You should receive a welcome")
        print("message at %s when this is complete." % email)
        print()
        input()
        data = dict(name = name, email = email, username = username,
                    password = password1, password2 = password2,
                    url = url, zap = "rowsdower")
        data = urllib.parse.urlencode(data)
        hub_url = "https://hub.yt-project.org/create_user"
        req = urllib.request.Request(hub_url, data)
        try:
            urllib.request.urlopen(req).read()
        except urllib.error.HTTPError as exc:
            if exc.code == 400:
                print("Sorry, the Hub couldn't create your user.")
                print("You can't register duplicate users, which is the most")
                print("common cause of this error.  All values for username,")
                print("name, and email must be unique in our system.")
                sys.exit(1)
        except urllib.URLError as exc:
            print("Something has gone wrong.  Here's the error message.")
            raise exc
        print()
        print("SUCCESS!")
        print()


class YTInstInfoCmd(YTCommand):
    name = ["instinfo", "version"]
    args = (
            dict(short="-u", longname="--update-source", action="store_true",
                 default = False,
                 help="Update the yt installation, if able"),
            dict(short="-o", longname="--output-version", action="store",
                  default = None, dest="outputfile",
                  help="File into which the current revision number will be" +
                       "stored")
           )
    description = \
        """
        Get some information about the yt installation

        """

    def __call__(self, opts):
        import pkg_resources
        yt_provider = pkg_resources.get_provider("yt")
        path = os.path.dirname(yt_provider.module_path)
        vstring = _print_installation_information(path)
        if vstring is not None:
            print("This installation CAN be automatically updated.")
            if opts.update_source:
                update_hg(path)
                print("Updated successfully.")
                _get_yt_stack_date()
        elif opts.update_source:
            _print_failed_source_update()
        if vstring is not None and opts.outputfile is not None:
            open(opts.outputfile, "w").write(vstring)

class YTLoadCmd(YTCommand):
    name = "load"
    description = \
        """
        Load a single dataset into an IPython instance

        """

    args = ("ds", )

    def __call__(self, args):
        if args.ds is None:
            print("Could not load file.")
            sys.exit()
        import yt.mods
        import yt

        import IPython
        from distutils import version
        if version.LooseVersion(IPython.__version__) <= version.LooseVersion('0.10'):
            api_version = '0.10'
        else:
            api_version = '0.11'

        local_ns = yt.mods.__dict__.copy()
        local_ns['ds'] = args.ds
        local_ns['pf'] = args.ds
        local_ns['yt'] = yt

        if api_version == '0.10':
            shell = IPython.Shell.IPShellEmbed()
            shell(local_ns = local_ns,
                  header =
                  "\nHi there!  Welcome to yt.\n\nWe've loaded your dataset as 'ds'.  Enjoy!"
                  )
        else:
            try:
                from traitlets.config.loader import Config
            except ImportError:
                from IPython.config.loader import Config
            import sys
            cfg = Config()
            # prepend sys.path with current working directory
            sys.path.insert(0,'')
            IPython.embed(config=cfg,user_ns=local_ns)

class YTMapserverCmd(YTCommand):
    args = ("proj", "field", "weight",
            dict(short="-a", longname="--axis", action="store", type=int,
                 dest="axis", default=0, help="Axis (4 for all three)"),
            dict(short ="-o", longname="--host", action="store", type=str,
                   dest="host", default=None, help="IP Address to bind on"),
            "ds",
            )

    name = "mapserver"
    description = \
        """
        Serve a plot in a GMaps-style interface

        """

    def __call__(self, args):
        if sys.version_info >= (3,0,0):
            print("yt mapserver is disabled for Python 3.")
            return -1
        ds = args.ds
        if args.axis == 4:
            print("Doesn't work with multiple axes!")
            return
        if args.projection:
            p = ProjectionPlot(ds, args.axis, args.field, weight_field=args.weight)
        else:
            p = SlicePlot(ds, args.axis, args.field)
        from yt.visualization.mapserver.pannable_map import PannableMapServer
        PannableMapServer(p.data_source, args.field)
        import yt.extern.bottle as bottle
        bottle.debug(True)
        bottle_dir = os.path.dirname(bottle.__file__)
        sys.path.append(bottle_dir)
        if args.host is not None:
            colonpl = args.host.find(":")
            if colonpl >= 0:
                port = int(args.host.split(":")[-1])
                args.host = args.host[:colonpl]
            else:
                port = 8080
            bottle.run(server='rocket', host=args.host, port=port)
        else:
            bottle.run(server='rocket')
        sys.path.remove(bottle_dir)


class YTPastebinCmd(YTCommand):
    name = "pastebin"
    args = (
             dict(short="-l", longname="--language", action="store",
                  default = None, dest="language",
                  help="Use syntax highlighter for the file in language"),
             dict(short="-L", longname="--languages", action="store_true",
                  default = False, dest="languages",
                  help="Retrive a list of supported languages"),
             dict(short="-e", longname="--encoding", action="store",
                  default = 'utf-8', dest="encoding",
                  help="Specify the encoding of a file (default is "
                        "utf-8 or guessing if available)"),
             dict(short="-b", longname="--open-browser", action="store_true",
                  default = False, dest="open_browser",
                  help="Open the paste in a web browser"),
             dict(short="-p", longname="--private", action="store_true",
                  default = False, dest="private",
                  help="Paste as private"),
             dict(short="-c", longname="--clipboard", action="store_true",
                  default = False, dest="clipboard",
                  help="File to output to; else, print."),
             dict(short="file", type=str),
            )
    description = \
        """
        Post a script to an anonymous pastebin

        """

    def __call__(self, args):
        import yt.utilities.lodgeit as lo
        lo.main(args.file, languages=args.languages, language=args.language,
                 encoding=args.encoding, open_browser=args.open_browser,
                 private=args.private, clipboard=args.clipboard)

class YTPastebinGrabCmd(YTCommand):
    args = (dict(short="number", type=str),)
    name = "pastebin_grab"
    description = \
        """
        Print an online pastebin to STDOUT for local use.
        """

    def __call__(self, args):
        import yt.utilities.lodgeit as lo
        lo.main( None, download=args.number )

class YTNotebookUploadCmd(YTCommand):
    args = (dict(short="file", type=str),)
    description = \
        """
        Upload an IPython notebook to hub.yt-project.org.
        """

    name = "upload_notebook"
    def __call__(self, args):
        filename = args.file
        if not os.path.isfile(filename):
            raise IOError(filename)
        if not filename.endswith(".ipynb"):
            print("File must be an IPython notebook!")
            return 1
        import json
        try:
            t = json.loads(open(filename).read())['metadata']['name']
        except (ValueError, KeyError):
            print("File does not appear to be an IPython notebook.")
        if len(t) == 0:
            t = filename.strip(".ipynb")
        from yt.utilities.minimal_representation import MinimalNotebook
        mn = MinimalNotebook(filename, t)
        rv = mn.upload()
        print("Upload successful!")
        print()
        print("To access your raw notebook go here:")
        print()
        print("  %s" % (rv['url']))
        print()
        print("To view your notebook go here:")
        print()
        print("  %s" % (rv['url'].replace("/go/", "/nb/")))
        print()

class YTPlotCmd(YTCommand):
    args = ("width", "unit", "bn", "proj", "center", "zlim", "axis", "field",
            "weight", "skip", "cmap", "output", "grids", "time", "ds", "max",
            "log", "linear")

    name = "plot"

    description = \
        """
        Create a set of images

        """

    def __call__(self, args):
        ds = args.ds
        center = args.center
        if args.center == (-1,-1,-1):
            mylog.info("No center fed in; seeking.")
            v, center = ds.find_max("density")
        if args.max:
            v, center = ds.find_max("density")
        elif args.center is None:
            center = 0.5*(ds.domain_left_edge + ds.domain_right_edge)
        center = np.array(center)
        if ds.dimensionality < 3:
            dummy_dimensions = np.nonzero(ds.index.grids[0].ActiveDimensions <= 1)
            axes = ensure_list(dummy_dimensions[0][0])
        elif args.axis == 4:
            axes = range(3)
        else:
            axes = [args.axis]

        unit = args.unit
        if unit is None:
            unit = 'unitary'
        if args.width is None:
            width = None
        else:
            width = (args.width, args.unit)

        for ax in axes:
            mylog.info("Adding plot for axis %i", ax)
            if args.projection:
                plt = ProjectionPlot(ds, ax, args.field, center=center,
                                     width=width,
                                     weight_field=args.weight)
            else:
                plt = SlicePlot(ds, ax, args.field, center=center,
                                width=width)
            if args.grids:
                plt.annotate_grids()
            if args.time:
                time = ds.current_time.in_units("yr")
                plt.annotate_text((0.2,0.8), 't = %5.2e yr'%time)

            plt.set_cmap(args.field, args.cmap)
            plt.set_log(args.field, args.takelog)
            if args.zlim:
                plt.set_zlim(args.field,*args.zlim)
            ensure_dir_exists(args.output)
            plt.save(os.path.join(args.output,"%s" % (ds)))


class YTRPDBCmd(YTCommand):
    name = "rpdb"
    description = \
        """
        Connect to a currently running (on localhost) rpd session.

        Commands run with --rpdb will trigger an rpdb session with any
        uncaught exceptions.

        """
    args = (
            dict(short="-t", longname="--task", action="store",
                 default = 0, dest='task',
                 help="Open a web browser."),
           )

    def __call__(self, args):
        from . import rpdb
        rpdb.run_rpdb(int(args.task))

class YTNotebookCmd(YTCommand):
    name = ["notebook"]
    args = (
            dict(short="-o", longname="--open-browser", action="store_true",
                 default = False, dest='open_browser',
                 help="Open a web browser."),
            dict(short="-p", longname="--port", action="store",
                 default = 0, dest='port',
                 help="Port to listen on; defaults to auto-detection."),
            dict(short="-prof", longname="--profile", action="store",
                 default = None, dest="profile",
                 help="The IPython profile to use when lauching the kernel."),
            dict(short="-n", longname="--no-password", action="store_true",
                 default = False, dest='no_password',
                 help="If set, do not prompt or use a password."),
            )
    description = \
        """
        Run the IPython Notebook
        """
    def __call__(self, args):
        kwargs = {}
        try:
            # IPython 1.0+
            from IPython.html.notebookapp import NotebookApp
        except ImportError:
            # pre-IPython v1.0
            from IPython.frontend.html.notebook.notebookapp import NotebookApp
        print("You must choose a password so that others cannot connect to " \
              "your notebook.")
        pw = ytcfg.get("yt", "notebook_password")
        if len(pw) == 0 and not args.no_password:
            import IPython.lib
            pw = IPython.lib.passwd()
            print("If you would like to use this password in the future,")
            print("place a line like this inside the [yt] section in your")
            print("yt configuration file at ~/.yt/config")
            print()
            print("notebook_password = %s" % pw)
            print()
        elif args.no_password:
            pw = None
        if args.port != 0:
            kwargs['port'] = int(args.port)
        if args.profile is not None:
            kwargs['profile'] = args.profile
        if pw is not None:
            kwargs['password'] = pw
        app = NotebookApp(open_browser=args.open_browser,
                          **kwargs)
        app.initialize(argv=[])
        print()
        print("***************************************************************")
        print()
        print("The notebook is now live at:")
        print()
        print("     http://127.0.0.1:%s/" % app.port)
        print()
        print("Recall you can create a new SSH tunnel dynamically by pressing")
        print("~C and then typing -L%s:localhost:%s" % (app.port, app.port))
        print("where the first number is the port on your local machine. ")
        print()
        print("If you are using %s on your machine already, try " \
              "-L8889:localhost:%s" % (app.port, app.port))
        print()
        print("***************************************************************")
        print()
        app.start()



class YTStatsCmd(YTCommand):
    args = ('outputfn','bn','skip','ds','field',
            dict(longname="--max", action='store_true', default=False,
                 dest='max', help="Display maximum of field requested through -f option."),
            dict(longname="--min", action='store_true', default=False,
                 dest='min', help="Display minimum of field requested through -f option."))
    name = "stats"
    description = \
        """
        Print stats and max/min value of a given field (if requested),
        for one or more datasets

        (default field is density)

        """

    def __call__(self, args):
        ds = args.ds
        ds.print_stats()
        vals = {}
        if args.field in ds.derived_field_list:
            if args.max is True:
                vals['min'] = ds.find_max(args.field)
                print("Maximum %s: %0.5e at %s" % (args.field,
                    vals['min'][0], vals['min'][1]))
            if args.min is True:
                vals['max'] = ds.find_min(args.field)
                print("Minimum %s: %0.5e at %s" % (args.field,
                    vals['max'][0], vals['max'][1]))
        if args.output is not None:
            t = ds.current_time * ds['years']
            with open(args.output, "a") as f:
                f.write("%s (%0.5e years)\n" % (ds, t))
                if 'min' in vals:
                    f.write('Minimum %s is %0.5e at %s\n' % (
                        args.field, vals['min'][0], vals['min'][1]))
                if 'max' in vals:
                    f.write('Maximum %s is %0.5e at %s\n' % (
                        args.field, vals['max'][0], vals['max'][1]))

class YTUpdateCmd(YTCommand):
    args = ("all", )
    name = "update"
    description = \
        """
        Update the yt installation to the most recent version

        """

    def __call__(self, opts):
        import pkg_resources
        yt_provider = pkg_resources.get_provider("yt")
        path = os.path.dirname(yt_provider.module_path)
        vstring = _print_installation_information(path)
        if vstring is not None:
            print()
            print("This installation CAN be automatically updated.")
            update_hg(path, skip_rebuild=opts.reinstall)
            print("Updated successfully.")
            _get_yt_stack_date()
            if opts.reinstall:
                _update_yt_stack(path)
        else:
            _print_failed_source_update(opts.reinstall)


class YTDeleteImageCmd(YTCommand):
    args = (dict(short="delete_hash", type=str),)
    description = \
        """
        Delete image from imgur.com.

        """
    name = "delete_image"
    def __call__(self, args):
        headers = {'Authorization':
            'Client-ID {}'.format(ytcfg.get("yt", "imagebin_api_key"))}

        delete_url = ytcfg.get("yt", "imagebin_delete_url")
        req = urllib.request.Request(
            delete_url.format(delete_hash=args.delete_hash),
            headers=headers, method='DELETE')
        try:
            response = urllib.request.urlopen(req).read().decode()
        except urllib.error.HTTPError as e:
            print("ERROR", e)
            return {'deleted': False}

        rv = json.loads(response)
        if 'success' in rv and rv["success"]:
            print("\nImage successfully deleted!\n")
        else:
            print()
            print("Something has gone wrong!  Here is the server response:")
            print()
            pprint.pprint(rv)


class YTUploadImageCmd(YTCommand):
    args = (dict(short="file", type=str),)
    description = \
        """
        Upload an image to imgur.com.  Must be PNG.

        """
    name = "upload_image"
    def __call__(self, args):
        filename = args.file
        if not filename.endswith(".png"):
            print("File must be a PNG file!")
            return 1
        headers = {'Authorization':
            'Client-ID {}'.format(ytcfg.get("yt", "imagebin_api_key"))}

        image_data = base64.b64encode(open(filename, 'rb').read())
        parameters = {'image': image_data, type: 'base64',
                      'name': filename,
                      'title': "%s uploaded by yt" % filename}
        data = urllib.parse.urlencode(parameters).encode('utf-8')
        req = urllib.request.Request(
            ytcfg.get("yt", "imagebin_upload_url"), data=data, headers=headers)
        try:
            response = urllib.request.urlopen(req).read().decode()
        except urllib.error.HTTPError as e:
            print("ERROR", e)
            return {'uploaded':False}
        rv = json.loads(response)
        if 'data' in rv and 'link' in rv['data']:
            print()
            print("Image successfully uploaded!  You can find it at:")
            print("    %s" % (rv['data']['link']))
            print()
            print("If you'd like to delete it, use the following")
            print("    yt delete_image %s" % rv['data']['deletehash'])
            print()
        else:
            print()
            print("Something has gone wrong!  Here is the server response:")
            print()
            pprint.pprint(rv)

class YTSearchCmd(YTCommand):
    args = (dict(short="-o", longname="--output",
                 action="store", type=str,
                 dest="output", default="yt_index.json",
                 help="File in which to place output"),
            dict(longname="--check-all", short="-a",
                 help="Attempt to load every file",
                 action="store_true", default=False,
                 dest="check_all"),
            dict(longname="--full", short="-f",
                 help="Output full contents of parameter file",
                 action="store_true", default=False,
                 dest="full_output"),
            )
    description = \
        """
        Attempt to find outputs that yt can recognize in directories.
        """
    name = "search"
    def __call__(self, args):
        from yt.utilities.parameter_file_storage import \
            output_type_registry
        candidates = []
        for base, dirs, files in os.walk(".", followlinks=True):
            print("(% 10i candidates) Examining %s" % (len(candidates), base))
            recurse = []
            if args.check_all:
                candidates.extend([os.path.join(base, _) for _ in files])
            for _, otr in sorted(output_type_registry.items()):
                c, r = otr._guess_candidates(base, dirs, files)
                candidates.extend([os.path.join(base, _) for _ in c])
                recurse.append(r)
            if len(recurse) > 0 and not all(recurse):
                del dirs[:]
        # Now we have a ton of candidates.  We're going to do something crazy
        # and try to load each one.
        records = []
        for i, c in enumerate(sorted(candidates)):
            print("(% 10i/% 10i) Evaluating %s" % (i, len(candidates), c))
            try:
                record = get_metadata(c, args.full_output)
            except YTOutputNotIdentified:
                continue
            records.append(record)
        with open(args.output, "w") as f:
            json.dump(records, f, indent=4)
        print("Identified %s records output to %s" % (
              len(records), args.output))

def run_main():
    args = parser.parse_args()
    # The following is a workaround for a nasty Python 3 bug:
    # http://bugs.python.org/issue16308
    # http://bugs.python.org/issue9253
    try:
        getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
        
    args.func(args)

if __name__ == "__main__": run_main()
