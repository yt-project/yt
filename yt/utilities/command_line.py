"""
A means of running standalone commands with a shared set of options.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2008-2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from yt.mods import *
from yt.funcs import *
import cmdln as cmdln
import optparse, os, os.path, math, sys, time, subprocess, getpass, tempfile
import urllib, urllib2, base64

def _fix_pf(arg):
    if os.path.isdir("%s" % arg) and \
        os.path.exists("%s/%s" % (arg,arg)):
        pf = load("%s/%s" % (arg,arg))
    elif os.path.isdir("%s.dir" % arg) and \
        os.path.exists("%s.dir/%s" % (arg,arg)):
        pf = load("%s.dir/%s" % (arg,arg))
    elif arg.endswith(".hierarchy"):
        pf = load(arg[:-10])
    else:
        pf = load(arg)
    if pf is None:
        raise IOError
    return pf

_common_options = dict(
    axis    = dict(short="-a", long="--axis",
                   action="store", type="int",
                   dest="axis", default=4,
                   help="Axis (4 for all three)"),
    log     = dict(short="-l", long="--log",
                   action="store_true",
                   dest="takelog", default=True,
                   help="Take the log of the field?"),
    text    = dict(short="-t", long="--text",
                   action="store", type="string",
                   dest="text", default=None,
                   help="Textual annotation"),
    field   = dict(short="-f", long="--field",
                   action="store", type="string",
                   dest="field", default="Density",
                   help="Field to color by"),
    weight  = dict(short="-g", long="--weight",
                   action="store", type="string",
                   dest="weight", default=None,
                   help="Field to weight projections with"),
    cmap    = dict(short="", long="--colormap",
                   action="store", type="string",
                   dest="cmap", default="jet",
                   help="Colormap name"),
    zlim    = dict(short="-z", long="--zlim",
                   action="store", type="float",
                   dest="zlim", default=None,
                   nargs=2,
                   help="Color limits (min, max)"),
    dex     = dict(short="", long="--dex",
                   action="store", type="float",
                   dest="dex", default=None,
                   nargs=1,
                   help="Number of dex above min to display"),
    width   = dict(short="-w", long="--width",
                   action="store", type="float",
                   dest="width", default=1.0,
                   help="Width in specified units"),
    unit    = dict(short="-u", long="--unit",
                   action="store", type="string",
                   dest="unit", default='unitary',
                   help="Desired units"),
    center  = dict(short="-c", long="--center",
                   action="store", type="float",
                   dest="center", default=None,
                   nargs=3,
                   help="Center, space separated (-1 -1 -1 for max)"),
    bn      = dict(short="-b", long="--basename",
                   action="store", type="string",
                   dest="basename", default=None,
                   help="Basename of parameter files"),
    output  = dict(short="-o", long="--output",
                   action="store", type="string",
                   dest="output", default="frames/",
                   help="Folder in which to place output images"),
    outputfn= dict(short="-o", long="--output",
                   action="store", type="string",
                   dest="output", default=None,
                   help="File in which to place output"),
    skip    = dict(short="-s", long="--skip",
                   action="store", type="int",
                   dest="skip", default=1,
                   help="Skip factor for outputs"),
    proj    = dict(short="-p", long="--projection",
                   action="store_true", 
                   dest="projection", default=False,
                   help="Use a projection rather than a slice"),
    maxw    = dict(short="", long="--max-width",
                   action="store", type="float",
                   dest="max_width", default=1.0,
                   help="Maximum width in code units"),
    minw    = dict(short="", long="--min-width",
                   action="store", type="float",
                   dest="min_width", default=50,
                   help="Minimum width in units of smallest dx (default: 50)"),
    nframes = dict(short="-n", long="--nframes",
                   action="store", type="int",
                   dest="nframes", default=100,
                   help="Number of frames to generate"),
    slabw   = dict(short="", long="--slab-width",
                   action="store", type="float",
                   dest="slab_width", default=1.0,
                   help="Slab width in specified units"),
    slabu   = dict(short="-g", long="--slab-unit",
                   action="store", type="string",
                   dest="slab_unit", default='1',
                   help="Desired units for the slab"),
    ptype   = dict(short="", long="--particle-type",
                   action="store", type="int",
                   dest="ptype", default=2,
                   help="Particle type to select"),
    agecut  = dict(short="", long="--age-cut",
                   action="store", type="float",
                   dest="age_filter", default=None,
                   nargs=2,
                   help="Bounds for the field to select"),
    uboxes  = dict(short="", long="--unit-boxes",
                   action="store_true",
                   dest="unit_boxes",
                   help="Display helpful unit boxes"),
    thresh  = dict(short="", long="--threshold",
                   action="store", type="float",
                   dest="threshold", default=None,
                   help="Density threshold"),
    dm_only = dict(short="", long="--all-particles",
                   action="store_false", 
                   dest="dm_only", default=True,
                   help="Use all particles"),
    grids   = dict(short="", long="--show-grids",
                   action="store_true",
                   dest="grids", default=False,
                   help="Show the grid boundaries"),
    time    = dict(short="", long="--time",
                   action="store_true",
                   dest="time", default=False,
                   help="Print time in years on image"),
    contours    = dict(short="", long="--contours",
                   action="store",type="int",
                   dest="contours", default=None,
                   help="Number of Contours for Rendering"),
    contour_width  = dict(short="", long="--contour_width",
                   action="store",type="float",
                   dest="contour_width", default=None,
                   help="Width of gaussians used for rendering."),
    enhance   = dict(short="", long="--enhance",
                   action="store_true",
                   dest="enhance", default=False,
                   help="Enhance!"),
    valrange  = dict(short="-r", long="--range",
                   action="store", type="float",
                   dest="valrange", default=None,
                   nargs=2,
                   help="Range, space separated"),
    up  = dict(short="", long="--up",
                   action="store", type="float",
                   dest="up", default=None,
                   nargs=3,
                   help="Up, space separated"),
    viewpoint  = dict(short="", long="--viewpoint",
                   action="store", type="float",
                   dest="viewpoint", default=[1., 1., 1.],
                   nargs=3,
                   help="Viewpoint, space separated"),
    pixels    = dict(short="", long="--pixels",
                   action="store",type="int",
                   dest="pixels", default=None,
                   help="Number of Pixels for Rendering"),
    halos   = dict(short="", long="--halos",
                   action="store", type="string",
                   dest="halos",default="multiple",
                   help="Run halo profiler on a 'single' halo or 'multiple' halos."),
    halo_radius = dict(short="", long="--halo_radius",
                       action="store", type="float",
                       dest="halo_radius",default=0.1,
                       help="Constant radius for profiling halos if using hop output files with no radius entry. Default: 0.1."),
    halo_radius_units = dict(short="", long="--halo_radius_units",
                             action="store", type="string",
                             dest="halo_radius_units",default="1",
                             help="Units for radius used with --halo_radius flag. Default: '1' (code units)."),
    halo_hop_style = dict(short="", long="--halo_hop_style",
                          action="store", type="string",
                          dest="halo_hop_style",default="new",
                          help="Style of hop output file.  'new' for yt_hop files and 'old' for enzo_hop files."),
    halo_parameter_file = dict(short="", long="--halo_parameter_file",
                               action="store", type="string",
                               dest="halo_parameter_file",default=None,
                               help="HaloProfiler parameter file."),
    make_profiles = dict(short="", long="--make_profiles",
                         action="store_true", default=False,
                         help="Make profiles with halo profiler."),
    make_projections = dict(short="", long="--make_projections",
                            action="store_true", default=False,
                            help="Make projections with halo profiler.")

    )

def _add_options(parser, *options):
    for opt in options:
        oo = _common_options[opt].copy()
        parser.add_option(oo.pop("short"), oo.pop("long"), **oo)

def _get_parser(*options):
    parser = optparse.OptionParser()
    _add_options(parser, *options)
    return parser

def add_cmd_options(options):
    opts = []
    for option in options:
        vals = _common_options[option].copy()
        opts.append(([vals.pop("short"), vals.pop("long")],
                      vals))
    def apply_options(func):
        for args, kwargs in opts:
            func = cmdln.option(*args, **kwargs)(func)
        return func
    return apply_options

def check_args(func):
    @wraps(func)
    def arg_iterate(self, subcmd, opts, *args):
        if len(args) == 1:
            pfs = args
        elif len(args) == 2 and opts.basename is not None:
            pfs = ["%s%04i" % (opts.basename, r)
                   for r in range(int(args[0]), int(args[1]), opts.skip) ]
        else: pfs = args
        for arg in pfs:
            func(self, subcmd, opts, arg)
    return arg_iterate

def _update_hg(path, skip_rebuild = False):
    from mercurial import hg, ui, commands
    f = open(os.path.join(path, "yt_updater.log"), "a")
    u = ui.ui()
    u.pushbuffer()
    config_fn = os.path.join(path, ".hg", "hgrc")
    print "Reading configuration from ", config_fn
    u.readconfig(config_fn)
    repo = hg.repository(u, path)
    commands.pull(u, repo)
    f.write(u.popbuffer())
    f.write("\n\n")
    u.pushbuffer()
    commands.identify(u, repo)
    if "+" in u.popbuffer():
        print "Can't rebuild modules by myself."
        print "You will have to do this yourself.  Here's a sample commands:"
        print
        print "    $ cd %s" % (path)
        print "    $ hg up"
        print "    $ %s setup.py develop" % (sys.executable)
        return 1
    print "Updating the repository"
    f.write("Updating the repository\n\n")
    commands.update(u, repo, check=True)
    if skip_rebuild: return
    f.write("Rebuilding modules\n\n")
    p = subprocess.Popen([sys.executable, "setup.py", "build_ext", "-i"], cwd=path,
                        stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    stdout, stderr = p.communicate()
    f.write(stdout)
    f.write("\n\n")
    if p.returncode:
        print "BROKEN: See %s" % (os.path.join(path, "yt_updater.log"))
        sys.exit(1)
    f.write("Successful!\n")
    print "Updated successfully."

def _get_hg_version(path):
    from mercurial import hg, ui, commands
    u = ui.ui()
    u.pushbuffer()
    repo = hg.repository(u, path)
    commands.identify(u, repo)
    return u.popbuffer()

def get_yt_version():
    import pkg_resources
    yt_provider = pkg_resources.get_provider("yt")
    path = os.path.dirname(yt_provider.module_path)
    version = _get_hg_version(path)[:12]
    return version

# This code snippet is modified from Georg Brandl
def bb_apicall(endpoint, data, use_pass = True):
    uri = 'https://api.bitbucket.org/1.0/%s/' % endpoint
    # since bitbucket doesn't return the required WWW-Authenticate header when
    # making a request without Authorization, we cannot use the standard urllib2
    # auth handlers; we have to add the requisite header from the start
    if data is not None:
        data = urllib.urlencode(data)
    req = urllib2.Request(uri, data)
    if use_pass:
        username = raw_input("Bitbucket Username? ")
        password = getpass.getpass()
        upw = '%s:%s' % (username, password)
        req.add_header('Authorization', 'Basic %s' % base64.b64encode(upw).strip())
    return urllib2.urlopen(req).read()

def _get_yt_supp(uu):
    supp_path = os.path.join(os.environ["YT_DEST"], "src",
                             "yt-supplemental")
    # Now we check that the supplemental repository is checked out.
    from mercurial import hg, ui, commands
    if not os.path.isdir(supp_path):
        print
        print "*** The yt-supplemental repository is not checked ***"
        print "*** out.  I can do this for you, but because this ***"
        print "*** is a delicate act, I require you to respond   ***"
        print "*** to the prompt with the word 'yes'.            ***"
        print
        response = raw_input("Do you want me to try to check it out? ")
        if response != "yes":
            print
            print "Okay, I understand.  You can check it out yourself."
            print "This command will do it:"
            print
            print "$ hg clone http://hg.yt-project.org/yt-supplemental/ ",
            print "%s" % (supp_path)
            print
            sys.exit(1)
        rv = commands.clone(uu,
                "http://hg.yt-project.org/yt-supplemental/", supp_path)
        if rv:
            print "Something has gone wrong.  Quitting."
            sys.exit(1)
    # Now we think we have our supplemental repository.
    return supp_path

class YTCommands(cmdln.Cmdln):
    name="yt"

    def __init__(self, *args, **kwargs):
        cmdln.Cmdln.__init__(self, *args, **kwargs)
        cmdln.Cmdln.do_help.aliases.append("h")

    def do_update(self, subcmd, opts):
        """
        Update the yt installation to the most recent version

        ${cmd_usage}
        ${cmd_option_list}
        """
        import pkg_resources
        yt_provider = pkg_resources.get_provider("yt")
        path = os.path.dirname(yt_provider.module_path)
        print
        print "yt module located at:"
        print "    %s" % (path)
        update_supp = False
        if "YT_DEST" in os.environ:
            spath = os.path.join(
                     os.environ["YT_DEST"], "src", "yt-supplemental")
            if os.path.isdir(spath):
                print "The supplemental repositories are located at:"
                print "    %s" % (spath)
                update_supp = True
        vstring = None
        if "site-packages" not in path:
            vstring = _get_hg_version(path)
            print
            print "The current version of the code is:"
            print
            print "---"
            print vstring.strip()
            print "---"
            print
            print "This installation CAN be automatically updated."
            _update_hg(path)
            print "Updated successfully."
        else:
            print
            print "YT site-packages not in path, so you must"
            print "update this installation manually by committing and"
            print "merging your modifications to the code before"
            print "updating to the newest changeset."
            print

    @cmdln.option("-u", "--update-source", action="store_true",
                  default = False,
                  help="Update the yt installation, if able")
    @cmdln.option("-o", "--output-version", action="store",
                  default = None, dest="outputfile",
                  help="File into which the current revision number will be stored")
    def do_instinfo(self, subcmd, opts):
        """
        Get some information about the yt installation

        ${cmd_usage}
        ${cmd_option_list}
        """
        import pkg_resources
        yt_provider = pkg_resources.get_provider("yt")
        path = os.path.dirname(yt_provider.module_path)
        print
        print "yt module located at:"
        print "    %s" % (path)
        update_supp = False
        if "YT_DEST" in os.environ:
            spath = os.path.join(
                     os.environ["YT_DEST"], "src", "yt-supplemental")
            if os.path.isdir(spath):
                print "The supplemental repositories are located at:"
                print "    %s" % (spath)
                update_supp = True
        vstring = None
        if "site-packages" not in path:
            vstring = _get_hg_version(path)
            print
            print "The current version of the code is:"
            print
            print "---"
            print vstring.strip()
            print "---"
            print
            print "This installation CAN be automatically updated."
            if opts.update_source:  
                _update_hg(path)
            print "Updated successfully."
        elif opts.update_source:
            print
            print "YT site-packages not in path, so you must"
            print "update this installation manually by committing and"
            print "merging your modifications to the code before"
            print "updating to the newest changeset."
            print
        if vstring is not None and opts.outputfile is not None:
            open(opts.outputfile, "w").write(vstring)

    def do_load(self, subcmd, opts, arg):
        """
        Load a single dataset into an IPython instance

        ${cmd_option_list}
        """
        try:
            pf = _fix_pf(arg)
        except IOError:
            print "Could not load file."
            sys.exit()
        import yt.mods

        import IPython
        if IPython.__version__.startswith("0.10"):
            api_version = '0.10'
        elif IPython.__version__.startswith("0.11") or \
             IPython.__version__.startswith("0.12"):
            api_version = '0.11'

        local_ns = yt.mods.__dict__.copy()
        local_ns['pf'] = pf

        if api_version == '0.10':
            shell = IPython.Shell.IPShellEmbed()
            shell(local_ns = local_ns,
                  header =
                  "\nHi there!  Welcome to yt.\n\nWe've loaded your parameter file as 'pf'.  Enjoy!"
                  )
        else:
            from IPython.config.loader import Config
            cfg = Config()
            cfg.InteractiveShellEmbed.local_ns = local_ns
            IPython.embed(config=cfg)
            from IPython.frontend.terminal.embed import InteractiveShellEmbed
            ipshell = InteractiveShellEmbed(config=cfg)

    @add_cmd_options(['outputfn','bn','thresh','dm_only','skip'])
    @check_args
    def do_hop(self, subcmd, opts, arg):
        """
        Run HOP on one or more datasets

        ${cmd_option_list}
        """
        pf = _fix_pf(arg)
        kwargs = {'dm_only' : opts.dm_only}
        if opts.threshold is not None: kwargs['threshold'] = opts.threshold
        hop_list = HaloFinder(pf, **kwargs)
        if opts.output is None: fn = "%s.hop" % pf
        else: fn = opts.output
        hop_list.write_out(fn)

    @add_cmd_options(['make_profiles','make_projections','halo_parameter_file',
                      'halos','halo_hop_style','halo_radius','halo_radius_units'])
    def do_halos(self, subcmd, opts, arg):
        """
        Run HaloProfiler on one dataset

        ${cmd_option_list}
        """
        import yt.analysis_modules.halo_profiler.api as HP
        kwargs = {'halos': opts.halos,
                  'halo_radius': opts.halo_radius,
                  'radius_units': opts.halo_radius_units}

        hp = HP.HaloProfiler(arg,opts.halo_parameter_file,**kwargs)
        if opts.make_profiles:
            hp.make_profiles()
        if opts.make_projections:
            hp.make_projections()

    @add_cmd_options(["width", "unit", "bn", "proj", "center",
                      "zlim", "axis", "field", "weight", "skip",
                      "cmap", "output", "grids", "time"])
    @check_args
    def do_plot(self, subcmd, opts, arg):
        """
        Create a set of images

        ${cmd_usage}
        ${cmd_option_list}
        """
        pf = _fix_pf(arg)
        center = opts.center
        if opts.center == (-1,-1,-1):
            mylog.info("No center fed in; seeking.")
            v, center = pf.h.find_max("Density")
        elif opts.center is None:
            center = 0.5*(pf.domain_left_edge + pf.domain_right_edge)
        center = na.array(center)
        pc=PlotCollection(pf, center=center)
        if opts.axis == 4:
            axes = range(3)
        else:
            axes = [opts.axis]
        for ax in axes:
            mylog.info("Adding plot for axis %i", ax)
            if opts.projection: pc.add_projection(opts.field, ax,
                                    weight_field=opts.weight, center=center)
            else: pc.add_slice(opts.field, ax, center=center)
            if opts.grids: pc.plots[-1].modify["grids"]()
            if opts.time: 
                time = pf.current_time*pf['Time']*pf['years']
                pc.plots[-1].modify["text"]((0.2,0.8), 't = %5.2e yr'%time)
        pc.set_width(opts.width, opts.unit)
        pc.set_cmap(opts.cmap)
        if opts.zlim: pc.set_zlim(*opts.zlim)
        if not os.path.isdir(opts.output): os.makedirs(opts.output)
        pc.save(os.path.join(opts.output,"%s" % (pf)))

    @add_cmd_options(["proj", "field", "weight"])
    @cmdln.option("-a", "--axis", action="store", type="int",
                   dest="axis", default=0, help="Axis (4 for all three)")
    @cmdln.option("-o", "--host", action="store", type="string",
                   dest="host", default=None, help="IP Address to bind on")
    @check_args
    def do_mapserver(self, subcmd, opts, arg):
        """
        Serve a plot in a GMaps-style interface

        ${cmd_usage}
        ${cmd_option_list}
        """
        pf = _fix_pf(arg)
        pc=PlotCollection(pf, center=0.5*(pf.domain_left_edge +
                                          pf.domain_right_edge))
        if opts.axis == 4:
            print "Doesn't work with multiple axes!"
            return
        if opts.projection:
            p = pc.add_projection(opts.field, opts.axis, weight_field=opts.weight)
        else:
            p = pc.add_slice(opts.field, opts.axis)
        from yt.gui.reason.pannable_map import PannableMapServer
        mapper = PannableMapServer(p.data, opts.field)
        import yt.utilities.bottle as bottle
        bottle.debug(True)
        if opts.host is not None:
            colonpl = opts.host.find(":")
            if colonpl >= 0:
                port = int(opts.host.split(":")[-1])
                opts.host = opts.host[:colonpl]
            else:
                port = 8080
            bottle.run(server='rocket', host=opts.host, port=port)
        else:
            bottle.run(server='rocket')

    def do_rpdb(self, subcmd, opts, task):
        """
        Connect to a currently running (on localhost) rpd session.

        Commands run with --rpdb will trigger an rpdb session with any
        uncaught exceptions.

        ${cmd_usage} 
        ${cmd_option_list}
        """
        import rpdb
        rpdb.run_rpdb(int(task))

    @add_cmd_options(['outputfn','bn','skip'])
    @check_args
    def do_stats(self, subcmd, opts, arg):
        """
        Print stats and maximum density for one or more datasets

        ${cmd_option_list}
        """
        pf = _fix_pf(arg)
        pf.h.print_stats()
        if "Density" in pf.h.field_list:
            v, c = pf.h.find_max("Density")
        print "Maximum density: %0.5e at %s" % (v, c)
        if opts.output is not None:
            t = pf.current_time * pf['years']
            open(opts.output, "a").write(
                "%s (%0.5e years): %0.5e at %s\n" % (pf, t, v, c))

    @add_cmd_options([])
    def _do_analyze(self, subcmd, opts, arg):
        """
        Produce a set of analysis for a given output.  This includes
        HaloProfiler results with r200, as per the recipe file in the cookbook,
        profiles of a number of fields, projections of average Density and
        Temperature, and distribution functions for Density and Temperature.

        ${cmd_option_list}
        """
        # We will do the following things:
        #   Halo profiling (default parameters ONLY)
        #   Projections: Density, Temperature
        #   Full-box distribution functions
        import yt.analysis_modules.halo_profiler.api as HP
        hp = HP.HaloProfiler(arg)
        # Add a filter to remove halos that have no profile points with overdensity
        # above 200, and with virial masses less than 1e14 solar masses.
        # Also, return the virial mass and radius to be written out to a file.
        hp.add_halo_filter(HP.VirialFilter,must_be_virialized=True,
                           overdensity_field='ActualOverdensity',
                           virial_overdensity=200, virial_filters=[],
                           virial_quantities=['TotalMassMsun','RadiusMpc'])

        # Add profile fields.
        pf = hp.pf
        all_fields = pf.h.field_list + pf.h.derived_field_list
        for field, wv, acc in HP.standard_fields:
            if field not in all_fields: continue
            hp.add_profile(field, wv, acc)
        hp.make_profiles(filename="FilteredQuantities.out")

        # Add projection fields.
        hp.add_projection('Density',weight_field=None)
        hp.add_projection('Temperature',weight_field='Density')
        if "Metallicity" in all_fields:
            hp.add_projection('Metallicity',weight_field='Density')

        # Make projections for all three axes using the filtered halo list and
        # save data to hdf5 files.
        hp.make_projections(save_cube=True,save_images=True,
                            halo_list='filtered',axes=[0,1,2])

        # Now we make full-box projections.
        pf = EnzoStaticOutput(arg)
        c = 0.5*(pf.domain_right_edge + pf.domain_left_edge)
        pc = PlotCollection(pf, center=c)
        for ax in range(3):
            pc.add_projection("Density", ax, "Density")
            pc.add_projection("Temperature", ax, "Density")
            pc.plots[-1].set_cmap("hot")

        # Time to add some phase plots
        dd = pf.h.all_data()
        ph = pc.add_phase_object(dd, ["Density", "Temperature", "CellMassMsun"],
                            weight=None)
        pc_dummy = PlotCollection(pf, center=c)
        pr = pc_dummy.add_profile_object(dd, ["Density", "Temperature"],
                            weight="CellMassMsun")
        ph.modify["line"](pr.field_data["Density"], pr.field_data["Temperature"])
        pc.save()

    @cmdln.option("-d", "--desc", action="store",
                  default = None, dest="desc",
                  help="Description for this pasteboard entry")
    def do_pasteboard(self, subcmd, opts, arg):
        """
        Place a file into your pasteboard.
        """
        if opts.desc is None: raise RuntimeError
        from yt.utilities.pasteboard import PostInventory
        pp = PostInventory()
        pp.add_post(arg, desc=opts.desc)

    @cmdln.option("-l", "--language", action="store",
                  default = None, dest="language",
                  help="Use syntax highlighter for the file in language")
    @cmdln.option("-L", "--languages", action="store_true",
                  default = False, dest="languages",
                  help="Retrive a list of supported languages")
    @cmdln.option("-e", "--encoding", action="store",
                  default = 'utf-8', dest="encoding",
                  help="Specify the encoding of a file (default is "
                        "utf-8 or guessing if available)")
    @cmdln.option("-b", "--open-browser", action="store_true",
                  default = False, dest="open_browser",
                  help="Open the paste in a web browser")
    @cmdln.option("-p", "--private", action="store_true",
                  default = False, dest="private",
                  help="Paste as private")
    @cmdln.option("-c", "--clipboard", action="store_true",
                  default = False, dest="clipboard",
                  help="File to output to; else, print.")
    def do_pastebin(self, subcmd, opts, arg):
        """
        Post a script to an anonymous pastebin

        Usage: yt pastebin [options] <script>

        ${cmd_option_list}
        """
        import yt.utilities.lodgeit as lo
        lo.main( arg, languages=opts.languages, language=opts.language,
                 encoding=opts.encoding, open_browser=opts.open_browser,
                 private=opts.private, clipboard=opts.clipboard)

    def do_pastebin_grab(self, subcmd, opts, arg):
        """
        Print an online pastebin to STDOUT for local use. Paste ID is 
        the number at the end of the url.  So to locally access pastebin:
        http://paste.yt-project.org/show/1688/

        Usage: yt pastebin_grab <Paste ID> 
        Ex: yt pastebin_grab 1688 > script.py

        """
        import yt.utilities.lodgeit as lo
        lo.main( None, download=arg )

    @cmdln.option("-o", "--output", action="store",
                  default = None, dest="output_fn",
                  help="File to output to; else, print.")
    def do_pasteboard_grab(self, subcmd, opts, username, paste_id):
        """
        Download from your or another user's pasteboard

        ${cmd_usage} 
        ${cmd_option_list}
        """
        from yt.utilities.pasteboard import retrieve_pastefile
        retrieve_pastefile(username, paste_id, opts.output_fn)

    def do_bugreport(self, subcmd, opts):
        """
        Report a bug in yt

        ${cmd_usage} 
        ${cmd_option_list}
        """
        print "==============================================================="
        print
        print "Hi there!  Welcome to the yt bugreport taker."
        print
        print "==============================================================="
        print "At any time in advance of the upload of the bug, you should feel free"
        print "to ctrl-C out and submit the bug report manually by going here:"
        print "   http://hg.yt-project.org/yt/issues/new"
        print 
        print "Also, in order to submit a bug through this interface, you"
        print "need a Bitbucket account. If you don't have one, exit this "
        print "bugreport now and run the 'yt bootstrap_dev' command to create one."
        print
        print "Have you checked the existing bug reports to make"
        print "sure your bug has not already been recorded by someone else?"
        print "   http://hg.yt-project.org/yt/issues?status=new&status=open"
        print
        print "Finally, are you sure that your bug is, in fact, a bug? It might"
        print "simply be a misunderstanding that could be cleared up by"
        print "visiting the yt irc channel or getting advice on the email list:"
        print "   http://yt-project.org/irc.html"
        print "   http://lists.spacepope.org/listinfo.cgi/yt-users-spacepope.org"
        print
        summary = raw_input("Press <enter> if you remain firm in your conviction to continue.")
        print
        print
        print "Okay, sorry about that. How about a nice, pithy ( < 12 words )"
        print "summary of the bug?  (e.g. 'Particle overlay problem with parallel "
        print "projections')"
        print
        try:
            current_version = get_yt_version()
        except:
            current_version = "Unavailable"
        summary = raw_input("Summary? ")
        bugtype = "bug"
        data = dict(title = summary, type=bugtype)
        print
        print "Okay, now let's get a bit more information."
        print
        print "Remember that if you want to submit a traceback, you can run"
        print "any script with --paste or --detailed-paste to submit it to"
        print "the pastebin and then include the link in this bugreport."
        if "EDITOR" in os.environ:
            print
            print "Press enter to spawn your editor, %s" % os.environ["EDITOR"]
            loki = raw_input()
            tf = tempfile.NamedTemporaryFile(delete=False)
            fn = tf.name
            tf.close()
            popen = subprocess.call("$EDITOR %s" % fn, shell = True)
            content = open(fn).read()
            try:
                os.unlink(fn)
            except:
                pass
        else:
            print
            print "Couldn't find an $EDITOR variable.  So, let's just take"
            print "take input here.  Type up your summary until you're ready"
            print "to be done, and to signal you're done, type --- by itself"
            print "on a line to signal your completion."
            print
            print "(okay, type now)"
            print
            lines = []
            while 1:
                line = raw_input()
                if line.strip() == "---": break
                lines.append(line)
            content = "\n".join(lines)
        content = "Reporting Version: %s\n\n%s" % (current_version, content)
        endpoint = "repositories/yt_analysis/yt/issues"
        data['content'] = content
        print
        print "==============================================================="
        print 
        print "Okay, we're going to submit with this:"
        print
        print "Summary: %s" % (data['title'])
        print
        print "---"
        print content
        print "---"
        print
        print "==============================================================="
        print
        print "Is that okay?  If not, hit ctrl-c.  Otherwise, enter means"
        print "'submit'.  Next we'll ask for your Bitbucket Username."
        print "If you don't have one, run the 'yt bootstrap_dev' command."
        print
        loki = raw_input()
        retval = bb_apicall(endpoint, data, use_pass=True)
        import json
        retval = json.loads(retval)
        url = "http://hg.yt-project.org/yt/issue/%s" % retval['local_id']
        print 
        print "==============================================================="
        print
        print "Thanks for your bug report!  Together we'll make yt totally bug free!"
        print "You can view bug report here:"
        print "   %s" % url
        print
        print "Keep in touch!"
        print

    def do_bootstrap_dev(self, subcmd, opts):
        """
        Bootstrap a yt development environment

        ${cmd_usage} 
        ${cmd_option_list}
        """
        from mercurial import hg, ui, commands
        import imp
        import getpass
        import json
        uu = ui.ui()
        print
        print "Hi there!  Welcome to the yt development bootstrap tool."
        print
        print "This should get you started with mercurial as well as a few"
        print "other handy things, like a pasteboard of your very own."
        print
        # We have to do a couple things.
        # First, we check that YT_DEST is set.
        if "YT_DEST" not in os.environ:
            print
            print "*** You must set the environment variable YT_DEST ***"
            print "*** to point to the installation location!        ***"
            print
            sys.exit(1)
        supp_path = _get_yt_supp(uu)
        print
        print "I have found the yt-supplemental repository at %s" % (supp_path)
        print
        print "Let's load up and check what we need to do to get up and"
        print "running."
        print
        print "There are three stages:"
        print
        print " 1. Setting up your ~/.hgrc to have a username."
        print " 2. Setting up your bitbucket user account and the hgbb"
        print "    extension."
        print " 3. Setting up a new pasteboard repository."
        print
        firstname = lastname = email_address = bbusername = repo_list = None
        # Now we try to import the cedit extension.
        try:
            result = imp.find_module("cedit", [supp_path])
        except ImportError:
            print "I was unable to find the 'cedit' module in %s" % (supp_path)
            print "This may be due to a broken checkout."
            print "Sorry, but I'm going to bail."
            sys.exit(1)
        cedit = imp.load_module("cedit", *result)
        try:
            result = imp.find_module("hgbb", [supp_path + "/hgbb"])
        except ImportError:
            print "I was unable to find the 'hgbb' module in %s" % (supp_path)
            print "This may be due to a broken checkout."
            print "Sorry, but I'm going to bail."
            sys.exit(1)
        hgbb = imp.load_module("hgbb", *result)
        if uu.config("ui","username",None) is None:
            print "You don't have a username specified in your ~/.hgrc."
            print "Let's set this up.  If you would like to quit at any time,"
            print "hit Ctrl-C."
            print
            firstname = raw_input("What is your first name? ")
            lastname = raw_input("What is your last name? ")
            email_address = raw_input("What is your email address? ")
            print
            print "Thanks.  I will now add a username of this form to your"
            print "~/.hgrc file:"
            print
            print "    %s %s <%s>" % (firstname, lastname, email_address)
            print
            loki = raw_input("Press enter to go on, Ctrl-C to exit.")
            print
            cedit.setuser(uu, name="%s %s" % (firstname, lastname),
                              email="%s" % (email_address),
                              local=False, username=False)
            print
        else:
            print "Looks like you already have a username!"
            print "We'll skip that step, then."
            print
        print "Now we'll set up BitBucket user.  If you would like to do this"
        print "yourself, please visit:"
        print " https://bitbucket.org/account/signup/?plan=5_users"
        print "for a free account."
        print
        loki = raw_input("Do you have a BitBucket.org user already? [yes/no]")
        if loki.strip().upper() == "YES":
            bbusername = raw_input("Okay, cool.  What is your username?  ").strip()
            # Now we get information about the username.
            if firstname is None or lastname is None:
                rv = hgbb._bb_apicall(uu, "users/%s" % bbusername, None, False)
                rv = json.loads(rv)
                firstname = rv['user']["first_name"]
                lastname = rv['user']["last_name"]
                repo_list = rv['repositories']
                print "Retrieved your info:"
                print "  username:   %s" % (bbusername)
                print "  first name: %s" % (firstname)
                print "  last name:  %s" % (lastname)
        elif loki.strip().upper() == "NO":
            print "Okay, we can set you up with one.  It's probably better for"
            print "it to be all lowercase letter."
            print
            bbusername = raw_input("What is your desired username? ").strip()
            if firstname is None:
                firstname = raw_input("What's your first name? ").strip()
            if lastname is None:
                lastname = raw_input("What's your last name? ").strip()
            if email_address is None:
                email_address = raw_input("What's your email address? ").strip()
            print
            print "Okay, I'll see if I can create a user with this information:"
            print "  username:   %s" % (bbusername)
            print "  first name: %s" % (firstname)
            print "  last name:  %s" % (lastname)
            print "  email:      %s" % (email_address)
            print
            print "Now, I'm going to ask for a password.  This password will"
            print "be transmitted over HTTPS (not HTTP) and will not be stored"
            print "in any local file.  But, it will be stored in memory for"
            print "the duration of the user-creation process."
            print
            while 1:
                password1 = getpass.getpass("Password? ")
                password2 = getpass.getpass("Confirm? ")
                if password1 == password2: break
                print "Sorry, they didn't match!  Let's try again."
                print
            rv = hgbb._bb_apicall(uu, "newuser",
                                    dict(username=bbusername,
                                         password=password1,
                                         email=email_address,
                                         first_name = firstname,
                                         last_name = lastname),
                                   False)
            del password1, password2
            if str(json.loads(rv)['username']) == bbusername:
                print "Successful!  You probably just got an email asking you"
                print "to confirm this."
            else:
                print "Okay, something is wrong.  Quitting!"
                sys.exit(1)
        else:
            print "Not really sure what you replied with.  Quitting!"
            sys.exit(1)
        # We're now going to do some modification of the hgrc.
        # We need an hgrc first.
        hgrc_path = [cedit.config.defaultpath("user", uu)]
        hgrc_path = cedit.config.verifypaths(hgrc_path)
        # Now we set up the hgbb extension
        if uu.config("extensions","config",None) is None:
            # cedit is a module, but hgbb is a file.  So we call dirname here.
            cedit_path = os.path.dirname(cedit.__file__)
            print "Now we're going to turn on the cedit extension in:"
            print "    ", hgrc_path
            print "This will enable you to edit your configuration from the"
            print "command line.  Mainly, this is useful to use the command"
            print "'addsource', which will let you add a new source to a given"
            print "mercurial repository -- like a fork, or your own fork."
            print
            print "This constitutes adding the path to the cedit extension,"
            print "which will look like this:"
            print
            print "   [extensions]"
            print "   config=%s" % cedit_path
            print
            loki = raw_input("Press enter to go on, Ctrl-C to exit.")
            cedit.config.setoption(uu, hgrc_path, "extensions.config=%s" % cedit_path)
        if uu.config("extensions","hgbb",None) is None:
            hgbb_path = hgbb.__file__
            if hgbb_path.endswith(".pyc"): hgbb_path = hgbb_path[:-1]
            print "Now we're going to turn on the hgbb extension in:"
            print "    ", hgrc_path
            print "This will enable you to access BitBucket more easily, as well"
            print "as create repositories from the command line."
            print
            print "This constitutes adding the path to the hgbb extension,"
            print "which will look like this:"
            print
            print "   [extensions]"
            print "   hgbb=%s" % hgbb_path
            print
            loki = raw_input("Press enter to go on, Ctrl-C to exit.")
            cedit.config.setoption(uu, hgrc_path, "extensions.hgbb=%s" % hgbb_path)
        if uu.config("bb","username", None) is None:
            print "We'll now set up your username for BitBucket."
            print "We will add this:"
            print
            print "   [bb]"
            print "   username = %s" % (bbusername)
            print
            loki = raw_input("Press enter to go on, Ctrl-C to exit.")
            cedit.config.setoption(uu, hgrc_path, "bb.username=%s" % bbusername)
        bb_fp = "81:2b:08:90:dc:d3:71:ee:e0:7c:b4:75:ce:9b:6c:48:94:56:a1:fe"
        if uu.config("hostfingerprints", "bitbucket.org", None) is None:
            print "Let's also add bitbucket.org to the known hosts, so hg"
            print "doesn't warn us about bitbucket."
            print "We will add this:"
            print
            print "   [hostfingerprints]"
            print "   bitbucket.org = %s" % (bb_fp)
            print
            loki = raw_input("Press enter to go on, Ctrl-C to exit.")
            cedit.config.setoption(uu, hgrc_path,
                                   "hostfingerprints.bitbucket.org=%s" % bb_fp)

        # We now reload the UI's config file so that it catches the [bb]
        # section changes.
        uu.readconfig(hgrc_path[0])
        # Now the only thing remaining to do is to set up the pasteboard
        # repository.
        # This is, unfortunately, the most difficult.
        print
        print "We are now going to set up a pasteboard. This is a mechanism"
        print "for versioned posting of snippets, collaboration and"
        print "discussion."
        print
        # Let's get the full list of repositories
        pasteboard_name = "%s.bitbucket.org" % (bbusername.lower())
        if repo_list is None:
            rv = hgbb._bb_apicall(uu, "users/%s" % bbusername, None, False)
            rv = json.loads(rv)
            repo_list = rv['repositories']
        create = True
        for repo in repo_list:
            if repo['name'] == pasteboard_name:
                create = False
        if create:
            # Now we first create the repository, but we
            # will only use the creation API, not the bbcreate command.
            print
            print "I am now going to create the repository:"
            print "    ", pasteboard_name
            print "on BitBucket.org.  This will set up the domain"
            print "     http://%s" % (pasteboard_name)
            print "which will point to the current contents of the repo."
            print
            loki = raw_input("Press enter to go on, Ctrl-C to exit.")
            data = dict(name=pasteboard_name)
            hgbb._bb_apicall(uu, 'repositories', data)
        # Now we clone
        pasteboard_path = os.path.join(os.environ["YT_DEST"], "src",
                                       pasteboard_name)
        if os.path.isdir(pasteboard_path):
            print "Found an existing clone of the pasteboard repo:"
            print "    ", pasteboard_path
        else:
            print
            print "I will now clone a copy of your pasteboard repo."
            print
            loki = raw_input("Press enter to go on, Ctrl-C to exit.")
            commands.clone(uu, "https://%s@bitbucket.org/%s/%s" % (
                             bbusername, bbusername, pasteboard_name),
                           pasteboard_path)
            pbtemplate_path = os.path.join(supp_path, "pasteboard_template")
            pb_hgrc_path = os.path.join(pasteboard_path, ".hg", "hgrc")
            cedit.config.setoption(uu, [pb_hgrc_path],
                                   "paths.pasteboard = " + pbtemplate_path)
            if create:
                # We have to pull in the changesets from the pasteboard.
                pb_repo = hg.repository(uu, pasteboard_path)
                commands.pull(uu, pb_repo,
                              os.path.join(supp_path, "pasteboard_template"))
        if ytcfg.get("yt","pasteboard_repo") != pasteboard_path:
            print
            print "Now setting the pasteboard_repo option in"
            print "~/.yt/config to point to %s" % (pasteboard_path)
            print
            loki = raw_input("Press enter to go on, Ctrl-C to exit.")
            dotyt_path = os.path.expanduser("~/.yt")
            if not os.path.isdir(dotyt_path):
                print "There's no directory:"
                print "    ", dotyt_path
                print "I will now create it."
                print
                loki = raw_input("Press enter to go on, Ctrl-C to exit.")
                os.mkdir(dotyt_path)
            ytcfg_path = os.path.expanduser("~/.yt/config")
            cedit.config.setoption(uu, [ytcfg_path],
                        "yt.pasteboard_repo=%s" % (pasteboard_path))
        try:
            import pygments
            install_pygments = False
        except ImportError:
            install_pygments = True
        if install_pygments:
            print "You are missing the Pygments package.  Installing."
            import pip
            rv = pip.main(["install", "pygments"])
            if rv == 1:
                print "Unable to install Pygments.  Please report this bug to yt-users."
                sys.exit(1)
        try:
            import lxml
            install_lxml = False
        except ImportError:
            install_lxml = True
        if install_lxml:
            print "You are missing the lxml package.  Installing."
            import pip
            rv = pip.main(["install", "lxml"])
            if rv == 1:
                print "Unable to install lxml.  Please report this bug to yt-users."
                sys.exit(1)
        print
        print "All done!"
        print
        print "You're now set up to use the 'yt pasteboard' command"
        print "as well as develop using Mercurial and BitBucket."
        print
        print "Good luck!"

    @cmdln.option("-o", "--open-browser", action="store_true",
                  default = False, dest='open_browser',
                  help="Open a web browser.")
    @cmdln.option("-p", "--port", action="store",
                  default = 0, dest='port',
                  help="Port to listen on")
    @cmdln.option("-f", "--find", action="store_true",
                  default = False, dest="find",
                  help="At startup, find all *.hierarchy files in the CWD")
    @cmdln.option("-d", "--debug", action="store_true",
                  default = False, dest="debug",
                  help="Add a debugging mode for cell execution")
    def do_serve(self, subcmd, opts):
        """
        Run the Web GUI Reason
        """
        # We have to do a couple things.
        # First, we check that YT_DEST is set.
        if "YT_DEST" not in os.environ:
            print
            print "*** You must set the environment variable YT_DEST ***"
            print "*** to point to the installation location!        ***"
            print
            sys.exit(1)
        if opts.port == 0:
            # This means, choose one at random.  We do this by binding to a
            # socket and allowing the OS to choose the port for that socket.
            import socket
            sock = socket.socket()
            sock.bind(('', 0))
            opts.port = sock.getsockname()[-1]
            del sock
        elif opts.port == '-1':
            port = raw_input("Desired yt port? ")
            try:
                opts.port = int(port)
            except ValueError:
                print "Please try a number next time."
                return 1
        base_extjs_path = os.path.join(os.environ["YT_DEST"], "src")
        if not os.path.isfile(os.path.join(base_extjs_path, "ext-resources", "ext-all.js")):
            print
            print "*** You are missing the ExtJS support files. You  ***"
            print "*** You can get these by either rerunning the     ***"
            print "*** install script installing, or downloading     ***"
            print "*** them manually.                                ***"
            print
            sys.exit(1)
        from yt.config import ytcfg;ytcfg["yt","__withinreason"]="True"
        import yt.utilities.bottle as bottle
        from yt.gui.reason.extdirect_repl import ExtDirectREPL
        from yt.gui.reason.bottle_mods import uuid_serve_functions, PayloadHandler
        hr = ExtDirectREPL(base_extjs_path)
        hr.debug = PayloadHandler.debug = opts.debug
        if opts.find:
            # We just have to find them and store references to them.
            command_line = ["pfs = []"]
            for fn in sorted(glob.glob("*/*.hierarchy")):
                command_line.append("pfs.append(load('%s'))" % fn[:-10])
            hr.execute("\n".join(command_line))
        bottle.debug()
        uuid_serve_functions(open_browser=opts.open_browser,
                    port=int(opts.port), repl=hr)

    
    def _do_remote(self, subcmd, opts):
        import getpass, sys, socket, time, webbrowser
        import yt.utilities.pexpect as pex

        host = raw_input('Hostname: ')
        user = raw_input('User: ')
        password = getpass.getpass('Password: ')

        sock = socket.socket()
        sock.bind(('', 0))
        port = sock.getsockname()[-1]
        del sock

        child = pex.spawn('ssh -L %s:localhost:%s -l %s %s'%(port, port, user, host))
        ssh_newkey = 'Are you sure you want to continue connecting'
        i = child.expect([pex.TIMEOUT, ssh_newkey, 'password: '])
        if i == 0: # Timeout
            print 'ERROR!'
            print 'SSH could not login. Here is what SSH said:'
            print child.before, child.after
            return 1
        if i == 1: # SSH does not have the public key. Just accept it.
            child.sendline ('yes')
            child.expect ('password: ')
            i = child.expect([pex.TIMEOUT, 'password: '])
            if i == 0: # Timeout
                print 'ERROR!'
                print 'SSH could not login. Here is what SSH said:'
                print child.before, child.after
                return 1
        print "Sending password"
        child.sendline(password)
        del password
        print "Okay, sending serving command"
        child.sendline('yt serve -p -1')
        print "Waiting ..."
        child.expect('Desired yt port?')
        child.sendline("%s" % port)
        child.expect('     http://localhost:([0-9]*)/(.+)/\r')
        print "Got:", child.match.group(1), child.match.group(2)
        port, urlprefix = child.match.group(1), child.match.group(2)
        print "Sleeping one second and opening browser"
        time.sleep(1)
        webbrowser.open("http://localhost:%s/%s/" % (port, urlprefix))
        print "Press Ctrl-C to terminate session"
        child.readlines()
        while 1:
            time.sleep(1)

    @cmdln.option("-R", "--repo", action="store", type="string",
                  dest="repo", default=".", help="Repository to upload")
    def do_hubsubmit(self, subcmd, opts):
        """
        Submit a mercurial repository to the yt Hub
        (http://hub.yt-project.org/), creating a BitBucket repo in the process
        if necessary.

        ${cmd_usage}
        ${cmd_option_list}
        """
        import imp
        from mercurial import hg, ui, commands, error, config
        uri = "http://hub.yt-project.org/3rdparty/API/api.php"
        uu = ui.ui()
        supp_path = _get_yt_supp(uu)
        try:
            result = imp.find_module("cedit", [supp_path])
        except ImportError:
            print "I was unable to find the 'cedit' module in %s" % (supp_path)
            print "This may be due to a broken checkout."
            print "Sorry, but I'm going to bail."
            sys.exit(1)
        cedit = imp.load_module("cedit", *result)
        try:
            result = imp.find_module("hgbb", [supp_path + "/hgbb"])
        except ImportError:
            print "I was unable to find the 'hgbb' module in %s" % (supp_path)
            print "This may be due to a broken checkout."
            print "Sorry, but I'm going to bail."
            sys.exit(1)
        hgbb = imp.load_module("hgbb", *result)
        try:
            repo = hg.repository(uu, opts.repo)
            conf = config.config()
            if os.path.exists(os.path.join(opts.repo,".hg","hgrc")):
                conf.read(os.path.join(opts.repo, ".hg", "hgrc"))
            needs_bb = True
            if "paths" in conf.sections():
                default = conf['paths'].get("default", "")
                if default.startswith("bb://") or "bitbucket.org" in default:
                    needs_bb = False
                    bb_url = default
                else:
                    for alias, value in conf["paths"].items():
                        if value.startswith("bb://") or "bitbucket.org" in value:
                            needs_bb = False
                            bb_url = value
                            break
        except error.RepoError:
            print "Unable to find repo at:"
            print "   %s" % (os.path.abspath(opts.repo))
            print
            print "Would you like to initialize one?  If this message"
            print "surprises you, you should perhaps press Ctrl-C to quit."
            print "Otherwise, type 'yes' at the prompt."
            print
            loki = raw_input("Create repo? ")
            if loki.upper().strip() != "YES":
                print "Okay, rad -- we'll let you handle it and get back to",
                print " us."
                return 1
            commands.init(uu, dest=opts.repo)
            repo = hg.repository(uu, opts.repo)
            commands.add(uu, repo)
            commands.commit(uu, repo, message="Initial automated import by yt")
            needs_bb = True
        if needs_bb:
            print
            print "Your repository is not yet on BitBucket, as near as I can tell."
            print "Would you like to create a repository there and upload to it?"
            print "Without this, I don't know what URL to submit!"
            print
            print "Type 'yes' to accept."
            print
            loki = raw_input("Upload to BitBucket? ")
            if loki.upper().strip() != "YES": return 1
            hgrc_path = [cedit.config.defaultpath("user", uu)]
            hgrc_path = cedit.config.verifypaths(hgrc_path)
            uu.readconfig(hgrc_path[0])
            bb_username = uu.config("bb", "username", None)
            if bb_username is None:
                print "Can't find your Bitbucket username.  Run the command:"
                print
                print "$ yt bootstrap_dev"
                print
                print "to get set up and ready to go."
                return 1
            bb_repo_name = os.path.basename(os.path.abspath(opts.repo))
            print
            print "I am now going to create the repository:"
            print "    ", bb_repo_name
            print "on BitBucket.org and upload this repository to that."
            print "If that is not okay, please press Ctrl-C to quit."
            print
            loki = raw_input("Press Enter to continue.")
            data = dict(name=bb_repo_name)
            hgbb._bb_apicall(uu, 'repositories', data)
            print
            print "Created repository!  Now I will set this as the default path."
            print
            bb_url = "https://%s@bitbucket.org/%s/%s" % (
                        bb_username, bb_username, bb_repo_name)
            cedit.config.addsource(uu, repo, "default", bb_url)
            commands.push(uu, repo, bb_url)
        if bb_url.startswith("bb://"):
            bb_username, bb_repo_name = bb_url.split("/")[-2:]
            bb_url = "https://%s@bitbucket.org/%s/%s" % (
                bb_username, bb_username, bb_repo_name)
        # Now we can submit
        import xml.etree.ElementTree as etree
        print
        print "Okay.  Now we're ready to submit to the Hub."
        print "Remember, you can go to the Hub at any time at"
        print " http://hub.yt-project.org/"
        print
        print "(Especially if you don't have a user yet!  We can wait.)"
        print
        hub_username = raw_input("What is your Hub username? ")
        hub_password = getpass.getpass("What is your Hub password? ")
        data = urllib.urlencode(dict(fn = "list",
                                     username=hub_username,
                                     password=hub_password))
        req = urllib2.Request(uri, data)
        rv = urllib2.urlopen(req).read()
        try:
            cats = etree.fromstring(rv)
        except:
            print "I think you entered your password wrong.  Please check!"
            return

        categories = {}

        for cat in cats.findall("./cate"):
            cat_id = int(cat.findall("id")[0].text)
            cat_name = cat.findall("name")[0].text
            categories[cat_id] = cat_name

        print
        for i, n in sorted(categories.items()):
            print "%i. %s" % (i, n)
        print
        cat_id = int(raw_input("Which category number does your script fit into? "))
        print
        print "What is the title of your submission? (Usually a repository name) "
        title = raw_input("Title? ")
        print
        print "What tags should be applied to this submission?  Separate with commas."
        print "(e.g., enzo, flash, gadget, ramses, nyx, yt, visualization, analysis,"
        print " utility, cosmology)"
        print
        tags = raw_input("Tags? ")
        print
        print "Give us a very brief summary of the project -- enough to get someone"
        print "interested enough to click the link and see what it's about.  This"
        print "should be a few sentences at most."
        print
        summary = raw_input("Summary? ")
        print
        print "Okay, we're going to submit!  Press enter to submit, Ctrl-C to back out."
        print
        loki = raw_input()

        data = urllib.urlencode(dict(fn = "post",
                                     username=hub_username, password=hub_password,
                                     url = bb_url, category = cat_id, title = title,
                                     content = summary, tags = tags))
        req = urllib2.Request(uri, data)
        rv = urllib2.urlopen(req).read()
        print rv

    def do_upload_image(self, subcmd, opts, filename):
        """
        Upload an image to imgur.com.  Must be PNG.

        ${cmd_usage} 
        ${cmd_option_list}
        """
        if not filename.endswith(".png"):
            print "File must be a PNG file!"
            return 1
        import base64, json, pprint
        image_data = base64.b64encode(open(filename).read())
        api_key = 'f62d550859558f28c4c214136bc797c7'
        parameters = {'key':api_key, 'image':image_data, type:'base64',
                      'caption': "",
                      'title': "%s uploaded by yt" % filename}
        data = urllib.urlencode(parameters)
        req = urllib2.Request('http://api.imgur.com/2/upload.json', data)
        try:
            response = urllib2.urlopen(req).read()
        except urllib2.HTTPError as e:
            print "ERROR", e
            return {'uploaded':False}
        rv = json.loads(response)
        if 'upload' in rv and 'links' in rv['upload']:
            print
            print "Image successfully uploaded!  You can find it at:"
            print "    %s" % (rv['upload']['links']['imgur_page'])
            print
            print "If you'd like to delete it, visit this page:"
            print "    %s" % (rv['upload']['links']['delete_page'])
            print
        else:
            print
            print "Something has gone wrong!  Here is the server response:"
            print
            pprint.pprint(rv)

    @add_cmd_options(["width", "unit", "center","enhance",'outputfn',
                      "field", "cmap", "contours", "viewpoint",
                      "pixels","up","valrange","log","contour_width"])
    @check_args
    def do_render(self, subcmd, opts, arg):
        """
        Create a simple volume rendering

        ${cmd_usage}
        ${cmd_option_list}
        """
        pf = _fix_pf(arg)
        center = opts.center
        if opts.center == (-1,-1,-1):
            mylog.info("No center fed in; seeking.")
            v, center = pf.h.find_max("Density")
        elif opts.center is None:
            center = 0.5*(pf.domain_left_edge + pf.domain_right_edge)
        center = na.array(center)

        L = opts.viewpoint
        if L is None:
            L = [1.]*3
        L = na.array(opts.viewpoint)

        unit = opts.unit
        if unit is None:
            unit = '1'
        width = opts.width
        if width is None:
            width = 0.5*(pf.domain_right_edge - pf.domain_left_edge)
        width /= pf[unit]

        N = opts.pixels
        if N is None:
            N = 512 
        
        up = opts.up
        if up is None:
            up = [0.,0.,1.]
            
        field = opts.field
        if field is None:
            field = 'Density'
        
        log = opts.takelog
        if log is None:
            log = True

        myrange = opts.valrange
        if myrange is None:
            roi = pf.h.region(center, center-width, center+width)
            mi, ma = roi.quantities['Extrema'](field)[0]
            if log:
                mi, ma = na.log10(mi), na.log10(ma)
        else:
            mi, ma = myrange[0], myrange[1]

        n_contours = opts.contours
        if n_contours is None:
            n_contours = 7

        contour_width = opts.contour_width

        cmap = opts.cmap
        if cmap is None:
            cmap = 'jet'
        tf = ColorTransferFunction((mi-2, ma+2))
        tf.add_layers(n_contours,w=contour_width,col_bounds = (mi,ma), colormap=cmap)

        cam = pf.h.camera(center, L, width, (N,N), transfer_function=tf)
        image = cam.snapshot()

        if opts.enhance:
            for i in range(3):
                image[:,:,i] = image[:,:,i]/(image[:,:,i].mean() + 5.*image[:,:,i].std())
            image[image>1.0]=1.0
            
        save_name = opts.output
        if save_name is None:
            save_name = "%s"%pf+"_"+field+"_rendering.png"
        if not '.png' in save_name:
            save_name += '.png'
        if cam.comm.rank != -1:
            write_bitmap(image,save_name)
        

def run_main():
    for co in ["--parallel", "--paste"]:
        if co in sys.argv: del sys.argv[sys.argv.index(co)]
    YT = YTCommands()
    sys.exit(YT.main())

if __name__ == "__main__": run_main()
