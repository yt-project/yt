"""
A means of running standalone commands with a shared set of options.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Matthew Turk.  All Rights Reserved.

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
import optparse, os, os.path, math, sys, time, subprocess

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
                   dest="unit", default='1',
                   help="Desired units"),
    center  = dict(short="-c", long="--center",
                   action="store", type="float",
                   dest="center", default=[0.5, 0.5, 0.5],
                   nargs=3,
                   help="Center (-1,-1,-1 for max)"),
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

class YTCommands(cmdln.Cmdln):
    name="yt"

    def __init__(self, *args, **kwargs):
        cmdln.Cmdln.__init__(self, *args, **kwargs)
        cmdln.Cmdln.do_help.aliases.append("h")

    def do_loop(self, subcmd, opts, *args):
        """
        Interactive loop

        ${cmd_option_list}
        """
        self.cmdloop()

    @cmdln.option("-u", "--update-source", action="store_true",
                  default = False,
                  help="Update the yt installation, if able")
    @cmdln.option("-o", "--output-version", action="store",
                  default = None, dest="outputfile",
                  help="File into which the current revision number will be stored")
    def do_instinfo(self, subcmd, opts):
        """
        ${cmd_name}: Get some information about the yt installation

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
            print "You have to update this installation yourself."
            print
        if vstring is not None and opts.outputfile is not None:
            open(opts.outputfile, "w").write(vstring)

    def do_load(self, subcmd, opts, arg):
        """
        Load a single dataset into an IPython instance.

        ${cmd_option_list}
        """
        try:
            pf = _fix_pf(arg)
        except IOError:
            print "Could not load file."
            sys.exit()
        import yt.mods
        from IPython.Shell import IPShellEmbed
        local_ns = yt.mods.__dict__.copy()
        local_ns['pf'] = pf
        shell = IPShellEmbed()
        shell(local_ns = local_ns,
              header =
            "\nHi there!  Welcome to yt.\n\nWe've loaded your parameter file as 'pf'.  Enjoy!"
             )

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
        Run HaloProfiler on one dataset.

        ${cmd_option_list}
        """
        import yt.analysis_modules.halo_profiler.api as HP
        kwargs = {'halos': opts.halos,
                  'hop_style': opts.halo_hop_style,
                  'radius': opts.halo_radius,
                  'radius_units': opts.halo_radius_units}

        hp = HP.HaloProfiler(arg,opts.halo_parameter_file,**kwargs)
        if opts.make_profiles:
            hp.makeProfiles()
        if opts.make_projections:
            hp.makeProjections()

    @add_cmd_options(["maxw", "minw", "proj", "axis", "field", "weight",
                      "zlim", "nframes", "output", "cmap", "uboxes", "dex",
                      "text"])
    def do_zoomin(self, subcmd, opts, arg):
        """
        Create a set of zoomin frames

        ${cmd_option_list}
        """
        pf = _fix_pf(arg)
        min_width = opts.min_width * pf.h.get_smallest_dx()
        if opts.axis == 4:
            axes = range(3)
        else:
            axes = [opts.axis]
        pc = PlotCollection(pf)
        for ax in axes: 
            if opts.projection: p = pc.add_projection(opts.field, ax,
                                    weight_field=opts.weight)
            else: p = pc.add_slice(opts.field, ax)
            if opts.unit_boxes: p.modify["units"](factor=8)
            if opts.text is not None:
                p.modify["text"](
                    (0.02, 0.05), opts.text.replace(r"\n", "\n"),
                    text_args=dict(size="medium", color="w"))
        pc.set_width(opts.max_width,'1')
        # Check the output directory
        if not os.path.isdir(opts.output):
            os.mkdir(opts.output)
        # Figure out our zoom factor
        # Recall that factor^nframes = min_width / max_width
        # so factor = (log(min/max)/log(nframes))
        mylog.info("min_width: %0.3e max_width: %0.3e nframes: %0.3e",
                   min_width, opts.max_width, opts.nframes)
        factor=10**(math.log10(min_width/opts.max_width)/opts.nframes)
        mylog.info("Zoom factor: %0.3e", factor)
        w = 1.0
        for i in range(opts.nframes):
            mylog.info("Setting width to %0.3e", w)
            mylog.info("Saving frame %06i",i)
            pc.set_width(w,"1")
            if opts.zlim:
                pc.set_zlim(*opts.zlim)
            if opts.dex:
                pc.set_zlim('min', None, opts.dex)
            pc.set_cmap(opts.cmap)
            pc.save(os.path.join(opts.output,"%s_frame%06i" % (pf,i)))
            w = factor**i

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
        v, c = pf.h.find_max("Density")
        print "Maximum density: %0.5e at %s" % (v, c)
        if opts.output is not None:
            t = pf.current_time * pf['years']
            open(opts.output, "a").write(
                "%s (%0.5e years): %0.5e at %s\n" % (pf, t, v, c))

    @add_cmd_options([])
    def do_analyze(self, subcmd, opts, arg):
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
        hp.add_profile('CellVolume',weight_field=None,accumulation=True)
        hp.add_profile('TotalMassMsun',weight_field=None,accumulation=True)
        hp.add_profile('Density',weight_field=None,accumulation=False)
        hp.add_profile('Temperature',weight_field='CellMassMsun',accumulation=False)
        hp.make_profiles(filename="FilteredQuantities.out")

        # Add projection fields.
        hp.add_projection('Density',weight_field=None)
        hp.add_projection('Temperature',weight_field='Density')
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
        ph.modify["line"](pr.data["Density"], pr.data["Temperature"])
        pc.save()

    @cmdln.option("-d", "--desc", action="store",
                  default = None, dest="desc",
                  help="Description for this pasteboard entry")
    def do_pasteboard(self, subcmd, opts, arg):
        """
        Place a file into the user's pasteboard
        """
        if opts.desc is None: raise RuntimeError
        from yt.utilities.pasteboard import PostInventory
        pp = PostInventory()
        pp.add_post(arg, desc=opts.desc)

    @cmdln.option("-o", "--output", action="store",
                  default = None, dest="output_fn",
                  help="File to output to; else, print.")
    def do_pastegrab(self, subcmd, opts, username, paste_id):
        from yt.utilities.pasteboard import retrieve_pastefile
        retrieve_pastefile(username, paste_id, opts.output_fn)

    def do_bootstrap_dev(self, subcmd, opts):
        """
        Bootstrap a yt development environment
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
        supp_path = os.path.join(os.environ["YT_DEST"], "src",
                                 "yt-supplemental")
        # Now we check that the supplemental repository is checked out.
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
                print "$ hg clone http://hg.enzotools.org/yt-supplemental/ ",
                print "%s" % (supp_path)
                print
                sys.exit(1)
            rv = commands.clone(uu,
                    "http://hg.enzotools.org/yt-supplemental/", supp_path)
            if rv:
                print "Something has gone wrong.  Quitting."
                sys.exit(1)
        # Now we think we have our supplemental repository.
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
    def do_serve(self, subcmd, opts):
        """
        Run the Web GUI
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
        import yt.gui.reason.bottle as bottle
        from yt.gui.reason.extdirect_repl import ExtDirectREPL
        from yt.gui.reason.bottle_mods import uuid_serve_functions

        hr = ExtDirectREPL(base_extjs_path)
        bottle.debug()
        uuid_serve_functions(open_browser=opts.open_browser,
                    port=int(opts.port))


def run_main():
    for co in ["--parallel", "--paste"]:
        if co in sys.argv: del sys.argv[sys.argv.index(co)]
    YT = YTCommands()
    sys.exit(YT.main())

if __name__ == "__main__": run_main()
