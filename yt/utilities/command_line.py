import argparse
import base64
import getpass
import json
import os
import pprint
import sys
import textwrap
import urllib
import urllib.request
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from more_itertools import always_iterable
from tqdm import tqdm

from yt.config import ytcfg
from yt.funcs import (
    download_file,
    enable_plugins,
    ensure_dir,
    ensure_dir_exists,
    get_git_version,
    mylog,
    update_git,
)
from yt.loaders import load
from yt.utilities.exceptions import YTFieldNotParseable, YTUnidentifiedDataType
from yt.utilities.metadata import get_metadata
from yt.visualization.plot_window import ProjectionPlot, SlicePlot

# isort: off
# This needs to be set before importing startup_tasks
ytcfg["yt", "internals", "command_line"] = True  # isort: skip
from yt.startup_tasks import parser, subparsers  # isort: skip # noqa: E402

# isort: on

# loading field plugins for backward compatibility, since this module
# used to do "from yt.mods import *"

try:
    enable_plugins()
except FileNotFoundError:
    pass

_default_colormap = ytcfg.get("yt", "default_colormap")


def _fix_ds(arg, *args, **kwargs):
    if os.path.isdir(f"{arg}") and os.path.exists(f"{arg}/{arg}"):
        ds = load(f"{arg}/{arg}", *args, **kwargs)
    elif os.path.isdir(f"{arg}.dir") and os.path.exists(f"{arg}.dir/{arg}"):
        ds = load(f"{arg}.dir/{arg}", *args, **kwargs)
    elif arg.endswith(".index"):
        ds = load(arg[:-10], *args, **kwargs)
    else:
        ds = load(arg, *args, **kwargs)
    return ds


def _add_arg(sc, arg):
    if isinstance(arg, str):
        arg = _common_options[arg].copy()
    elif isinstance(arg, tuple):
        exclusive, *args = arg
        if exclusive:
            grp = sc.add_mutually_exclusive_group()
        else:
            grp = sc.add_argument_group()
        for arg in args:
            _add_arg(grp, arg)
        return
    argc = dict(arg.items())
    argnames = []
    if "short" in argc:
        argnames.append(argc.pop("short"))
    if "longname" in argc:
        argnames.append(argc.pop("longname"))
    sc.add_argument(*argnames, **argc)


def _print_failed_source_update(reinstall=False):
    print()
    print("The yt package is not installed from a git repository,")
    print("so you must update this installation manually.")
    if "Continuum Analytics" in sys.version or "Anaconda" in sys.version:
        # see http://stackoverflow.com/a/21318941/1382869 for why we need
        # to check both Continuum *and* Anaconda
        print()
        print("Since it looks like you are using a python installation")
        print("that is managed by conda, you may want to do:")
        print()
        print("    $ conda update yt")
        print()
        print("to update your yt installation.")
        if reinstall:
            print()
            print("To update all of your packages, you can do:")
            print()
            print("    $ conda update --all")
    else:
        print("If you manage your python dependencies with pip, you may")
        print("want to do:")
        print()
        print("    $ python -m pip install -U yt")
        print()
        print("to update your yt installation.")


def _print_installation_information(path):
    import yt

    print()
    print("yt module located at:")
    print(f"    {path}")
    if "YT_DEST" in os.environ:
        spath = os.path.join(os.environ["YT_DEST"], "src", "yt-supplemental")
        if os.path.isdir(spath):
            print("The supplemental repositories are located at:")
            print(f"    {spath}")
    print()
    print("The current version of yt is:")
    print()
    print("---")
    print(f"Version = {yt.__version__}")
    vstring = get_git_version(path)
    if vstring is not None:
        print(f"Changeset = {vstring.strip()}")
    print("---")
    return vstring


class FileStreamer:
    final_size = None
    next_sent = 0
    chunksize = 100 * 1024

    def __init__(self, f, final_size=None):
        location = f.tell()
        f.seek(0, os.SEEK_END)
        self.final_size = f.tell() - location
        f.seek(location)
        self.f = f

    def __iter__(self):
        with tqdm(
            total=self.final_size, desc="Uploading file", unit="B", unit_scale=True
        ) as pbar:
            while self.f.tell() < self.final_size:
                yield self.f.read(self.chunksize)
                pbar.update(self.chunksize)


_subparsers = {None: subparsers}
_subparsers_description = {
    "config": "Get and set configuration values for yt",
}


class YTCommandSubtype(type):
    def __init__(cls, name, b, d):
        type.__init__(cls, name, b, d)
        if cls.name is None:
            return
        if cls.subparser not in _subparsers:
            try:
                description = _subparsers_description[cls.subparser]
            except KeyError:
                description = cls.subparser
            parent_parser = argparse.ArgumentParser(add_help=False)
            p = subparsers.add_parser(
                cls.subparser,
                help=description,
                description=description,
                parents=[parent_parser],
            )
            _subparsers[cls.subparser] = p.add_subparsers(
                title=cls.subparser, dest=cls.subparser
            )
        sp = _subparsers[cls.subparser]
        for name in always_iterable(cls.name):
            sc = sp.add_parser(name, description=cls.description, help=cls.description)
            sc.set_defaults(func=cls.run)
            for arg in cls.args:
                _add_arg(sc, arg)


class YTCommand(metaclass=YTCommandSubtype):
    args: Tuple[Union[str, Dict[str, Any]], ...] = ()
    name: Optional[Union[str, List[str]]] = None
    description: str = ""
    aliases = ()
    ndatasets: int = 1
    subparser: Optional[str] = None

    @classmethod
    def run(cls, args):
        self = cls()
        # Check for some things we know; for instance, comma separated
        # field names should be parsed as tuples.
        if getattr(args, "field", None) is not None and "," in args.field:
            if args.field.count(",") > 1:
                raise YTFieldNotParseable(args.field)
            args.field = tuple(_.strip() for _ in args.field.split(","))
        if getattr(args, "weight", None) is not None and "," in args.weight:
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
                args.ds = getattr(args, "ds", [None])[0]
                self(args)


class GetParameterFiles(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if len(values) == 1:
            datasets = values
        elif len(values) == 2 and namespace.basename is not None:
            datasets = [
                "%s%04i" % (namespace.basename, r)
                for r in range(int(values[0]), int(values[1]), namespace.skip)
            ]
        else:
            datasets = values
        namespace.ds = [_fix_ds(ds) for ds in datasets]


_common_options = dict(
    all=dict(
        longname="--all",
        dest="reinstall",
        default=False,
        action="store_true",
        help=(
            "Reinstall the full yt stack in the current location. "
            "This option has been deprecated and will not have any "
            "effect."
        ),
    ),
    ds=dict(short="ds", action=GetParameterFiles, nargs="+", help="datasets to run on"),
    ods=dict(
        action=GetParameterFiles,
        dest="ds",
        nargs="*",
        help="(Optional) datasets to run on",
    ),
    axis=dict(
        short="-a",
        longname="--axis",
        action="store",
        type=int,
        dest="axis",
        default=4,
        help="Axis (4 for all three)",
    ),
    log=dict(
        short="-l",
        longname="--log",
        action="store_true",
        dest="takelog",
        default=True,
        help="Use logarithmic scale for image",
    ),
    linear=dict(
        longname="--linear",
        action="store_false",
        dest="takelog",
        help="Use linear scale for image",
    ),
    text=dict(
        short="-t",
        longname="--text",
        action="store",
        type=str,
        dest="text",
        default=None,
        help="Textual annotation",
    ),
    field=dict(
        short="-f",
        longname="--field",
        action="store",
        type=str,
        dest="field",
        default="density",
        help=("Field to color by, use a comma to separate field tuple values"),
    ),
    weight=dict(
        short="-g",
        longname="--weight",
        action="store",
        type=str,
        dest="weight",
        default=None,
        help=(
            "Field to weight projections with, "
            "use a comma to separate field tuple values"
        ),
    ),
    cmap=dict(
        longname="--colormap",
        action="store",
        type=str,
        dest="cmap",
        default=_default_colormap,
        help="Colormap name",
    ),
    zlim=dict(
        short="-z",
        longname="--zlim",
        action="store",
        type=float,
        dest="zlim",
        default=None,
        nargs=2,
        help="Color limits (min, max)",
    ),
    dex=dict(
        longname="--dex",
        action="store",
        type=float,
        dest="dex",
        default=None,
        nargs=1,
        help="Number of dex above min to display",
    ),
    width=dict(
        short="-w",
        longname="--width",
        action="store",
        type=float,
        dest="width",
        default=None,
        help="Width in specified units",
    ),
    unit=dict(
        short="-u",
        longname="--unit",
        action="store",
        type=str,
        dest="unit",
        default="1",
        help="Desired axes units",
    ),
    center=dict(
        short="-c",
        longname="--center",
        action="store",
        type=float,
        dest="center",
        default=None,
        nargs=3,
        help="Center, space separated (-1 -1 -1 for max)",
    ),
    max=dict(
        short="-m",
        longname="--max",
        action="store_true",
        dest="max",
        default=False,
        help="Center the plot on the density maximum",
    ),
    bn=dict(
        short="-b",
        longname="--basename",
        action="store",
        type=str,
        dest="basename",
        default=None,
        help="Basename of datasets",
    ),
    output=dict(
        short="-o",
        longname="--output",
        action="store",
        type=str,
        dest="output",
        default="frames/",
        help="Folder in which to place output images",
    ),
    outputfn=dict(
        short="-o",
        longname="--output",
        action="store",
        type=str,
        dest="output",
        default=None,
        help="File in which to place output",
    ),
    skip=dict(
        short="-s",
        longname="--skip",
        action="store",
        type=int,
        dest="skip",
        default=1,
        help="Skip factor for outputs",
    ),
    proj=dict(
        short="-p",
        longname="--projection",
        action="store_true",
        dest="projection",
        default=False,
        help="Use a projection rather than a slice",
    ),
    maxw=dict(
        longname="--max-width",
        action="store",
        type=float,
        dest="max_width",
        default=1.0,
        help="Maximum width in code units",
    ),
    minw=dict(
        longname="--min-width",
        action="store",
        type=float,
        dest="min_width",
        default=50,
        help="Minimum width in units of smallest dx (default: 50)",
    ),
    nframes=dict(
        short="-n",
        longname="--nframes",
        action="store",
        type=int,
        dest="nframes",
        default=100,
        help="Number of frames to generate",
    ),
    slabw=dict(
        longname="--slab-width",
        action="store",
        type=float,
        dest="slab_width",
        default=1.0,
        help="Slab width in specified units",
    ),
    slabu=dict(
        short="-g",
        longname="--slab-unit",
        action="store",
        type=str,
        dest="slab_unit",
        default="1",
        help="Desired units for the slab",
    ),
    ptype=dict(
        longname="--particle-type",
        action="store",
        type=int,
        dest="ptype",
        default=2,
        help="Particle type to select",
    ),
    agecut=dict(
        longname="--age-cut",
        action="store",
        type=float,
        dest="age_filter",
        default=None,
        nargs=2,
        help="Bounds for the field to select",
    ),
    uboxes=dict(
        longname="--unit-boxes",
        action="store_true",
        dest="unit_boxes",
        help="Display heldsul unit boxes",
    ),
    thresh=dict(
        longname="--threshold",
        action="store",
        type=float,
        dest="threshold",
        default=None,
        help="Density threshold",
    ),
    dm_only=dict(
        longname="--all-particles",
        action="store_false",
        dest="dm_only",
        default=True,
        help="Use all particles",
    ),
    grids=dict(
        longname="--show-grids",
        action="store_true",
        dest="grids",
        default=False,
        help="Show the grid boundaries",
    ),
    time=dict(
        longname="--time",
        action="store_true",
        dest="time",
        default=False,
        help="Print time in years on image",
    ),
    contours=dict(
        longname="--contours",
        action="store",
        type=int,
        dest="contours",
        default=None,
        help="Number of Contours for Rendering",
    ),
    contour_width=dict(
        longname="--contour_width",
        action="store",
        type=float,
        dest="contour_width",
        default=None,
        help="Width of gaussians used for rendering.",
    ),
    enhance=dict(
        longname="--enhance",
        action="store_true",
        dest="enhance",
        default=False,
        help="Enhance!",
    ),
    valrange=dict(
        short="-r",
        longname="--range",
        action="store",
        type=float,
        dest="valrange",
        default=None,
        nargs=2,
        help="Range, space separated",
    ),
    up=dict(
        longname="--up",
        action="store",
        type=float,
        dest="up",
        default=None,
        nargs=3,
        help="Up, space separated",
    ),
    viewpoint=dict(
        longname="--viewpoint",
        action="store",
        type=float,
        dest="viewpoint",
        default=[1.0, 1.0, 1.0],
        nargs=3,
        help="Viewpoint, space separated",
    ),
    pixels=dict(
        longname="--pixels",
        action="store",
        type=int,
        dest="pixels",
        default=None,
        help="Number of Pixels for Rendering",
    ),
    halos=dict(
        longname="--halos",
        action="store",
        type=str,
        dest="halos",
        default="multiple",
        help="Run halo profiler on a 'single' halo or 'multiple' halos.",
    ),
    halo_radius=dict(
        longname="--halo_radius",
        action="store",
        type=float,
        dest="halo_radius",
        default=0.1,
        help="Constant radius for profiling halos if using hop output files with no "
        + "radius entry. Default: 0.1.",
    ),
    halo_radius_units=dict(
        longname="--halo_radius_units",
        action="store",
        type=str,
        dest="halo_radius_units",
        default="1",
        help="Units for radius used with --halo_radius flag. "
        + "Default: '1' (code units).",
    ),
    halo_hop_style=dict(
        longname="--halo_hop_style",
        action="store",
        type=str,
        dest="halo_hop_style",
        default="new",
        help="Style of hop output file. "
        + "'new' for yt_hop files and 'old' for enzo_hop files.",
    ),
    halo_dataset=dict(
        longname="--halo_dataset",
        action="store",
        type=str,
        dest="halo_dataset",
        default=None,
        help="HaloProfiler dataset.",
    ),
    make_profiles=dict(
        longname="--make_profiles",
        action="store_true",
        default=False,
        help="Make profiles with halo profiler.",
    ),
    make_projections=dict(
        longname="--make_projections",
        action="store_true",
        default=False,
        help="Make projections with halo profiler.",
    ),
)

# This code snippet is modified from Georg Brandl
def bb_apicall(endpoint, data, use_pass=True):
    uri = f"https://api.bitbucket.org/1.0/{endpoint}/"
    # since bitbucket doesn't return the required WWW-Authenticate header when
    # making a request without Authorization, we cannot use the standard urllib2
    # auth handlers; we have to add the requisite header from the start
    if data is not None:
        data = urllib.parse.urlencode(data)
    req = urllib.request.Request(uri, data)
    if use_pass:
        username = input("Bitbucket Username? ")
        password = getpass.getpass()
        upw = f"{username}:{password}"
        req.add_header("Authorization", f"Basic {base64.b64encode(upw).strip()}")
    return urllib.request.urlopen(req).read()


class YTInstInfoCmd(YTCommand):
    name = ["instinfo", "version"]
    args = (
        dict(
            short="-u",
            longname="--update-source",
            action="store_true",
            default=False,
            help="Update the yt installation, if able",
        ),
        dict(
            short="-o",
            longname="--output-version",
            action="store",
            default=None,
            dest="outputfile",
            help="File into which the current revision number will be stored",
        ),
    )
    description = """
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
                update_git(path)
        elif opts.update_source:
            _print_failed_source_update()
        if vstring is not None and opts.outputfile is not None:
            open(opts.outputfile, "w").write(vstring)


class YTLoadCmd(YTCommand):
    name = "load"
    description = """
        Load a single dataset into an IPython instance

        """

    args = ("ds",)

    def __call__(self, args):
        if args.ds is None:
            print("Could not load file.")
            sys.exit()
        import IPython

        import yt
        import yt.mods

        local_ns = yt.mods.__dict__.copy()
        local_ns["ds"] = args.ds
        local_ns["pf"] = args.ds
        local_ns["yt"] = yt

        try:
            from traitlets.config.loader import Config
        except ImportError:
            from IPython.config.loader import Config
        import sys

        cfg = Config()
        # prepend sys.path with current working directory
        sys.path.insert(0, "")
        IPython.embed(config=cfg, user_ns=local_ns)


class YTMapserverCmd(YTCommand):
    args = (
        "proj",
        "field",
        "weight",
        "linear",
        "center",
        "width",
        "cmap",
        dict(
            short="-a",
            longname="--axis",
            action="store",
            type=int,
            dest="axis",
            default=0,
            help="Axis",
        ),
        dict(
            short="-o",
            longname="--host",
            action="store",
            type=str,
            dest="host",
            default=None,
            help="IP Address to bind on",
        ),
        dict(short="ds", nargs=1, type=str, help="The dataset to load."),
    )

    name = "mapserver"
    description = """
        Serve a plot in a GMaps-style interface

        """

    def __call__(self, args):
        from yt.frontends.ramses.data_structures import RAMSESDataset
        from yt.visualization.mapserver.pannable_map import PannableMapServer

        # For RAMSES datasets, use the bbox feature to make the dataset load faster
        if RAMSESDataset._is_valid(args.ds) and args.center and args.width:
            kwa = dict(
                bbox=[
                    [c - args.width / 2 for c in args.center],
                    [c + args.width / 2 for c in args.center],
                ]
            )
        else:
            kwa = dict()

        ds = _fix_ds(args.ds, **kwa)
        if args.center and args.width:
            center = args.center
            width = args.width
            ad = ds.box(
                left_edge=[c - args.width / 2 for c in args.center],
                right_edge=[c + args.width / 2 for c in args.center],
            )
        else:
            center = [0.5] * 3
            width = 1.0
            ad = ds.all_data()

        if args.axis >= 4:
            print("Doesn't work with multiple axes!")
            return
        if args.projection:
            p = ProjectionPlot(
                ds,
                args.axis,
                args.field,
                weight_field=args.weight,
                data_source=ad,
                center=center,
                width=width,
            )
        else:
            p = SlicePlot(
                ds, args.axis, args.field, data_source=ad, center=center, width=width
            )
        p.set_log("all", args.takelog)
        p.set_cmap("all", args.cmap)

        PannableMapServer(p.data_source, args.field, args.takelog, args.cmap)
        try:
            import bottle
        except ImportError as e:
            raise ImportError(
                "The mapserver functionality requires the bottle "
                "package to be installed. Please install using `pip "
                "install bottle`."
            ) from e
        bottle.debug(True)
        if args.host is not None:
            colonpl = args.host.find(":")
            if colonpl >= 0:
                port = int(args.host.split(":")[-1])
                args.host = args.host[:colonpl]
            else:
                port = 8080
            bottle.run(server="auto", host=args.host, port=port)
        else:
            bottle.run(server="auto")


class YTPastebinCmd(YTCommand):
    name = "pastebin"
    args = (
        dict(
            short="-l",
            longname="--language",
            action="store",
            default=None,
            dest="language",
            help="Use syntax highlighter for the file in language",
        ),
        dict(
            short="-L",
            longname="--languages",
            action="store_true",
            default=False,
            dest="languages",
            help="Retrieve a list of supported languages",
        ),
        dict(
            short="-e",
            longname="--encoding",
            action="store",
            default="utf-8",
            dest="encoding",
            help="Specify the encoding of a file (default is "
            "utf-8 or guessing if available)",
        ),
        dict(
            short="-b",
            longname="--open-browser",
            action="store_true",
            default=False,
            dest="open_browser",
            help="Open the paste in a web browser",
        ),
        dict(
            short="-p",
            longname="--private",
            action="store_true",
            default=False,
            dest="private",
            help="Paste as private",
        ),
        dict(
            short="-c",
            longname="--clipboard",
            action="store_true",
            default=False,
            dest="clipboard",
            help="File to output to; else, print.",
        ),
        dict(short="file", type=str),
    )
    description = """
        Post a script to an anonymous pastebin

        """

    def __call__(self, args):
        from yt.utilities import lodgeit as lo

        lo.main(
            args.file,
            languages=args.languages,
            language=args.language,
            encoding=args.encoding,
            open_browser=args.open_browser,
            private=args.private,
            clipboard=args.clipboard,
        )


class YTPastebinGrabCmd(YTCommand):
    args = (dict(short="number", type=str),)
    name = "pastebin_grab"
    description = """
        Print an online pastebin to STDOUT for local use.
        """

    def __call__(self, args):
        from yt.utilities import lodgeit as lo

        lo.main(None, download=args.number)


class YTPlotCmd(YTCommand):
    args = (
        "width",
        "unit",
        "bn",
        "proj",
        "center",
        "zlim",
        "axis",
        "field",
        "weight",
        "skip",
        "cmap",
        "output",
        "grids",
        "time",
        "ds",
        "max",
        "log",
        "linear",
        dict(
            short="-fu",
            longname="--field-unit",
            action="store",
            type=str,
            dest="field_unit",
            default=None,
            help="Desired field units",
        ),
        dict(
            longname="--show-scale-bar",
            action="store_true",
            help="Annotate the plot with the scale",
        ),
    )

    name = "plot"

    description = """
        Create a set of images

        """

    def __call__(self, args):
        ds = args.ds
        center = args.center
        if args.center == (-1, -1, -1):
            mylog.info("No center fed in; seeking.")
            v, center = ds.find_max("density")
        if args.max:
            v, center = ds.find_max("density")
        elif args.center is None:
            center = 0.5 * (ds.domain_left_edge + ds.domain_right_edge)
        center = np.array(center)
        if ds.dimensionality < 3:
            dummy_dimensions = np.nonzero(ds.index.grids[0].ActiveDimensions <= 1)
            axes = dummy_dimensions[0][0]
        elif args.axis == 4:
            axes = range(3)
        else:
            axes = args.axis

        unit = args.unit
        if unit is None:
            unit = "unitary"
        if args.width is None:
            width = None
        else:
            width = (args.width, args.unit)

        for ax in always_iterable(axes):
            mylog.info("Adding plot for axis %i", ax)
            if args.projection:
                plt = ProjectionPlot(
                    ds,
                    ax,
                    args.field,
                    center=center,
                    width=width,
                    weight_field=args.weight,
                )
            else:
                plt = SlicePlot(ds, ax, args.field, center=center, width=width)
            if args.grids:
                plt.annotate_grids()
            if args.time:
                plt.annotate_timestamp()
            if args.show_scale_bar:
                plt.annotate_scale()

            if args.field_unit:
                plt.set_unit(args.field, args.field_unit)

            plt.set_cmap(args.field, args.cmap)
            plt.set_log(args.field, args.takelog)
            if args.zlim:
                plt.set_zlim(args.field, *args.zlim)
            ensure_dir_exists(args.output)
            plt.save(os.path.join(args.output, f"{ds}"))


class YTRPDBCmd(YTCommand):
    name = "rpdb"
    description = """
        Connect to a currently running (on localhost) rpd session.

        Commands run with --rpdb will trigger an rpdb session with any
        uncaught exceptions.

        """
    args = (
        dict(
            short="-t",
            longname="--task",
            action="store",
            default=0,
            dest="task",
            help="Open a web browser.",
        ),
    )

    def __call__(self, args):
        from . import rpdb

        rpdb.run_rpdb(int(args.task))


class YTNotebookCmd(YTCommand):
    name = ["notebook"]
    args = (
        dict(
            short="-o",
            longname="--open-browser",
            action="store_true",
            default=False,
            dest="open_browser",
            help="Open a web browser.",
        ),
        dict(
            short="-p",
            longname="--port",
            action="store",
            default=0,
            dest="port",
            help="Port to listen on; defaults to auto-detection.",
        ),
        dict(
            short="-prof",
            longname="--profile",
            action="store",
            default=None,
            dest="profile",
            help="The IPython profile to use when launching the kernel.",
        ),
        dict(
            short="-n",
            longname="--no-password",
            action="store_true",
            default=False,
            dest="no_password",
            help="If set, do not prompt or use a password.",
        ),
    )
    description = """
        Start the Jupyter Notebook locally.
        """

    def __call__(self, args):
        kwargs = {}
        from notebook.notebookapp import NotebookApp

        print(
            "You must choose a password so that others cannot connect to "
            "your notebook."
        )
        pw = ytcfg.get("yt", "notebook_password")
        if len(pw) == 0 and not args.no_password:
            import IPython.lib

            pw = IPython.lib.passwd()
            print("If you would like to use this password in the future,")
            print("place a line like this inside the [yt] section in your")
            print("yt configuration file at ~/.config/yt/yt.toml")
            print()
            print(f"notebook_password = {pw}")
            print()
        elif args.no_password:
            pw = None
        if args.port != 0:
            kwargs["port"] = int(args.port)
        if args.profile is not None:
            kwargs["profile"] = args.profile
        if pw is not None:
            kwargs["password"] = pw
        app = NotebookApp(open_browser=args.open_browser, **kwargs)
        app.initialize(argv=[])
        print()
        print("***************************************************************")
        print()
        print("The notebook is now live at:")
        print()
        print(f"     http://127.0.0.1:{app.port}/")
        print()
        print("Recall you can create a new SSH tunnel dynamically by pressing")
        print(f"~C and then typing -L{app.port}:localhost:{app.port}")
        print("where the first number is the port on your local machine. ")
        print()
        print(
            "If you are using %s on your machine already, try "
            "-L8889:localhost:%s" % (app.port, app.port)
        )
        print()
        print("***************************************************************")
        print()
        app.start()


class YTStatsCmd(YTCommand):
    args = (
        "outputfn",
        "bn",
        "skip",
        "ds",
        "field",
        dict(
            longname="--max",
            action="store_true",
            default=False,
            dest="max",
            help="Display maximum of field requested through -f option.",
        ),
        dict(
            longname="--min",
            action="store_true",
            default=False,
            dest="min",
            help="Display minimum of field requested through -f option.",
        ),
    )
    name = "stats"
    description = """
        Print stats and max/min value of a given field (if requested),
        for one or more datasets

        (default field is density)

        """

    def __call__(self, args):
        ds = args.ds
        ds.print_stats()
        vals = {}
        field = ds._get_field_info(args.field)
        if args.max:
            vals["max"] = ds.find_max(field)
            print(f"Maximum {field.name}: {vals['max'][0]:0.5e} at {vals['max'][1]}")
        if args.min:
            vals["min"] = ds.find_min(field)
            print(f"Minimum {field.name}: {vals['min'][0]:0.5e} at {vals['min'][1]}")
        if args.output is not None:
            t = ds.current_time.to("yr")
            with open(args.output, "a") as f:
                f.write(f"{ds} ({t:0.5e})\n")
                if "min" in vals:
                    f.write(
                        "Minimum %s is %0.5e at %s\n"
                        % (field.name, vals["min"][0], vals["min"][1])
                    )
                if "max" in vals:
                    f.write(
                        "Maximum %s is %0.5e at %s\n"
                        % (field.name, vals["max"][0], vals["max"][1])
                    )


class YTUpdateCmd(YTCommand):
    args = ("all",)
    name = "update"
    description = """
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
            update_git(path)
        else:
            _print_failed_source_update(opts.reinstall)


class YTDeleteImageCmd(YTCommand):
    args = (dict(short="delete_hash", type=str),)
    description = """
        Delete image from imgur.com.

        """
    name = "delete_image"

    def __call__(self, args):
        headers = {"Authorization": f"Client-ID {ytcfg.get('yt', 'imagebin_api_key')}"}

        delete_url = ytcfg.get("yt", "imagebin_delete_url")
        req = urllib.request.Request(
            delete_url.format(delete_hash=args.delete_hash),
            headers=headers,
            method="DELETE",
        )
        try:
            response = urllib.request.urlopen(req).read().decode()
        except urllib.error.HTTPError as e:
            print("ERROR", e)
            return {"deleted": False}

        rv = json.loads(response)
        if "success" in rv and rv["success"]:
            print("\nImage successfully deleted!\n")
        else:
            print()
            print("Something has gone wrong!  Here is the server response:")
            print()
            pprint.pprint(rv)


class YTUploadImageCmd(YTCommand):
    args = (dict(short="file", type=str),)
    description = """
        Upload an image to imgur.com.  Must be PNG.

        """
    name = "upload_image"

    def __call__(self, args):
        filename = args.file
        if not filename.endswith(".png"):
            print("File must be a PNG file!")
            return 1
        headers = {"Authorization": f"Client-ID {ytcfg.get('yt', 'imagebin_api_key')}"}

        image_data = base64.b64encode(open(filename, "rb").read())
        parameters = {
            "image": image_data,
            type: "base64",
            "name": filename,
            "title": f"{filename} uploaded by yt",
        }
        data = urllib.parse.urlencode(parameters).encode("utf-8")
        req = urllib.request.Request(
            ytcfg.get("yt", "imagebin_upload_url"), data=data, headers=headers
        )
        try:
            response = urllib.request.urlopen(req).read().decode()
        except urllib.error.HTTPError as e:
            print("ERROR", e)
            return {"uploaded": False}
        rv = json.loads(response)
        if "data" in rv and "link" in rv["data"]:
            print()
            print("Image successfully uploaded!  You can find it at:")
            print(f"    {rv['data']['link']}")
            print()
            print("If you'd like to delete it, use the following")
            print(f"    yt delete_image {rv['data']['deletehash']}")
            print()
        else:
            print()
            print("Something has gone wrong!  Here is the server response:")
            print()
            pprint.pprint(rv)


class YTUploadFileCmd(YTCommand):
    args = (dict(short="file", type=str),)
    description = """
        Upload a file to yt's curldrop.

        """
    name = "upload"

    def __call__(self, args):
        from yt.utilities.on_demand_imports import _requests as requests

        fs = iter(FileStreamer(open(args.file, "rb")))
        upload_url = ytcfg.get("yt", "curldrop_upload_url")
        r = requests.put(upload_url + "/" + os.path.basename(args.file), data=fs)
        print()
        print(r.text)


class YTConfigLocalConfigHandler:
    def load_config(self, args):
        import os

        from yt.config import YTConfig
        from yt.utilities.configure import CONFIG

        local_config_file = YTConfig.get_local_config_file()
        global_config_file = YTConfig.get_global_config_file()

        local_exists = os.path.exists(local_config_file)
        global_exists = os.path.exists(global_config_file)

        local_arg_exists = hasattr(args, "local")
        global_arg_exists = hasattr(args, "global")

        if getattr(args, "local", False):
            config_file = local_config_file
        elif getattr(args, "global", False):
            config_file = global_config_file
        else:
            config_file: Optional[str] = None
            if local_exists and global_exists:
                s = (
                    "Yt detected a local and a global configuration file, refusing "
                    "to proceed.\n"
                    f"Local config file: {local_config_file}\n"
                    f"Global config file: {global_config_file}"
                )
                # Only print the info about "--global" and "--local" if they exist
                if local_arg_exists and global_arg_exists:
                    s += (
                        "\n"  # missing eol from previous string
                        "Specify which one you want to use using the `--local` or the "
                        "`--global` flags."
                    )
                sys.exit(s)
            elif local_exists:
                config_file = local_config_file
            elif global_exists:
                config_file = global_config_file

            if config_file is None:
                print("WARNING: no configuration file installed.", file=sys.stderr)
            else:
                print(
                    f"INFO: reading configuration file: {config_file}", file=sys.stderr
                )
        CONFIG.read(config_file)

        self.config_file = config_file


_global_local_args = (
    dict(
        short="--local",
        action="store_true",
        help="Store the configuration in the local configuration file.",
    ),
    dict(
        short="--global",
        action="store_true",
        help="Store the configuration in the global configuration file.",
    ),
)


class YTConfigGetCmd(YTCommand, YTConfigLocalConfigHandler):
    subparser = "config"
    name = "get"
    description = "get a config value"
    args = (
        dict(short="section", help="The section containing the option."),
        dict(short="option", help="The option to retrieve."),
        *_global_local_args,
    )

    def __call__(self, args):
        from yt.utilities.configure import get_config

        self.load_config(args)

        print(get_config(args.section, args.option))


class YTConfigSetCmd(YTCommand, YTConfigLocalConfigHandler):
    subparser = "config"
    name = "set"
    description = "set a config value"
    args = (
        dict(short="section", help="The section containing the option."),
        dict(short="option", help="The option to set."),
        dict(short="value", help="The value to set the option to."),
        *_global_local_args,
    )

    def __call__(self, args):
        from yt.utilities.configure import set_config

        self.load_config(args)
        if self.config_file is None:
            self.config_file = os.path.join(os.getcwd(), "yt.toml")
            print(
                f"INFO: configuration will be written to {self.config_file}",
                file=sys.stderr,
            )
        set_config(args.section, args.option, args.value, self.config_file)


class YTConfigRemoveCmd(YTCommand, YTConfigLocalConfigHandler):
    subparser = "config"
    name = "rm"
    description = "remove a config option"
    args = (
        dict(short="section", help="The section containing the option."),
        dict(short="option", help="The option to remove."),
        *_global_local_args,
    )

    def __call__(self, args):
        from yt.utilities.configure import rm_config

        self.load_config(args)

        rm_config(args.section, args.option, self.config_file)


class YTConfigListCmd(YTCommand, YTConfigLocalConfigHandler):
    subparser = "config"
    name = "list"
    description = "show the config content"
    args = _global_local_args

    def __call__(self, args):
        from yt.utilities.configure import write_config

        self.load_config(args)

        write_config(sys.stdout)


class YTConfigPrintPath(YTCommand, YTConfigLocalConfigHandler):
    subparser = "config"
    name = "print-path"
    description = "show path to the config file"
    args = _global_local_args

    def __call__(self, args):
        self.load_config(args)

        print(self.config_file)


class YTSearchCmd(YTCommand):
    args = (
        dict(
            short="-o",
            longname="--output",
            action="store",
            type=str,
            dest="output",
            default="yt_index.json",
            help="File in which to place output",
        ),
        dict(
            longname="--check-all",
            short="-a",
            help="Attempt to load every file",
            action="store_true",
            default=False,
            dest="check_all",
        ),
        dict(
            longname="--full",
            short="-f",
            help="Output full contents of parameter file",
            action="store_true",
            default=False,
            dest="full_output",
        ),
    )
    description = """
        Attempt to find outputs that yt can recognize in directories.
        """
    name = "search"

    def __call__(self, args):
        from yt.utilities.object_registries import output_type_registry

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
            except YTUnidentifiedDataType:
                continue
            records.append(record)
        with open(args.output, "w") as f:
            json.dump(records, f, indent=4)
        print(f"Identified {len(records)} records output to {args.output}")


class YTDownloadData(YTCommand):

    args = (
        dict(
            short="filename",
            action="store",
            type=str,
            help="The name of the file to download",
            nargs="?",
            default="",
        ),
        dict(
            short="location",
            action="store",
            type=str,
            nargs="?",
            help="The location in which to place the file, can be "
            '"supp_data_dir", "test_data_dir", or any valid '
            "path on disk. ",
            default="",
        ),
        dict(
            longname="--overwrite",
            short="-c",
            help="Overwrite existing file.",
            action="store_true",
            default=False,
        ),
        dict(
            longname="--list",
            short="-l",
            help="Display all available files.",
            action="store_true",
            default=False,
        ),
    )
    description = """
        Download a file from http://yt-project.org/data and save it to a
        particular location. Files can be saved to the locations provided
        by the "test_data_dir" or "supp_data_dir" configuration entries, or
        any valid path to a location on disk.
        """
    name = "download"

    def __call__(self, args):
        if args.list:
            self.get_list()
            return
        if not args.filename:
            raise RuntimeError(
                "You need to provide a filename. See --help "
                "for details or use --list to get available "
                "datasets."
            )
        elif not args.location:
            raise RuntimeError(
                "You need to specify download location. See --help for details."
            )
        data_url = f"http://yt-project.org/data/{args.filename}"
        if args.location in ["test_data_dir", "supp_data_dir"]:
            data_dir = ytcfg.get("yt", args.location)
            if data_dir == "/does/not/exist":
                raise RuntimeError(f"'{args.location}' is not configured!")
        else:
            data_dir = args.location
        if not os.path.exists(data_dir):
            print(f"The directory '{data_dir}' does not exist. Creating...")
            ensure_dir(data_dir)
        data_file = os.path.join(data_dir, args.filename)
        if os.path.exists(data_file) and not args.overwrite:
            raise OSError(f"File '{data_file}' exists and overwrite=False!")
        print(f"Attempting to download file: {args.filename}")
        fn = download_file(data_url, data_file)

        if not os.path.exists(fn):
            raise OSError(f"The file '{args.filename}' did not download!!")
        print(f"File: {args.filename} downloaded successfully to {data_file}")

    def get_list(self):
        data = (
            urllib.request.urlopen("http://yt-project.org/data/datafiles.json")
            .read()
            .decode("utf8")
        )
        data = json.loads(data)
        for key in data:
            for ds in data[key]:
                ds["fullname"] = ds["url"].replace("http://yt-project.org/data/", "")
                print("{fullname} ({size}) type: {code}".format(**ds))
                for line in textwrap.wrap(ds["description"]):
                    print("\t", line)


def run_main():
    args = parser.parse_args()
    # The following is a workaround for a nasty Python 3 bug:
    # http://bugs.python.org/issue16308
    # http://bugs.python.org/issue9253
    try:
        args.func
    except AttributeError:
        parser.print_help()
        sys.exit(0)

    args.func(args)


if __name__ == "__main__":
    run_main()
