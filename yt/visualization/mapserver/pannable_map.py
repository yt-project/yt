import os
from functools import wraps

import bottle
import numpy as np

from yt.fields.derived_field import ValidateSpatial
from yt.utilities.lib.misc_utilities import get_color_bounds
from yt.utilities.png_writer import write_png_to_string
from yt.visualization.fixed_resolution import FixedResolutionBuffer
from yt.visualization.image_writer import apply_colormap

local_dir = os.path.dirname(__file__)


def exc_writeout(f):
    import traceback

    @wraps(f)
    def func(*args, **kwargs):
        try:
            rv = f(*args, **kwargs)
            return rv
        except Exception:
            traceback.print_exc(None, open("temp.exc", "w"))
            raise

    return func


class PannableMapServer:
    _widget_name = "pannable_map"

    def __init__(self, data, field, takelog, cmap, route_prefix=""):
        self.data = data
        self.ds = data.ds
        self.field = field
        self.cmap = cmap

        bottle.route(f"{route_prefix}/map/:field/:L/:x/:y.png")(self.map)
        bottle.route(f"{route_prefix}/map/:field/:L/:x/:y.png")(self.map)
        bottle.route(f"{route_prefix}/")(self.index)
        bottle.route(f"{route_prefix}/:field")(self.index)
        bottle.route(f"{route_prefix}/index.html")(self.index)
        bottle.route(f"{route_prefix}/list", "GET")(self.list_fields)
        # This is a double-check, since we do not always mandate this for
        # slices:
        self.data[self.field] = self.data[self.field].astype("float64")
        bottle.route(f"{route_prefix}/static/:path", "GET")(self.static)

        self.takelog = takelog
        self._lock = False

        for unit in ["Gpc", "Mpc", "kpc", "pc"]:
            v = self.ds.domain_width[0].in_units(unit).value
            if v > 1:
                break
        self.unit = unit
        self.px2unit = self.ds.domain_width[0].in_units(unit).value / 256

    def lock(self):
        import time

        while self._lock:
            time.sleep(0.01)
        self._lock = True

    def unlock(self):
        self._lock = False

    def map(self, field, L, x, y):
        if "," in field:
            field = tuple(field.split(","))
        cmap = self.cmap
        dd = 1.0 / (2.0 ** (int(L)))
        relx = int(x) * dd
        rely = int(y) * dd
        DW = self.ds.domain_right_edge - self.ds.domain_left_edge
        xl = self.ds.domain_left_edge[0] + relx * DW[0]
        yl = self.ds.domain_left_edge[1] + rely * DW[1]
        xr = xl + dd * DW[0]
        yr = yl + dd * DW[1]
        try:
            self.lock()
            w = 256  # pixels
            data = self.data[field]
            frb = FixedResolutionBuffer(self.data, (xl, xr, yl, yr), (w, w))
            cmi, cma = get_color_bounds(
                self.data["px"],
                self.data["py"],
                self.data["pdx"],
                self.data["pdy"],
                data,
                self.ds.domain_left_edge[0],
                self.ds.domain_right_edge[0],
                self.ds.domain_left_edge[1],
                self.ds.domain_right_edge[1],
                dd * DW[0] / (64 * 256),
                dd * DW[0],
            )
        finally:
            self.unlock()

        if self.takelog:
            cmi = np.log10(cmi)
            cma = np.log10(cma)
            to_plot = apply_colormap(
                np.log10(frb[field]), color_bounds=(cmi, cma), cmap_name=cmap
            )
        else:
            to_plot = apply_colormap(
                frb[field], color_bounds=(cmi, cma), cmap_name=cmap
            )

        rv = write_png_to_string(to_plot)
        return rv

    def index(self, field=None):
        if field is not None:
            self.field = field
        return bottle.static_file(
            "map_index.html", root=os.path.join(local_dir, "html")
        )

    def static(self, path):
        if path[-4:].lower() in (".png", ".gif", ".jpg"):
            bottle.response.headers["Content-Type"] = f"image/{path[-3:].lower()}"
        elif path[-4:].lower() == ".css":
            bottle.response.headers["Content-Type"] = "text/css"
        elif path[-3:].lower() == ".js":
            bottle.response.headers["Content-Type"] = "text/javascript"
        full_path = os.path.join(os.path.join(local_dir, "html"), path)
        return open(full_path).read()

    def list_fields(self):
        d = {}

        # Add fluid fields (only gas for now)
        for ftype in self.ds.fluid_types:
            d[ftype] = []
            for f in self.ds.derived_field_list:
                if f[0] != ftype:
                    continue
                # Discard fields which need ghost zones for now
                df = self.ds.field_info[f]
                if any(isinstance(v, ValidateSpatial) for v in df.validators):
                    continue
                # Discard cutting plane fields
                if "cutting" in f[1]:
                    continue
                active = f[1] == self.field
                d[ftype].append((f, active))

        print(self.px2unit, self.unit)
        return dict(data=d, px2unit=self.px2unit, unit=self.unit, active=self.field)
