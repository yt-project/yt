# This extension is quite simple:
#  1. It accepts a script name
#  2. This script is added to the document in a literalinclude
#  3. Any _static images found will be added

import glob
import os
import shutil

from docutils.parsers.rst import Directive, directives

# Some of this magic comes from the matplotlib plot_directive.


def setup(app):
    app.add_directive("yt_colormaps", ColormapScript)
    setup.app = app
    setup.config = app.config
    setup.confdir = app.confdir

    retdict = dict(version="0.1", parallel_read_safe=True, parallel_write_safe=True)

    return retdict


class ColormapScript(Directive):
    required_arguments = 1
    optional_arguments = 0

    def run(self):
        rst_file = self.state_machine.document.attributes["source"]
        rst_dir = os.path.abspath(os.path.dirname(rst_file))
        script_fn = directives.path(self.arguments[0])
        script_bn = os.path.basename(script_fn)

        # This magic is from matplotlib
        dest_dir = os.path.abspath(
            os.path.join(setup.app.builder.outdir, os.path.dirname(script_fn))
        )
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)  # no problem here for me, but just use built-ins

        rel_dir = os.path.relpath(rst_dir, setup.confdir)
        place = os.path.join(dest_dir, rel_dir)
        if not os.path.isdir(place):
            os.makedirs(place)
        shutil.copyfile(
            os.path.join(rst_dir, script_fn), os.path.join(place, script_bn)
        )

        im_path = os.path.join(rst_dir, "_static")
        images = sorted(glob.glob(os.path.join(im_path, "*.png")))
        lines = []
        for im in images:
            im_name = os.path.join("_static", os.path.basename(im))
            lines.append(f".. image:: {im_name}")
            lines.append("   :width: 400")
            lines.append(f"   :target: ../../_images/{os.path.basename(im)}")
            lines.append("\n")
        lines.append("\n")
        self.state_machine.insert_input(lines, rst_file)
        return []
