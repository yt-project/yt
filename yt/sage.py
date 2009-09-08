"""
Convenience functions for operation inside a SAGE notebook

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2009 Matthew Turk.  All Rights Reserved.

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

import os

if "SAGE_ROOT" not in os.environ: raise ImportError

from yt.config import ytcfg

from sage.server.support import EMBEDDED_MODE
if not EMBEDDED_MODE: raise ImportError

ytcfg["yt","__insagenotebook"] = "True"

def width_slider(plot_obj):
    """
    This function accepts a plot object (*plot_obj*) and generates a width
    slider that will automatically save when released.
    """
    import sage.server.notebook.interact as interact

    if hasattr(plot_obj, 'save'):
        sfunc = plot_obj.save
    elif hasattr(plot_obj, 'save_image'):
        sfunc = plot_obj.save_image
    else:
        raise RuntimeError

    @interact.interact
    def do_slider(width=interact.slider(-3.0, 0.0, 
        plot_obj.set_width(10**width, 'unitary')
        sfunc("temp")


