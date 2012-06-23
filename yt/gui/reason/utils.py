"""
Utilities for Reason

Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: Michigan State University
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt-project.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

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

from .bottle_mods import PayloadHandler
import base64
import types

def load_script(filename):
    contents = open(filename).read()
    payload_handler = PayloadHandler()
    payload_handler.add_payload(
        {'type': 'cell_contents',
         'value': contents}
    )
    return

def deliver_image(im):
    if hasattr(im, 'read'):
        img_data = base64.b64encode(im.read())
    elif isinstance(im, types.StringTypes) and \
         im.endswith(".png"):
        img_data = base64.b64encode(open(im).read())
    elif isinstance(im, types.StringTypes):
        img_data = im
    else:
        raise RuntimeError
    ph = PayloadHandler()
    payload = {'type':'png_string',
               'image_data':img_data}
    ph.add_payload(payload)

