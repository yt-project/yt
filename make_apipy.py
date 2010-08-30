ss = """\"\"\"
API for %s

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Author: J.S. Oishi <jsoishi@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Author: Britton Smith <brittonsmith@gmail.com>
Affiliation: MSU
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

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

\"\"\"
"""

import sys, os
n = sys.argv[-1]

fn = os.path.join(os.path.abspath(n), "api.py")
#if os.path.exists(fn):
    #print "%s exists!  Not overwriting." % fn
bn = n.replace("/",".")

print "%s -> %s" % (bn, fn)

open(fn, "w").write(ss % (bn))
