"""
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: NSF / Columbia
Homepage: http://yt.enzotools.org/
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

from mercurial import ui, repo, commands, hg
import json
import os
import time
import uuid

from yt.config import ytcfg

class PostInventory(object):
    def __init__(self, uu = None, repo_fn = None):
        if uu is None: uu = ui.ui()
        if repo_fn is None: repo_fn = ytcfg.get("yt","pasteboard_repo")
        if repo_fn == '':
            raise KeyError("~/.yt/config:[yt]pasteboard_repo")
        self.repo_fn = repo_fn
        self.bbrepo = hg.repository(uu, repo_fn)
        config_fn = os.path.join(repo_fn, ".hg", "hgrc")
        uu.readconfig(config_fn)
        commands.pull(uu, self.bbrepo)
        commands.update(uu, self.bbrepo, clean=True)
        self.uu = uu

    def regenerate_posts(self):
        self.posts = []
        for file in self.bbrepo["tip"]:
            if file.startswith("posts/") and file.count("/") == 1:
                filectx = self.bbrepo["tip"][file]
                last_mod = filectx.filectx(filectx.filerev()).date()
                self.posts.append((last_mod[0] + last_mod[1], file))
        self.posts.sort()
        self.posts = self.posts[::-1]

    def add_post(self, filename, uu = None,
                 highlight = True, push = True):
        # We assume the post filename exists in the current space
        self.regenerate_posts()
        if uu is None: uu = self.uu
        prefix = uuid.uuid4()
        name = "%s-%s" % (prefix, os.path.basename(filename))
        name_noext = name.replace(".","-")
        hfn = "html/%s.html" % (name_noext)
        pfn = "posts/%s" % (name)
        abs_pfn = os.path.join(self.repo_fn, pfn)
        abs_hfn = os.path.join(self.repo_fn, hfn)
        self.posts.insert(0, (int(time.time()), "posts/%s" % name))
        if not os.path.exists(abs_pfn):
            open(abs_pfn,"w").write(open(filename).read())
        inv_fname = self.update_inventory()
        if highlight and not name.endswith(".html"):
            from pygments.cmdline import main as pygmain
            rv = pygmain(["pygmentize", "-o", abs_hfn,
                          "-O", "full", abs_pfn])
        if not highlight or rv:
            content = open(abs_pfn).read()
            open(abs_hfn, "w").write(
                "<HTML><BODY><PRE>" + content + "</PRE></BODY></HTML>")
        commands.add(uu, self.bbrepo, abs_pfn, abs_hfn)
        commands.commit(uu, self.bbrepo, abs_hfn, abs_pfn,
                        inv_fname, message="Adding %s" % name)
        if push: commands.push(uu, self.bbrepo)

    def update_inventory(self):
        tip = self.bbrepo["tip"]
        vals = []
        for t, pfn in self.posts:
            if pfn not in tip:
                d = open(os.path.join(self.repo_fn, pfn)).read()
            else:
                d = tip[pfn].data()
            if len(d) > 80: d = d[:77] + "..."
            name_noext = pfn[6:].replace(".","-")
            vals.append(dict(modified = time.ctime(t),
                             modtime = t,
                             fullname = pfn,
                             htmlname = "html/%s.html" % name_noext,
                             name = pfn[43:], # 6 for posts/ then 36 for UUID
                             descr = d)) 
        fn = os.path.join(self.repo_fn, "inventory.json")
        f = open(fn, "w")
        f.write("var inventory_data = ")
        json.dump(vals, f)
        f.write(";")
        return fn
