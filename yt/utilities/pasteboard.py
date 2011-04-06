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
import urllib

from yt.config import ytcfg

def _get_last_mod(filectx):
    rev = filectx.filectx(filectx.filerev())
    return rev

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
        if not os.path.exists(os.path.join(repo_fn, "posts")):
            os.makedirs(os.path.join(repo_fn, "posts"))
        if not os.path.exists(os.path.join(repo_fn, "html")):
            os.makedirs(os.path.join(repo_fn, "html"))
        self.uu = uu

    def regenerate_posts(self):
        self.posts = []
        for file in self.bbrepo["tip"]:
            if file.startswith("posts/") and file.count("/") == 1 \
               and not file.endswith(".desc"):
                filectx = self.bbrepo["tip"][file]
                last_mod = _get_last_mod(filectx).date()
                self.posts.append((last_mod[0] + last_mod[1], file))
        self.posts.sort()
        self.posts = self.posts[::-1]

    def add_post(self, filename, desc = None,
                 uu = None, highlight = True, push = True):
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
        if desc is not None:
            open(abs_pfn + ".desc", "w").write(desc)
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
        to_manage = [abs_pfn, abs_hfn]
        if desc is not None: to_manage.append(abs_pfn + ".desc")
        commands.add(uu, self.bbrepo, *to_manage)
        commands.commit(uu, self.bbrepo, *(to_manage + [inv_fname]),
                        message="Adding %s" % name)
        if push: commands.push(uu, self.bbrepo)

    def update_inventory(self):
        tip = self.bbrepo["tip"]
        vals = []
        for t, pfn in self.posts:
            dfn = pfn + ".desc"
            if dfn in tip:
                d = tip[dfn].data()
                last_mod =_get_last_mod(tip[dfn])
                last_hash = last_mod.hex()
                uname = last_mod.user()
            elif pfn not in tip:
                abs_pfn = os.path.join(self.repo_fn, pfn)
                uname = self.uu.config("ui","username")
                if os.path.exists(abs_pfn + ".desc"):
                    d = open(abs_pfn + ".desc").read()
                else:
                    d = open(abs_pfn).read()
                last_hash = "tip"
            else:
                d = tip[pfn].data()
                last_mod = _get_last_mod(tip[pfn])
                last_hash = last_mod.hex()
                uname = last_mod.user()
            if len(d) > 80: d = d[:77] + "..."
            name_noext = pfn[6:].replace(".","-")
            vals.append(dict(modified = time.ctime(t),
                             modtime = t,
                             lastmod_hash = last_hash,
                             fullname = pfn,
                             htmlname = "html/%s.html" % name_noext,
                             name = pfn[43:], # 6 for posts/ then 36 for UUID
                             username = uname,
                             descr = d)) 
        fn = os.path.join(self.repo_fn, "inventory.json")
        f = open(fn, "w")
        f.write("var inventory_data = ")
        json.dump(vals, f, indent = 1)
        f.write(";")
        return fn

def retrieve_pastefile(username, paste_id, output_fn = None):
    # First we get the username's inventory.json
    s = urllib.urlopen("http://%s.bitbucket.org/inventory.json" % (username))
    data = s.read()
    # This is an ugly, ugly hack for my lack of understanding of how best to
    # handle this JSON stuff.
    data = data[data.find("=")+1:data.rfind(";")] 
    #import pdb;pdb.set_trace()
    inv = json.loads(data)
    k = None
    if len(paste_id) == 36:
        # Then this is a UUID
        for k in inv:
            if k['fullname'][6:42] == paste_id: break
    elif len(paste_id) == 10:
        pp = int(paste_id)
        for k in inv:
            if k['modtime'] == pp: break
    if k is None: raise KeyError(k)
    # k is our key
    url = "http://%s.bitbucket.org/%s" % (username, k['fullname'])
    s = urllib.urlopen(url)
    data = s.read()
    if output_fn is not None:
        if os.path.exists(output_fn): raise IOError(output_fn)
        open(output_fn, "w").write(data)
    else:
        print data
