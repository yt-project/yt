# This custom importer for yt will set up some IPython notebook-specific
# helpers.  For instance, it will eventually add items to the menubar.

import __builtin__
if "__IPYTHON__" not in dir(__builtin__):
    raise ImportError

from IPython.core.interactiveshell import InteractiveShell
from IPython.core.display import display, display_html
inst = InteractiveShell.instance()
ip = inst.get_ipython()
ip.enable_pylab("inline", import_all=False)

from yt.config import ytcfg
ytcfg["yt", "ipython_notebook"] = "True"

from yt.mods import *
