"""
This module handles both the interactive component as well as the command-line
scripts.

@deprecated: IPython is the wave of the future, people!
@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

#from yt.fido import *  # Not sure we want this
import sys
import glob, os.path

from yt import ytcfg

from code import InteractiveConsole

ytBannerDefault="Welcome to yt!  For help, type ytHelp()"

def ytHelp():
    s = """
    This is a fully-functional Python prompt.  You should already have access
    to the modules lagos, raven, enki and fido, and the functions selectRun and
    selectPF.  Additionally, any run you have already fetched will be available
    as myRun.  (Type dir() to see what's available.)
    """
    print s

def runyt(my_dict = {}, ytBanner = ytBannerDefault):
    """
    This function spawns a shell for yt, with some modules already imported
    @keyword my_dict: namespace dictionary
    @type my_dict: dict
    """
    fs = ytShell(my_dict)
    pythonstartup = os.getenv("PYTHONSTARTUP")
    if pythonstartup:
        fs.runsource("execfile('%s')" % (pythonstartup), pythonstartup)
    if ytcfg.has_option("yt","ytrc"):
        # This is *highly* unsafe.  I think.
        ytrc = ytcfg.get("yt","ytrc")
        if os.path.exists(ytrc):
            fs.runsource("execfile('%s')" % (ytrc), ytrc)
    fs.interact(ytBanner)

class ytShell(InteractiveConsole):
    def __init__(self, my_dict = {}):
        import yt.fido, yt.lagos, yt.raven, yt.enki
        my_dict.update({'fido':yt.fido, \
                        'lagos':yt.lagos, \
                        'raven':yt.raven, \
                        'enki':yt.enki, \
                        'ytHelp':ytHelp, \
                        'selectRun':yt.fido.selectRun, \
                        'selectPF':yt.fido.selectPF})
        sys.ps1 = "yt> "
        try:
            import readline
        except ImportError:
            print "Module readline not available."
        else:
            import rlcompleter
            readline.parse_and_bind("tab: complete")
        # Put all of our readline completer actions here, I s'pose
        InteractiveConsole.__init__(self, locals=my_dict)
