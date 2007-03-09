"""
This module handles both the interactive component as well as the command-line
scripts.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

#from yt.fido import *  # Not sure we want this
import dialog, sys
import glob, os.path

import yt.lagos as lagos
from yt.fido import *

from code import InteractiveConsole

FidoBanner="Welcome to Fido!  For help, type fidoHelp()"

def fidoHelp():
    s = """
    This is a fully-functional Python prompt.  You should already have access
    to the modules lagos, raven, enki and fido, and the function selectRun.
    Additionally, any run you have already fetched will be available as myRun.
    """
    print s

def runBranch():
    if len(sys.argv) <= 2:
        print "You need to give me some more info, chuckles!"
        print "Like, what's the new metadata string, and which file do you want to split at?"
        sys.exit()
    new_md = sys.argv[1]
    pf = guessPF(sys.argv[2])
    guess = os.path.basename(os.getcwd())
    thisRun = fetchRun(guess)
    new_dir = ytcfg.get("Fido","WorkDir")
    branch(thisRun, pf, new_md, os.path.join(new_dir, new_md))
    sys.exit()

def runDigup():
    if len(sys.argv) <= 1:
        print "You need to give me some more info, chuckles!"
        print "Like, which is the parameter file you want to dig up?"
        sys.exit()
    new_pf = guessPF(sys.argv[1])
    guess = os.path.basename(os.getcwd())
    digUp(guess, new_pf)
    sys.exit()

def runRevert():
    if len(sys.argv) <= 1:
        print "You need to give me some more info, chuckles!"
        print "Like, which is the latest parameter file you want to keep?"
        sys.exit()
    new_pf = guessPF(sys.argv[1])
    guess = os.path.basename(os.getcwd())
    revert(guess, new_pf)

def runImport():
    myRun = lagos.EnzoRun(os.path.basename(os.getcwd()), classType=lagos.EnzoParameterFile)
    a = []
    dirs = glob.glob("*.dir")
    dirs.sort()
    for i in dirs:
        bn = os.path.join(i,i[:i.rfind(".")])
        a.append(bn)
    myRun.addOutputByFilename(a, lagos.EnzoParameterFile)
    submitRun(myRun)

def selectRun():
    d=dialog.Dialog()
    d.setBackgroundTitle('FIDO GOES FOR WALKIES!!!!')
    runs_ids, runs_mds = getRuns()
    runs = zip(runs_mds.keys(), runs_mds.keys())
    retCode,md = d.menu("Which run do you want?", choices=runs)
    if retCode == 1:
        return None
    myRun = fetchRun(md, classType = lagos.EnzoParameterFile)
    return myRun
    
def runFetch():
    myRun = selectRun()
    if myRun == None:
        sys.exit()
    runFido({'myRun':myRun})

def runFido(my_dict = {}):
    fs = FidoShell(my_dict)
    fs.interact(FidoBanner)

class FidoShell(InteractiveConsole):
    def __init__(self, my_dict = {}):
        import yt.fido, yt.lagos, yt.raven, yt.enki
        my_dict.update({'fido':yt.fido, \
                        'lagos':yt.lagos, \
                        'raven':yt.raven, \
                        'enki':yt.enki, \
                        'fidoHelp':fidoHelp, \
                        'selectRun':selectRun})
        InteractiveConsole.__init__(self, locals=my_dict)

def run():
    exeName = os.path.basename(sys.argv[0])
    if exeName == "fbranch":
        runBranch()
    elif exeName == "fdigup":
        runDigup()
    elif exeName == "frevert":
        runRevert()
    elif exeName == "fimport":
        runImport()
    elif exeName == "ffetch":
        runFetch()
    elif exeName == "fido":
        print "Hey!"
        runFido()
