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
import yt.shell
from yt.fido import *

def runDigup():
    guess = os.path.basename(os.getcwd())
    run = fetchRun(guess)
    if len(sys.argv) <= 1:
        new_pf = selectPF(run, prompt="Which run should I dig up?")
        if new_pf == None:
            sys.exit()
    else:
        new_pf = guessPF(sys.argv[1])
    digUp(run, new_pf)
    sys.exit()

def runRevert():
    guess = os.path.basename(os.getcwd())
    run = fetchRun(guess)
    if len(sys.argv) <= 1:
        new_pf = selectPF(run, prompt="Which is the last run to keep?")
        if new_pf == None:
            sys.exit()
    else:
        new_pf = guessPF(sys.argv[1])
    revert(run, new_pf)

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

def selectPF(run, prompt = "Which parameter file do you want?"):
    """
    Choose a parameterfile from a menu, with menu options populated by the database
    @param run: the run to select from
    @type run: L{EnzoRun<EnzoRun>}
    @param prompt: the prompt to give the user
    @type prompt: string
    """
    opts=run.keys()
    d=dialog.Dialog()
    d.setBackgroundTitle('FIDO GOES FOR WALKIES!!!!')
    pfs = zip(map(str,range(1,len(run)+1)),run.keys())
    retCode,pfI = d.menu(prompt, choices=pfs)
    if retCode == 1:
        return None
    myPF = run[int(pfI)-1]
    return myPF

def selectRun():
    """
    Choose a run from a menu, with menu options populated by the database
    """
    d=dialog.Dialog()
    d.setBackgroundTitle('FIDO GOES FOR WALKIES!!!!')
    runs_ids, runs_mds = getRuns()
    runs = zip(map(str,range(1,len(runs_mds)+1)),runs_mds.keys())
    print runs
    retCode,mdI = d.menu("Which run do you want?", choices=runs)
    if retCode == 1:
        return None
    myRun = fetchRun(runs[int(mdI)-1][1], classType = lagos.EnzoParameterFile)
    return myRun
    
def runFetch():
    myRun = selectRun()
    if myRun == None:
        sys.exit()
    newBanner = yt.shell.ytBannerDefault + \
                "\nmyRun = %s" % (myRun)
    yt.shell.runyt({'myRun':myRun}, ytBanner=newBanner)

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
