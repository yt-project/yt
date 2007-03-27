"""
This module handles both the interactive component as well as the command-line
scripts.

@todo: implement runBranch()

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

#from yt.fido import *  # Not sure we want this
import dialog, sys, getopt, shutil
import glob, os.path

import yt.lagos as lagos
import yt.shell
from yt.fido import *

def runDigup():
    guess = os.path.basename(os.getcwd())
    run = fetchRun(guess)
    if len(sys.argv) <= 1:
        new_pf = selectPF(run, prompt="Which output should I dig up?")
        if new_pf == None:
            sys.exit()
    else:
        new_pf = guessPF(sys.argv[1])
    digUp(run, new_pf)
    sys.exit()

def runGetNewEnzo(toPut = None):
    import yt.enki as enki
    opts,args = getopt.getopt(sys.argv[1:], 'j:cp:')
    j = 1
    clean = False
    for o,v in opts:
        if o == "-j":
            j = v
            mylog.info("Using %s procs", j)
        elif o == "-c":
            clean = True
            mylog.info("Cleaning")
        elif o == "-p":
            toPut = v
    target = None
    if len(args) > 0:
        target = " ".join(args)
    newExe = enki.compileNew(nProc=j, Clean=clean)
    if toPut:
        mylog.info("Copying new enzo exe to %s", toPut)
        shutil.copy2(newExe, toPut)

def runBranch():
    """
    @todo: Fix up the new directory setup.  Right now, is ugly.  UGGGGGLY.
    """
    guess = os.path.basename(os.getcwd())
    run = fetchRun(guess)
    if len(sys.argv) <= 1:
        new_pf = selectPF(run, prompt="Which is the output to branch at?")
        if new_pf == None:
            sys.exit()
    else:
        new_pf = guessPF(sys.argv[1])
    new_md=getNewMD(initial = "%s_1" % (guess))
    new_dir=getNewDirectory(new_md)
    # This is awful, but I need to do it
    branch(run, new_pf, new_md, new_dir)

def runRevert():
    guess = os.path.basename(os.getcwd())
    run = fetchRun(guess)
    if len(sys.argv) <= 1:
        new_pf = selectPF(run, prompt="Which is the last output to keep?")
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

def getNewDirectory(new_md = ""):
    """
    Ask the user for a new directory, guessing at first the current directory,
    gone up one step, and then with the new md appended
    """
    d=dialog.Dialog()
    d.setBackgroundTitle('FIDO GOES FOR WALKIES!!!!')
    guess=os.path.join(os.path.dirname(os.getcwd()),new_md)
    retCode,new_dir = d.inputbox("What is the new directory?", init=guess)
    if retCode == 1:
        return None
    return new_dir
    

def getNewMD(prompt = "What is the new metadata string?", initial = ""):
    """
    Ask the user for a new metadata string, prompting with the old one
    """
    d=dialog.Dialog()
    d.setBackgroundTitle('FIDO GOES FOR WALKIES!!!!')
    retCode,new_md = d.inputbox(prompt, init=initial)
    if retCode == 1:
        return None
    return new_md
    

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

def runBury():
    """
    This will just look in the current directory, bury everything, and then
    drop out
    """
    runOnce()

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
    elif exeName == "fbury":
        runBury()
    elif exeName == "frecomp":
        runGetNewEnzo()
