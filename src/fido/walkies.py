"""
This module handles both the interactive component as well as the command-line
scripts.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

#from yt.fido import *  # Not sure we want this
import dialog, sys, optparse, shutil
import glob, os, os.path

import yt.shell
from yt.fido import *

from IPython.Extensions.ipipe import *

from yt.config import ytcfg

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

def runFido():
    pp = optparse.OptionParser(description="Archives data and creates plots.", \
                               version="Your Mom")
    pp.add_option("-c","--config", action="store", dest="cfgFile", \
                  default="raven_config.xml", type="string")
    pp.add_option("-p","--path", action="store", dest="path", \
                  default="", help="Where to move files when archiving")
    pp.add_option("-d","--deliverator", action="store_true", dest="deliverator", \
                  default=False, \
                  help="Should we submit to the Deliverator?\n(This requires that httpPrefix be defined in [raven] in your yt config file!)")
    opts, args = pp.parse_args()
    imagePath = ytcfg.get("raven", "imagePath", raw=True)
    imageSkel = ytcfg.get("raven", "imageSkel", raw=True)
    import yt.fido as fido
    import yt.raven as raven
    def fidoPlots(h):
        md = h["MetaDataString"]
        if md.strip().rstrip() == "(null)":
            fido.error("MetaDataString undefined in parameter file; either use it, or write your own durn function!")
            raise KeyError
        mdd = {'md':md}
        if opts.deliverator:
            httpPrefix = ytcfg.get("raven", "httpPrefix")
            RunID, Response = deliverator_upload.SubmitRun(md, ytcfg.get("deliverator","user"))
            plot=raven.EnzoHippo(h, submitToDeliverator=RunID, httpPrefix=httpPrefix % mdd)
        else:
            plot=raven.EnzoHippo(h)
        imageDir = imagePath % mdd
        if not os.path.isdir(imageDir):
            os.makedirs(imageDir)
        prefix = os.path.join(imageDir, imageSkel)
        raven.MakePlots(plot, opts.cfgFile, prefix)
        return
    fido.watchDir(funcHandler=fidoPlots, newPrefix=opts.path) 

def runGetNewEnzo():
    pp = optparse.OptionParser(description="Recompiles Enzo based on configuration file.", \
                               version="Your Mom")
    pp.add_option("-j","--numProc", action="store", dest="nproc", \
                                    default=1, type="int")
    pp.add_option("-p","--path", action="store", dest="path", \
                help="Where to copy resultant executable (defaults to no copying)")
    pp.add_option("-c","--clean", action="store_true", dest="clean", \
                                  default=False)
    opts, args = pp.parse_args()
    import yt.enki as enki
    newExe = enki.compileNew(nProc=opts.nproc, Clean=opts.clean)
    if opts.path != None:
        mylog.info("Copying new enzo exe to %s", opts.path)
        shutil.copy2(newExe, opts.path)

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
    # Wait, why is it awful, March-Matt?  Sincerely, May-Matt

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
    import yt.lagos as lagos
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
    import yt.lagos as lagos
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
    elif exeName == "fido":
        runFido()

def runBrowse():
    """
    Trying out the new ipipe functionality of IPython
    """
    rids, rmds = getRuns()
    runs = []
    for md in rmds.keys():
        runs.append(fetchRun(md, getPFs=False))
    bw = ibrowse()
    runs | bw
    myRun = bw.display()
    return myRun

def outputBrowse(run):
    ll = list(run.outputs)
    bw = ibrowse()
    ll | bw
    return bw.display()

# Here is our spec sheet:
#
# fbranch:
#   Takes output to branch at, as well as new metadata and new directory.
#   Require these options!
# fdigup:
#   Which one to dig up
# frevert:
#   Which one to revert to
# ffetch:
#   We can disregard this one for now
# fbury:
#   Does this one even work?
# frecomp:
#   -j : number of processors
#   -c : clean first?
#   -p : where should we toss the executable?
#   -t : target to build?
