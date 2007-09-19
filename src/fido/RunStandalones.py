"""
@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
  Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

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

import optparse
import os, os.path, sys, re, time
from yt.fido import *
import exceptions

def selectOutput(collection = None):
    if not collection: collection = selectCollection()
    for i, o in enumerate(collection):
        print "%s\t %s" % (i+1,o)
    print
    loki = int(raw_input("Which output do you want?"))-1
    return collection[loki-1]

def selectCollection():
    cs = GrabCollections()
    for i, o in enumerate(cs):
        print "%s\t %s" % (i+1,o)
    print
    loki = int(raw_input("Which collection do you want?"))-1
    return cs[loki-1]

import optparse

class FidoAction:
    """
    This class defines the skeleton of every action's options.
    """
    def __init__(self):
        self.parser = optparse.OptionParser(
                                description=self.description,
                                version="SVN")

    def SetupParser(self):
        pass

    def GuessOC(self, createNew = False):
        tryTitle=os.path.basename(os.getcwd())
        gc = GrabCollections()
        for c in gc:
            if c.title == tryTitle:
                self.oc = c
                return
        self.oc = None
        if createNew:
            Import().PerformAction()
            self.GuessOC(createNew=False)

    def ParseArgs(self):
        self.SetupParser()
        self.opts, self.args = self.parser.parse_args()

    def PerformAction(self):
        pass

    def CleanUp(self):
        pass

class Bury(FidoAction):
    description = "Bury an output in a subdirectory"
    def __init__(self):
        FidoAction.__init__(self)
        self.ParseArgs()
        self.GuessOC()
        if not self.oc:
            # Should we import here?
            raise KeyError, "Needs to be imported before you can bury."

    def PerformAction(self):
        for bn in self.args:
            newName = buryOutput(bn)
            self.oc.addOutput(newName)
        self.oc.writeOut()

class DigUp(FidoAction):
    description = "Un-bury an output"
    def __init__(self):
        FidoAction.__init__(self)
        self.ParseArgs()
        self.GuessOC()
        if not self.oc:
            # Should we import here?
            raise KeyError, "Needs to be imported before you can bury."

    def SetupParser(self):
        self.parser.add_option("-d", "--dest",
                               action='store', type='string',
                               dest='newLocation', default=os.getcwd())

    def PerformAction(self):
        for bn in self.args:
            b = bn
            if b.endswith('.dir/'):
                b = os.path.join(b, b[:-5])
            elif b.endswith('.dir'):
                b = os.path.join(b, b[:-4])
            b = os.path.basename(b)
            print b
            toGlob = self.oc[b][0]
            digupOutput(toGlob, newLocation=self.opts.newLocation)
            del self.oc[b]
        self.oc.writeOut()

class Branch(FidoAction):
    description = "Un-bury an output"
    def __init__(self):
        FidoAction.__init__(self)
        self.ParseArgs()
        if self.opts.newLocation == None:
            print "You must supply a location for the branching."
            print "See --help ."
            print
            raise KeyError
        self.GuessOC()
        if not self.oc:
            # Should we import here?
            raise KeyError, "Needs to be imported before you can branch."

    def SetupParser(self):
        self.parser.add_option("-d", "--newdir",
                               action='store', type='string',
                               dest='newLocation', default=None)
        self.parser.add_option("-m", "--metadatastring",
                               action='store', type='string',
                               dest='md', default=None)

    def PerformAction(self):
        bn = self.args[-1] # only do the last one
        # First we copy, then we modify the MDS.  Note that we don't need to
        # create a new OC.
        b = bn
        if b.endswith('/'):
            b = os.path.basename(b)
        if b.endswith('.dir'):
            b = os.path.join(b, b[:-4])
        #print "Copying %s to %s" % (b, self.opts.newLocation)
        newName=copyOutput(b, self.opts.newLocation)
        print newName
        #Okay, now that we've copied it...
        #First we update the CurrentTimeIdentifier string
        expr=re.compile("CurrentTimeIdentifier\s*=\s*[0-9]*$",re.M)
        pfContents = open(newName).read()
        newId = time.time()
        newPf=expr.sub("CurrentTimeIdentifier = %i" % (newId), pfContents)
        expr=re.compile("MetaDataString\s*=\s*.*$",re.M)
        if self.opts.md == None: self.opts.md = os.path.basename(
                                                 os.path.dirname(
                                                  os.path.abspath(newName)))
        newPf=expr.sub("MetaDataString      = %s" % (self.opts.md),newPf)
        open(newName,"w").write(newPf)

class Import(FidoAction):
    description = "Import an existing set of buried outputs"
    def __init__(self):
        FidoAction.__init__(self)
        self.ParseArgs()
        self.GuessOC()
        self.title = os.path.basename(os.getcwd())
        if self.oc != None: self.title = self.oc.title
        else: self.oc=OutputCollection(self.title)
        open(NEW_OUTPUT_CREATED,"w").close()

    def PerformAction(self):
        for i in glob.glob("*.dir/*.hierarchy"):
            fn = i[:-10]
            self.oc.addOutput(os.path.abspath(fn))
        Giles = Watcher(title=self.title, oc=self.oc)
        Giles.run(True)
        self.oc.writeOut()

class FidoStandalone(FidoAction):
    description = "Run Fido, all by itself."
    def __init__(self):
        FidoAction.__init__(self)
        self.ParseArgs()
        self.GuessOC()
        self.title = os.path.basename(os.getcwd())
        if self.oc != None: self.title = self.oc.title

    def SetupParser(self):
        self.parser.add_option("-s", "--submit",
                               action='store_true', 
                               dest='submit', default=False)
        self.parser.add_option("-c", "--cfg",
                               action='store', type='string',
                               dest='cfgFile', default=None)
        self.parser.add_option("-p", "--newpath",
                               action='store', type='string',
                               dest='newPrefix', default=None)

    def PerformAction(self):
        if self.opts.cfgFile != None: func=self.MakePlots
        else: func = None
        Giles = Watcher(title=self.title, oc=self.oc, functionHandler=func,
                        newPrefix=self.opts.newPrefix)
        Giles.run()
        self.oc.writeOut()

    def MakePlots(self):
        import yt.lagos as lagos
        import yt.raven as raven
        import yt.raven.deliveration as deliverator
        imagePath = ytcfg.get("raven", "imagePath", raw=True)
        imageSkel = ytcfg.get("raven", "imageSkel", raw=True)
        def fidoPlots(fn):
            pf = lagos.EnzoStaticOutput(fn)
            md = pf["MetaDataString"].strip().rstrip()
            if md == "(null)":
                fido.error("MetaDataString undefined in parameter file; either use it, or write your own durn function!")
                raise KeyError
            if self.opts.cfgFile != None:
                RunID = -1
                if self.opts.submit:
                    httpPrefix = ytcfg.get("raven", "httpPrefix", raw=True)
                    RunID, Response = deliverator.SubmitRun(md, ytcfg.get("Deliverator","user"))
                imageDir = imagePath % pf
                if not os.path.isdir(imageDir):
                    os.makedirs(imageDir)
                prefix = os.path.join(imageDir, imageSkel)
                raven.MakePlots(pf, self.opts.cfgFile, prefix, RunID)
            return
        return fidoPlots

def runAction():
    pg = os.path.basename(sys.argv[0])
    acts = {'fbury':Bury, 'fdigup':DigUp,
            'fbranch':Branch, 'fimport':Import,
            'fido':FidoStandalone}
    #try:
    acts[pg]().PerformAction()
    #except:
        #print "Error: terminating."
