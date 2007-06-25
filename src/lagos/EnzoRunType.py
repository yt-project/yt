"""
Handling of sets of EnzoHierarchy objects

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

from yt.lagos import *
import ConfigParser, time, os, sys

class EnzoRun:
    """
    A class that is used to hold the information about an entire Enzo
    Simulation.  This includes all of the datadumps, the parameter file,
    and possibly even the initial conditions.
    """
    def __init__(self, metaData, outputs=[], runFilename=None, classType=EnzoHierarchy, timeID=None, getPFs = True):
        """
        We're going to try to avoid setting too many of the parameters here,
        as many will be changed on and off.  However, we can definitely set all
        of the parameters that are fixed from the initial conditions -- things
        like root grid conditions, problem type, and so on.  While these *may*
        change during the run (for instance, changing to SN problem type) they
        won't do so without quite a bit of effort on the part of the user, and
        we thus don't quite have to care about them too much.

        This is primarily a storage container.  If we open it up, and add to
        it all of our datadumps, we can snag lots of fun information from them
        all.

        @param metaData: text describing the run as a whole.
        @type metaData: string
        @keyword outputs: outputs to add to the run
        @type outputs: list of L{EnzoHierarchies<EnzoHierarchy>}
        """
        self.metaData = metaData
        self.classType = classType
        self.outputs = obj.array(outputs)       # Object array of EnzoHierarchies
        self.timesteps = na.zeros(shape=self.outputs.shape, dtype=nT.Float64) # Timesteps
        self.runFilename = runFilename
        self.gotRuns = False
        if runFilename != None and getPFs:
            self.Import(runFilename)
        if not timeID:
            timeID = int(time.time())
        self.timeID = timeID

    def promoteType(self, outputID):
        if hasattr(self.outputs[outputID], "grids"):
            return
        pp = self.outputs[outputID]
        pf_name = os.path.join(pp.fullpath, pp.basename)
        self.outputs[outputID] = EnzoHierarchy(pf_name)
        del pp

    def demoteType(self, outputID):
        if not hasattr(self.outputs[outputID], "grids"):
            return
        pp = self.outputs[outputID]
        pf_name = os.path.join(pp.fullpath, pp.basename)
        self.outputs[outputID] = EnzoParameterFile(pf_name)
        del pp

    def sortOutputs(self):
        """
        Sorts outputs, solely by the time at which they were created.

        Note that this may not be what we want -- we may actually want the time
        in the simulation at which they were dumped.
        """
        order = na.argsort(self.timesteps)
        self.outputs = self.outputs[order]
        self.timesteps = self.timesteps[order]

    def addOutput(self, hierarchy):
        """
        Add an output.  Also, sorts.

        @param hierarchy: either a single hierarchy or a list of hierarchies
        @type hierarchy: L{EnzoHierarchy<EnzoHierarchy>}
        """
        if not isinstance(hierarchy, types.ListType):
            hierarchy = [hierarchy]
        # Our arrays are both one-d, so we'll just extend them
        t = []
        for h in hierarchy:
            t.append(h["InitialTime"])
        self.outputs = obj.array(self.outputs.tolist() + hierarchy)
        self.timesteps=na.array(self.timesteps.tolist() + t,dtype=nT.Float64)
        self.sortOutputs()
        self.gotRuns = True

    def addOutputByFilename(self, filename, hdf_version=4):
        """
        Feed it a list of parameter files, andi t will add 'em all

        @param filename: either a single filename or a list of filenames
        @keyword hdf_version: version of hdf to use
        @type hdf_version: int
        """
        if not isinstance(filename, types.ListType):
            filename = [filename]
        k = []
        for fn in filename:
            mylog.debug("Adding %s to EnzoRun '%s'", fn, self.metaData)
            k.append(self.classType(fn, hdf_version=hdf_version))
            k[-1].run = self
        self.addOutput(k)

    def getCommandLine(self):
        return "./enzo_red_i9_r16"

    def getTime(self):
        return time.ctime(float(self.timeID))


    def runFunction(self, func, args, fmt_string = None):
        """
        Here we can call runFunction, feeding it a function and some arguments.

        The function will be called with every EnzoHierarchy as the first
        argument, and the *args as the rest of the arguments.

        @param func: function handler to be called for every EnzoHierarchy
                     note that it will be called as such:
                     func(EnzoHierarchy, [fmt_string % index], *args)
        @type func: function
        @param args: a list of arguments to pass in
        @type args: list
        @param fmt_string: this is an optional argument, which will be %'d
                           against the index of the parameter file.  Note that
                           it is only %'d against the *index*, so even if the
                           thing is called DataDump1042, if it's the first one, it'll get 0.
        @type fmt_string: string
        """
        if not isinstance(args, types.ListType):
            args = [args]
        for i in range(self.outputs.shape[0]):
            pid = os.fork()
            if pid:
                newpid, exit = os.wait()
                mylog.info("Exit status %s from PID %s", exit, newpid)
            else:
                mylog.debug("Forked process reporting for duty")
                self.promoteType(i)
                a = [self.outputs[i]]
                if fmt_string != None:
                    s = fmt_string % i
                    a += [s]
                a += args
                mylog.debug("Calling %s on %s", func.func_name, a[0].parameterFilename)
                func(*a)
                mylog.debug("Done calling %s, now dying", func.func_name)
                sys.exit()

    def getBefore(self, time):
        return na.where(self.timesteps <= time)[0]

    def getAfter(self, time):
        return na.where(self.timesteps > time)[0]

    def removeOutput(self, key):
        pp = self[key]
        id = self.index(key)
        # Okay, I know this sucks, but, it seems to be the only way, since
        # concatenate doesn't want to work with ObjectArrays
        newOutputs = obj.array(self.outputs[:id].tolist() + self.outputs[id+1:].tolist())
        newTimesteps = obj.array(self.timesteps[:id].tolist() + self.timesteps[id+1:].tolist())
        self.outputs = newOutputs
        self.timesteps = newTimesteps
        del pp, newOutputs, newTimesteps

    def removeOutputFiles(self, outputID):
        """
        This function handles things a bit  more simply than
        removeOutputFilesByGrid, in that it tosses everything over to fido.
        This *deletes* files.  Note that this can be DANGEROUS.  Use with care.

        @param outputID: the ID or basename of the parameter file to toast
        @type outputID: int or string
        """
        import yt.fido
        mylog.info("Deleting %s" % (self[outputID].fullpath))
        #yt.fido.deleteFiles(self.outputs[outputID].fullpath, \
                #os.path.basename(self.outputs[outputID].parameterFilename))
        yt.fido.deleteOutput(self[outputID].fullpath)

    def moveOutput(self, outputID, dest):
        """
        This moves a data directory tree to another location.  This assumes the
        data is already buried.

        @param outputID: the output we want to move
        @type outputID: integer
        @param dest: the basename of the new location
        @type dest: string
        """
        import yt.fido
        pf = self[outputID]
        mylog.info("Moving %s to %s", os.path.basename(pf.fullpath), \
                    dest)
        pf.fullpath = yt.fido.moveOutput(pf.fullpath, os.path.abspath(dest))

    def moveOutputFiles(self, outputID, dest):
        import yt.fido
        pf = self[outputID]
        mylog.info("Moving %s to %s", os.path.basename(pf.fullpath), \
                    dest)
        yt.fido.moveFiles(os.path.abspath(dest), \
            os.path.basename(pf.parameterFilename), \
            extraFiles = [], wd=pf.directory)

    def removeOutputFilesByGrid(self, outputID):
        """
        This function deletes all the files associated with a given parameter
        file, and then the instance itself.  I am of two minds here; do we want
        to execute the costly operation of removing each file associated with
        the hierarchy instance, or do we want to just kill 'em all off via a
        system command?
        
        @param outputID: the ID of the parameter file to toast
        @type outputID: int
        @deprecated: Use removeOutputFiles
        """
        run = self.outputs[outputID]
        for grid in run.grids:
            mylog.debug("Deleting %s", grid.filename)
            os.unlink(grid.filename)
        mylog.debug("Deleting parameter file %s", run.parameterFilename)
        os.unlink(run.parameterFilename)
        mylog.debug("Deleting hierarchy %s", run.hierarchyFilename)
        os.unlink(run.hierarchyFilename)
        mylog.debug("Deleting boundary %s", run.boundaryFilename)
        os.unlink(run.boundaryFilename)
        mylog.debug("Deleting boundary.hdf  %s", run.boundaryFilename+".hdf")
        os.unlink(run.boundaryFilename+".hdf")
        mylog.debug("Deleted everything, I think")

    def Export(self, nn):
        """
        This function returns a ConfigParser object that describes the enough
        that it can be resurrected at a later date.

        @return: ConfigParser
        """
        cfg = ConfigParser.ConfigParser()
        if os.path.exists(nn):
            cfg.read(nn)
        if not cfg.has_section("RunInformation"):
            cfg.add_section("RunInformation")
        cfg.set("RunInformation","MetaData",self.metaData)
        cfg.set("RunInformation","User", os.getenv("USER"))
        cfg.set("RunInformation","ExportTime", time.ctime())
        # Hm, I guess that's all we need for the metadata...
        # Now the question, really, is going to be how do we handle the data
        # listings?
        self.UpdateOutputs(cfg)
        return cfg

    def Import(self, filename):
        """
        Here we read in an EnzoRun file and import the parameter files
        @param filename: the run filename
        @type filename: string
        """
        cc = ConfigParser.ConfigParser()
        cc.read(filename)
        self.metaData = cc.get("RunInformation","MetaData")
        new_outputs = []
        for sec in cc.sections():
            if sec.startswith("Output"):
                # Okay, so it's an output, now what?
                d = cc.get(sec,"Directory")
                f = cc.get(sec,"Basename")
                new_outputs.append(os.path.join(d,f))
        self.gotRuns = True
        self.addOutputByFilename(new_outputs)

    def UpdateOutputs(self, cfg):
        """
        This function accepts a filename, and will update the ConfigParser
        object corresponding to that filename, ensuring that all outputs I{currently}
        associated with the run instance are in the file.  (It does this by
        I{deleting} the [Outputs] section, by the way, and re-initializing it.)
        @param cfg: The filename or ConfigParser object to update
        @type cfg: string or instantiated ConfigParser
        """
        if isinstance(cfg, types.StringType):
            # We open the file, and we will write to it when done
            cc = ConfigParser.ConfigParser()
            cc.read(cfg)
            toWrite = True
        else:
            # Assume it is a ConfigParser object, and let exception get thrown
            # otherwise.
            cc = cfg
            toWrite = False
        for sec in cc.sections():
            if sec.startswith("Output"):
                cc.remove_section(sec)
        num = self.outputs.shape[0]
        otherParams=ytcfg.get("Fido","OtherParamsToStore").split(",")
        for i in range(num):
            n = "Output%04i" % (i)
            cc.add_section(n)
            cc.set(n,"CurrentTimeIdentifier",self.outputs[i]["CurrentTimeIdentifier"])
            cc.set(n,"Directory",self.outputs[i].fullpath)
            cc.set(n,"Basename",self.outputs[i].basename)
            for op in otherParams:
                if self.outputs[i].has_key(op):
                    cc.set(n,op,self.outputs[i][op])
        if toWrite:
            cc.write(cfg)
        # We return *nothing*, because we either modify the ConfigParser in
        # place or we update and write out

    def __getitem__(self, key):
        """
        Based on the type of input, we either return based on index or
        basename.
        """
        if isinstance(key, types.StringType):
            index = self.index(key)
        elif isinstance(key, types.IntType):
            index = key
        return self.outputs[index]

    def index(self, key):
        t = os.path.basename(key)
        # Find out the index
        index = None
        for i in range(self.outputs.shape[0]):
            if self.outputs[i].basename == t:
                index = i
                break
        if index == None:
            raise KeyError
        return index

    def getOutputs(self):
        if not self.gotRuns:
            self.Import(self.runFilename)
        return self.outputs.tolist()

    def __xattrs__(self, mode="default"):
        return("metaData", "getTime()")

    def __repr__(self):
        return self.metaData

    def __len__(self):
        if not self.gotRuns:
            self.Import(self.runFilename)
        return self.outputs.shape[0]

    def __iter__(self):
        if not self.gotRuns:
            self.Import(self.runFilename)
        for i in range(len(self)):
            yield self.outputs[i]

    def keys(self):
        """
        Returns a list of parameterfile basenames
        """
        keys = []
        for i in range(self.outputs.shape[0]):
            keys.append(self.outputs[i].basename)
        return keys

    def exportAmira(self, basename, fields):
        a5 = basename + "-%(field)s.a5"
        sorted_times = []
        sorted_timesteps = []
        for field in fields:
            a5b=tables.openFile(a5 % {'field':field},"w")
            a5b.close()
        if (not iterable(fields)) or (isinstance(fields, types.StringType)):
            fields = [fields]
        for i in range(len(self)):
            mylog.info("Exporting %s to Amira format", self[i])
            bn = basename + "-%06i" % (i) + "-%(field)s.h5"
            numDone = 0
            self.promoteType(i)
            self[i].exportAmira(bn, fields, a5, i)
            rootDelta = na.array([self[i].grids[0].dx, self[i].grids[0].dy, self[i].grids[0].dz], dtype=nT.Float64)
            sorted_timesteps.append(i)
            sorted_times.append(self[i]["InitialTime"])
            self.demoteType(i)
            numDone += 1
        for field in fields:
            a5b=tables.openFile(a5 % {'field':field},"a")
            a5b.createGroup("/","globalMetaData")
            node=a5b.getNode("/","globalMetaData")
            node._f_setAttr("datatype",0)
            node._f_setAttr("fieldtype",1)
            node._f_setAttr("staggering",1)
            node._f_setAttr("maxTime",self.timesteps[-1])
            node._f_setAttr("maxTimeStep",len(self)-1)
            node._f_setAttr("minTime",self.timesteps[0])
            node._f_setAttr("minTimeStep",0)
            node._f_setAttr("numTimeSteps",len(self))
            node._f_setAttr("rootDelta",rootDelta)
            a5b.createArray("/","sorted_timesteps", na.array(sorted_timesteps, dtype=nT.Int32))
            a5b.createArray("/","sorted_times", na.array(sorted_times, dtype=nT.Float64))
            a5b.close()
