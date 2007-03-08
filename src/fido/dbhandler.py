"""
Not yet functional.  Better strategy needed first...

So the goal here will be to have a set of classes that we can use in a
flat-file database.  In fact, everything will be stored in flat-files in a
directory.  They probably won't be pickled.

Characteristics of a given run:
===============================

    0. Unique ID
    1. MetaData identifier string (including date instantiated)
    2. User
    3. Date instantiated (separate sorting field)
    4. All associated runs

Characteristics of a given parameter file:
==========================================

    0. Unique creation-time
    1. MetaData identifier string
    2. Join: I{Backwards} run-association (in filename?)
    3. Contents of parameter file
    4. User (for double-checking)
    5. (Unique) directory location of data

We only want to associate a parameter file with the run used to I{create} it.
This enables us to duplicate a run file, associate a new creation time, and
thusly I{branch} the run into a separate run.  Typically this would be done if
the underlying code were changed, for instance, or if comparison studies were
run between late-time evolution.  We have no real interest in creating forward
association of parameter files with runs.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}

@todo: Implement!
@todo: Create mappers for each object type.
@todo: Pickle objects?  Deep-pickling, or surface?  Or just do flat-file
storage?
"""

RUN_PREFIX = "runF"

from yt.fido import *
import yt.lagos as lagos
import glob, re, types

my_rundir = ytcfg.get("Fido","rundir")

def getRuns():
    """
    This returns two dictionaries, one with IDs as the keys, the other with
    metaData strings as keys
    """
    possible = glob.glob(os.path.join(my_rundir,RUN_PREFIX+"_*"))
    runs_ids = {}
    runs_md = {}
    for f in possible:
        run, id, md = f.split("_",2)
        if runs_ids.has_key(id):
            mylog.warning("Duplicate run ID found: %s (overwriting old)", id)
        runs_ids[id] = (md, f)
        if runs_md.has_key(md):
            mylog.warning("Duplicate run metaData found: %s (overwriting old)", md)
        runs_md[md] = (id, f)
        mylog.debug("Found: (%s) %s", id, md)
    return runs_ids, runs_md

def submitRun(run, additional = []):
    """
    This accepts an EnzoRun instance and puts it into the rundir, in a filename
    of the form %(id)_%(md) where id is generated as the current time.

    @param run: the EnzoRun instance
    @type run: L{EnzoRun<EnzoRun>}
    """
    #id = int(time.time())
    id = run.timeID
    md = run.metaData
    new_filename = "_".join([RUN_PREFIX,str(id),md])
    nn = os.path.join(my_rundir,new_filename)
    cfg=run.Export(nn)
    for set in additional:
        cfg.set("RunInformation",set[0],set[1])
    cfg.write(open(os.path.join(my_rundir,new_filename),"w"))
    mylog.info("Added '%s' to %s", md, new_filename)
    return

def fetchRun(md, classType = lagos.EnzoParameterFile):
    """
    This takes a metaData string and then returns the EnzoRun instance
    associated with that metaData string

    @param md: metaData
    @type md: string
    @return: L{EnzoRun<EnzoRun>}
    """
    ids, mds =getRuns()
    if not mds.has_key(md):
        raise KeyError, "MetaData string (%s) not found!" % (md)
    run = lagos.EnzoRun(md, runFilename=mds[md][1], classType=classType, timeID=mds[md][0])
    return run

def submitParameterFile(run, pf):
    """
    This inserts a parameterfile into the run_dir 'db'
    @param run: run instance
    @type run: L{EnzoRun<EnzoRun>}
    @param pf: parameterfile
    @type pf: L{EnzoHierarchy<EnzoHierarchy>}
    @deprecated: I don't think we want to ever, ever use this.
    """
    ids, mds = getRuns()
    if not mds.has_key(run.metaData):
        submitRun(run)
        ids, mds = getRuns()
    rid = mds[run.metaData]
    #pid = int(pf.parameters["CurrentTimeIdentifier"])
    #md = "_".join([pf.filename,pf["MetaDataString"]])
    #new_filename = "P%s_" + "_".join([id,pf.filename,pf["MetaDataString"]])
    new_filename = pf.GetUniqueFilename()
    mylog.info("Storing new parameter file %s", new_filename)
    cfg=run.Export()
    cfg.write(os.path.join(my_rundir,new_filename))
    return

def fetchParameterFile(id=None, filename=None):
    pass

def branch(run, pf, new_md, new_dir):
    """
    This will duplicate a parameter file and all of its consitutent data.
    Additionally, the parameter file gets a new ID.
    @param run: the run to branch
    @type run: L{EnzoRun<EnzoRun>}
    @param pf: the parameter file to branch at
    @type pf: integer or L{EnzoParameterFile<EnzoParameterFile>}
    @param new_md: the new metadata identifier
    @type new_md: string
    @param new_dir: the location to I{copy} the data to
    @type new_dir: string
    @return: the new L{EnzoRun<EnzoRun>}
    @note: The EnzoRun this contains does B{not} contain EnzoHierarchy objects,
    but instead EnzoParameterFile objects.  This should make things much less
    intense when dealing with large hierarchies.
    """
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    if isinstance(pf, types.StringType):
        pf = lagos.EnzoParameterFile(pf)
    # First we instantiate a new EnzoRun
    new_run = lagos.EnzoRun(metaData=new_md, classType=lagos.EnzoParameterFile)
    # Now we copy the old (buried) data to a new location
    mylog.info("Copyin'")
    ll = copyOutput(pf.directory, new_dir)
    # Now in the new parameter file, we update CurrentTimeIdentifier
    expr=re.compile("CurrentTimeIdentifier\s*=\s*[0-9]*$",re.M)
    pf_contents = open(os.path.join(pf.fullpath, pf.basename)).read()
    # Replace CurrentTimeIdentifier with current time
    new_id = time.time()
    new_pf=expr.sub("CurrentTimeIdentifier = %i" % (new_id), pf_contents)
    # Replace old MetaData with new
    expr=re.compile("MetaDataString\s*=\s*.*$",re.M)
    new_pf=expr.sub("MetaDataString      = %s" % (new_md),new_pf)
    new_filename = os.path.join(ll, pf.basename)
    mylog.info("Branched id %s to id %s with %s" % (pf["CurrentTimeIdentifier"],new_id, new_md))
    open(new_filename,"w").write(new_pf)
    # Now, add the things to the new run, and go
    new_run.addOutputByFilename(new_filename)
    # Now we cann the parameter file submitter, which also takes care of
    # submitting the run
    submitRun(new_run, additional=[("BranchedFrom","%s_%s" % (run.timeID, run.metaData))])
    # And now, it's branched!  Hooray.
    return new_run

def revert(md, maxTime):
    """
    When fed a maximum time, this will delete all the parameter files that come
    after that time.  Note that this is a "Greater-than" operation, not a
    "Greater-than-or-equal-to" operation.

    @param md: the run we're reverting
    @type md: string
    @param maxTime: the maximum time
    @type maxTime: float
    @note: Does this necessarily work with FP precision?  Do we need to do a
           +1e-30 or something to it?
    """
    thisRun = fetchRun(md)
    try:
        # Maybe they passed in a float?
        tt = float(maxTime)
    except ValueError:
        # Otherwise, we'll have to figure it out from the time of the output file
        for i in thisRun.outputs:
            if i.basename == os.path.basename(maxTime):
                # This is the one
                tt = i["InitialTime"]
    bad_outputs = thisRun.getAfter(tt)
    good_outputs = thisRun.getBefore(tt)
    for oI in bad_outputs:
        mylog.info( "Removing OutputFiles for %s (%s)", oI, thisRun.outputs[oI].fullpath)
        thisRun.removeOutputFiles(oI)
    newRun = lagos.EnzoRun(md, outputs=thisRun.outputs[(good_outputs,)], classType=lagos.EnzoParameterFile, timeID=thisRun.timeID)
    # Note that we're not reinstantiating the config parser object, since the
    # whole point is to I{delete} that and rewrite it!
    submitRun(newRun)

def guessRun(pf=None, wd="."):
    """
    Based on a parameter file and the working directory, it attempts to guess
    the appropriate run

    @param pf: the parameter file
    @type pf: string or EnzoParameterFile
    @param wd: the working directory
    @type wd: string
    """
    if pf:
        try:
            md = pf["MetaDataString"]
        except:
            pf = lagos.EnzoParameterFile(pf)

def digUp(md, pf):
    """
    When fed a parameter file, this will 'dig it up' into the current
    directory.  Note that this *dissociates* it from the EnzoRun.

    @param md: the run we're digging from
    @type md: string
    @param pf: the pf to dig up
    @type pf: float
    @note: Does this necessarily work with FP precision?  Do we need to do a
           +1e-30 or something to it?
    """
    thisRun = fetchRun(md)
    # Otherwise, we'll have to figure it out from the time of the output
    # file
    for i in range(thisRun.outputs.shape[0]):
        print pf, thisRun.outputs[i].basename, os.path.basename(pf)
        if thisRun.outputs[i].basename == os.path.basename(pf):
            goodI = i
            break
    thisRun.moveOutputFiles(goodI, dest=os.getcwd())
    thisRun.removeOutput(goodI)
    submitRun(thisRun)

def guessPF(pf):
    """
    Simple way of guessing a parameter file's name.
    """
    try_dir = os.path.join(pf+".dir",pf)
    if pf[-1] == "/":
        pf = pf[:-1]
    if pf[-3:] == "dir":
        # So we guess that we want the parameterfile with the same basename as
        # the dir
        new_pf = os.path.join(pf, pf[:pf.rfind(".dir")])
    elif not os.path.exists(pf) and os.path.exists(try_dir):
        new_pf = try_dir
    else:
        new_pf = pf
    print "new_pf:", new_pf
    return new_pf
