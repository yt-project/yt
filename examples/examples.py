"""
Examples for many common operations.  First we define the functions to use,
and then we call them if appropriate.
"""

import yt.lagos as lagos
import yt.fido as fido

fn = "LOCATION_OF_HIERARCHY"

def returnStaticOutput(fn=fn):
    """
    This is simply how one instantiates a static output, given a filename.
    
    I've also added a call to printStats, to give some more info.
    """
    o = lagos.EnzoStaticOutput(fn)
    o.hierarchy.printStats()
    return o

def getStarsForPartiview(hierarchy, fn):
    """
    This can be called on any type of particle, and the fields
    can be selected.  This is an example of a standard call.
    """
    hierarchy.exportParticlesPB(
        fn, filter=2, scale=100.0,
        fields=["metallicity_fraction", "particle mass"])

def getDMForPartiview(hierarchy, fn):
    a.exportParticlesPB(fn, scale=100.0, filter=1)
    
def getClusterFiles(basename):
    """
    This just gives an example of how to instantiate cluster outputs,
    and the concatenate them together.  It assumes you're using disk and species
    output, too.
    """
    cl1=lagos.AnalyzeClusterOutput(basename)
    cl2=lagos.AnalyzeClusterOutput(basename + ".Species")
    cl3=lagos.AnalyzeClusterOutput(basename + ".Disk")
    return cl1 + cl2 + cl3

def getWeightedProfile(hierarchy):
    """
    Simple weighted profile support is available.  You pass in the 
    fields you want to a call to EnzoSphere.makeProfile, as well as
    the inner and outer radii in code units.  It should mostly take
    care of automatically deciding on weights, as appropriate, for you.
    """
    v, c = hierarchy.findMax("Density")
    region = lagos.EnzoSphere(
        a, c, 100./a["au"], ["CellMass","Density","Volume"])
    bins,profiles=region.makeProfile(
        ["Density", "NumberDensity"],10,1./a["au"],100./a["au"])
    return bins, profiles

def runWatcher():
    def watcherFunc():
        # We do all the setup here,
        import yt.lagos as lagos, yt.raven as raven
        # then define our function,
        def runFunc(filename):
            a = lagos.EnzoStaticOutput(filename)
            pc = raven.PlotCollection(a)
            pc.addProjection("Density", 0, center=[0.5,0.5,0.5])
            pc.save("temp")
        # then return the function.
        # That way we only have to import lagos and raven if we need to to make
        # plots.  Simply running Fido doesn't require them.
        return runFunc
    Giles = fido.Watcher(functionHandler = watcherFunc)
    Giles.run()
