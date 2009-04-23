#!/usr/bin/python

# Stephen Skory
# sskory@physics.ucsd.edu
# created August 2006
# modified June 2008
# this python script builds a script for creating halo merger trees using hop inside yt
# see http://yt.enzotools.org/ for more info on yt and this script.

# the output of GraphViz at the very will have boxes connected by arrows.
# The arrow into a box is how much of the previous halo's DM & stars is
# contributing to that halo.
# Similarly the arrow leaving a box is how much of the DM & stars is leaving
# that halo for the next one.
# Inside the box, the top % is how much of the halo is accreted from the 
# 'ether'; particles that in the previous time step were assigned to no group.
# Similarly, the bottom % is how much of the halo leaves and is assigned to no
# group in the next time step.
# The middle of the box shows how many star/DM particles are in the halo,
# and the (x,y,z) position of the halo's center.

# to run this, set all the things below to your liking,
# and then 'python fastBuildMerge.py > mergeScript.py
# 'chmod u+x mergeScript.py'
# './mergeScript.py'

# the name of the python script
name = "merger_yt"
# the GraphViz file
outfile = "157-120.dot"
# the directory basename, no trailing slash needed
# if a data dump is /path/to/enzo/data/DD0243/data0243, put /path/to/enzo/data/DD below
dirbasename = "/path/to/enzo/data/DD"
# usually DD, maybe RD, or data
filebasename = "data"
# the maximum number of groups analyzed at one time step
hardmaxgroup = 2500
# the hop density threshold used for grouping
hopthreshold = 80.0
# dm only, True or False
dmonly = "True"

# The range of data dumps you wish to operate upon. It's a good idea to manually see how far back
# haloes exist before setting the end value.
# start > end
end = 120
start = 157

# A list of the haloes in the final data dump for which you want to find ancestors. Use the IDs provided
# by hop, so you'll probably have to run hop on the final data dump before deciding which haloes to
# put here, unless you know it's the largest halo, which is always ID=0. It can be a singleton.
indices = [0,1]

# change below to reflect where python is
# it is important to use python >= 2.5 as previous versions have a memory bug
# which affects this script, which is important if you're running large datasets.
print "#!/usr/bin/env /Library/Frameworks/Python.framework/Versions/Current/bin/python"

# ---------- you shouldn't have to change anything below. ---------
# -----------------------------------------------------------------
# -----------------------------------------------------------------

print "import %s" % name
print "import yt.raven as raven"
print "import yt.lagos as lagos"
print "from yt.lagos import mylog"
print "from yt.mods import *"

print """
def instantiateSnapshots(dirbasename,filebasename,start,end):
    redshift = {}
    snapshots = {}
    for i in range(start-end):
        j = start - i
        file = "%s%04d/%s%04d" % (dirbasename,j,filebasename,j)
        snapshot = lagos.EnzoStaticOutput(file)
        redshift[j] = snapshot["CosmologyCurrentRedshift"]
        # put the snapshot instance into the dict snapshots
        snapshots[j] = snapshot
    
    string = "Read in %d snapshots." % len(redshift)
    mylog.info(string)
    return(snapshots,redshift)
"""


print "positions = {}"
print "dirbasename = \"%s\"" % dirbasename
print "filebasename = \"%s\"" % filebasename
print "hardmaxgroup = %d" % hardmaxgroup
print "indices = " + str(indices)
print "start = %d" % start
print "end = %d" % end
print "dmonly = \"%s\"" % dmonly

# first let's read in the redshifts, and instantiate the snapshots
print "(snapshots,redshift) = instantiateSnapshots(dirbasename,filebasename,start,end-1)"

# open the dotfile
print "fopen = open('%s', 'w')" % outfile

# write the preambles for the dotfile
print "%s.writeTop(fopen)" % name
print "%s.writeNode(fopen)" % name

print "\n"

# loop over all the snapshots
for i in range(start-end):
    snapshot = start - i
    # snapshots should be in time-reverse order, so this loop will count *down*, so the next snapshot
    # to load has a lower index
    nextsnapshot = snapshot-1
    # if we're at the bottom, we're done!
    if (nextsnapshot < end):
        break
    
    # we need to run hop on the snapshot
    print "hop_results = HaloFinder.HOPHaloFinder(snapshots[%d], threshold=%f, dm_only=dmonly)" % \
        (nextsnapshot,hopthreshold)
    
    # get the group info for the next snapshot
    # get group positions
    print "positions = %s.convertPositions(hop_results,%d,positions,hardmaxgroup)" % (name,nextsnapshot)
    # get particle IDs
    print "g%04d = %s.convertGroups(hop_results,%d,hardmaxgroup)" % (nextsnapshot,name,nextsnapshot)
    # sort the particles in each group
    print "for g in g%04d: g.particles.sort()" % nextsnapshot
    # delete the sphere which contains all the particles, to clear up memory
    print "del sphere"
    # delete hop_results now that we're done with them, too
    print "del hop_results\n"
    
    # in the beginning we need to read in two groups before we build the links, so do everything above
    # again
    if (snapshot==start):
        print "hop_results = HaloFinder.HOPHaloFinder(snapshots[%d], threshold=%f, dm_only=dmonly)" % \
            (snapshot,hopthreshold)
        print "positions = %s.convertPositionsSelected(hop_results,%d,positions,indices)" % (name,snapshot)
        print "g%04d = %s.convertGroupsSelected(hop_results,%d,indices)" % (snapshot,name,snapshot)
        print "for g in g%04d: g.particles.sort()\n" % snapshot
        print "del sphere"
        print "del hop_results"
    
    
    # for just the last group (in time), a singleton, usually.
    if (snapshot==start):
        print "for g in g%04d: g.flag = 1\n" % start
    
    # build links which requires two neighboring groups at once.
    print "(links%04d%04d,g%04d,g%04d) = %s.buildLinks(g%04d,g%04d,positions)" % \
    (nextsnapshot,snapshot,nextsnapshot,snapshot,name,nextsnapshot,snapshot)
    
    # write out levels which keeps the boxes on one line rather than put how graphviz wants aesthetically
    
    print "%s.writeLevels(fopen, g%04d,redshift[%d])" % (name,nextsnapshot,nextsnapshot)
    if (snapshot==start):
        print "%s.writeLevels(fopen, g%04d,redshift[%d])" % (name,snapshot,snapshot)
    
    # print out label stuff
    
    print "%s.writeOpen(fopen)" % name
    
    # for the first (last in time) halo, write out the box labels using style 2
    if (snapshot==start):
        print "%s.writeLabels(fopen, g%04d, positions, 2)" % (name,snapshot)
    # middle time steps get type 1
    if (snapshot!=start):
        print "%s.writeLabels(fopen, g%04d, positions, 1)" % (name,snapshot)
    # the top, first in time get type 0
    if (nextsnapshot==end):
        print "%s.writeLabels(fopen, g%04d, positions, 0)" % (name,nextsnapshot)
    
    print "%s.writeClose(fopen)" % name
    
    print "%s.writeOpen(fopen)" % name
    print "%s.writeLinks(fopen,links%04d%04d, g%04d, g%04d)" % (name,nextsnapshot,snapshot,nextsnapshot,snapshot)
    print "%s.writeClose(fopen)" % name

    print "for g in g%04d: del g.particles\n\n" % snapshot

# close the dotfile
print "%s.writeClose(fopen)" % name
print "fopen.close()"
