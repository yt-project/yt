# Merger.py. Written by Stephen Skory & Rick Wagner, August 2006.

# figures the 'family tree,' and outputs in GraphViz format, the galaxy
# merger history for a series of files that list the particles & the group that
# they are assigned to.

# this file only contains functions that are used by the output of fastBuildMerge.py.

# June 2008, converted/upgraded for use in Matthew Turk's enzo analyzer yt by Stephen Skory

# define the Group class for the list of groups, which contains
# the group id, a list of the particles in that group, the numerical
# order of the group, the percentage of the groups particles that go to the 
# ether (no group), and the percentage of the group that comes from no group.
class Group:


    def __init__(self, id, particles, orderIndex, toEther, fromEther, flag):
        self.id = id
        self.particles = particles
        self.orderIndex = orderIndex
        # these are fields to track how much of each group comes and goes to the
        # ether
        self.toEther = toEther
        self.fromEther = fromEther
        # the flag says whether or not this group is linked to the final,
        # interesting group(s).
        self.flag = flag

# the heavy lifting. This takes in a parentGroup and sees how many of its particles
# go to the childGroup.
# DEPRECIATED, replaced by pushListCount. I've left it here because it's still functional and
# a lot simpler to understand than pushListCount. If you want to use it, uncomment out the line
# where it's used in buildLinks, and comment out the pushListCount line. This is a whole bunch slower
# than pushListCount, so I don't know why you would want to do that.

    def isParent(self, childGroup):

        n = 0
        for particleNumber in self.particles:
            if childGroup.particles.count(particleNumber):
                #return True
                n+=1

        return n

# this function takes two lists which have both been sorted in increasing order. This is a way to
# take the O(n^2) matching and make it roughly O(2n). Since most of the time most particles go from
# one group to the next, a step down one list will get the same next particle in the other list, so
# most of the time we can move down both lists simultaneously. When one entry is bigger than the other,
# we use findNextMatch to see how far to move the second list until it matches the first, or 
# surpasses the first (in which case we then move the first list until it catches up with the second, 
# ad naseum).
def pushListCount(list1, list2):


    count = 0
    i = 0
    j = 0
    
    while (i < len(list1)) and (j < len(list2)):
        #print 'list1['+str(i)+']='+str(list1[i])+' list2['+str(j)+']='+str(list2[j])
        if list1[i] == list2[j]:
            i += 1
            j += 1
            count += 1
            #print 'Matched first try'
        else:
            if (i <= len(list1) -1) and (j <= len(list2) -1):
                value = min(list1[i], list2[j])
                                #print str(value)
                if value == list1[i]:
                    #print 'First list smaller'
                    nextindex = findNextMatch(i+1, list2[j], list1)
                    if nextindex > -1:
                        if nextindex == 1.5: # special case
                            j += 1
                            i += 1
                        else:
                            i = nextindex + 1
                            j += 1
                            count += 1
                    if nextindex < -1:
                        i = -1 * nextindex
                        j += 1
                    if (j == len(list2)) or (i == len(list1)):
                        break
                    if nextindex == -1:
                        #print "Got minus one on second list, nextindex= " + str(nextindex)
                        break

                if value == list2[j]:
                    #print 'Second list smaller'
                    nextindex = findNextMatch(j+1, list1[i], list2)
                    if nextindex > -1:
                        if nextindex == 1.5: # special case
                            j += 1
                            i +=1
                        else:
                            j = nextindex +1
                            i += 1
                            count += 1
                    if nextindex < -1:
                        j = -1 * nextindex
                        i += 1
                    if nextindex == 1:
                        j = nextindex + 1
                        i += 1
                    if (j == len(list2)) or (i == len(list1)):
                        break
                    if nextindex == -1:
                        #print "Got minus one on first list, nextindex= " + str(nextindex)
                        break
        #print 'count: ' + str(count)

    return count

def findNextMatch(startIndex, value, thislist):
    #print "looking for " + str(value) + ', startindex= ' + str(startIndex)
    outval = -1
    for i in range(startIndex, len(thislist)):
        #print 'thislist[' + str(i) +'] = ' +str(thislist[i])
        if thislist[i] == value:
            outval = i
            break
        if thislist[i] > value:
            outval = -1 * i
            if outval == -1: # this in case i=1, which isn't a problem unless outval=-1 *is* the value we want
                         # this happens if the second item in each list matches
                outval = 1.5
            break
    #print 'found match for ' + str(value) + ' at index: ' + str(outval)

    return outval
    


# add yt.hop particles into the Group structure
def convertGroups(snapshot,orderIndex,hardmaxgroup):
    from yt.lagos import mylog
    # ytHopResults is the stuff from yt containing the results of hop.
    # orderIndex is the index for this snapshot
    # hardmaxgroup is some reasonable limit for the number of galaxies looked at for each
    # time step.
    
    groups = []
    particles = []
    # set up the percent fields for later...
    toEther = 100
    fromEther = 100
    pcount = 0 # counts how many particles are groups for that particular time step
    flag = 0 # zero means not assoc with final group, one means yes.
    
    # loop over the groups
    for i,group in enumerate(snapshot):
        # don't go past hardmaxgroup
        if (i >= hardmaxgroup):
            break
        # add it to our groups as a Group
        groups.append(Group(i,group["particle_index"],orderIndex,toEther,fromEther,flag))
        pcount += len(group["particle_index"])
    
    string = "There are " + str(pcount) + " particles in orderIndex " + str(orderIndex) + "."
    mylog.info(string)
    return(groups)


# from the yt.hop output, add to the positions dict the central position and group index for this
# snapshot
def convertPositions(snapshot,orderIndex,positions,hardmaxgroup):
    from yt.lagos import mylog
    # loop over the hop groups
    for i,group in enumerate(snapshot):
        # if we're trying to read in too many groups, we're done so we leave
        if (i >= hardmaxgroup):
            break
        center = group.maximum_density_location()
        positions[orderIndex,i] = (center[0],center[1],center[2])
    
    string = 'There are ' + str(len(snapshot)) + ' groups in orderIndex ' + str(orderIndex)
    mylog.info(string)
    
    return (positions)

# add yt.hop particles into the Group structure, but for selected groups, for the 
# last time step, which is the first inspected.
def convertGroupsSelected(snapshot,orderIndex,indices):
    from yt.lagos import mylog
    # ytHopResults is the stuff from yt containing the results of hop.
    # orderIndex is the index for this snapshot
    # hardmaxgroup is some reasonable limit for the number of galaxies looked at for each
    # time step.
    
    groups = []
    particles = []
    # set up the percent fields for later...
    toEther = 100
    fromEther = 100
    pcount = 0 # counts how many particles are groups for that particular time step
    flag = 0 # zero means not assoc with final group, one means yes.
    
    # loop over the selected groups
    for index in indices:
        # pick out the one group
        group = snapshot[index]
        
        # add it to our groups as a Group
        groups.append(Group(index,group["particle_index"],orderIndex,toEther,fromEther,flag))
        pcount += len(group["particle_index"])
    
    string = "There are " + str(pcount) + " particles in orderIndex " + str(orderIndex) + "."
    mylog.info(string)
    return(groups)

# from the yt.hop output, add to the positions dict the central position and group index for this
# snapshot, but only chosen groups, which is what we want at the final data dump
def convertPositionsSelected(snapshot,orderIndex,positions,indices):
    from yt.lagos import mylog
    
    # loop over the selected groups
    for index in indices:
        group = snapshot[index]
        center = group.maximum_density_location()
        positions[orderIndex,index] = (center[0],center[1],center[2])
    
    string = 'There are ' + str(len(indices)) +' groups in orderIndex ' + str(orderIndex)
    mylog.info(string)
    
    return (positions)



# this does three things:
# 1. Builds the links list, which contains the parentGroup.id and childGroup.id
# of the matched pair, along with the number of particles that are transferred.
# 2. Calculates what percentage of the parentGroup is transferred to the child
# groups. This is then subtracted from 100 to give how much of the group goes to
# the ether.
# 3. Similarly, calculates what percentage of the child group comes from the
# ether.
def buildLinks(parentGroups, childGroups,positions):
    from yt.lagos import mylog
    links = []
    childPercents = {}
    
    for parentGroup in parentGroups:
        parentPercent = 0
        for childGroup in childGroups:
            number = 0 # this will be the number of particles that go from the parent
            # group to the child group.
            # Below we calculate the distance between the groups, taking the periodic
            # boundary conditions into account.
            dist_x = abs(positions[parentGroup.orderIndex,parentGroup.id][0] - positions[childGroup.orderIndex,childGroup.id][0])
            if dist_x > 0.5:
                dist_x = 1 - dist_x
            dist_y = abs(positions[parentGroup.orderIndex,parentGroup.id][1] - positions[childGroup.orderIndex,childGroup.id][1])
            if dist_y > 0.5:
                dist_y = 1 - dist_y
            dist_z = abs(positions[parentGroup.orderIndex,parentGroup.id][2] - positions[childGroup.orderIndex,childGroup.id][2])
            if dist_z > 0.5:
                dist_z = 1 - dist_z
            dist = (dist_x**2 + dist_y**2 + dist_z**2)**0.5
            # If the childgroup isn't flagged, or if the two groups are too far away,
            # the groups don't get checked; number stays equal to zero; and
            # we move to the next parentgroup.
            if (childGroup.flag == 1) and (dist <= 0.1):
                string = 'checking parentGroup (' + str(len(parentGroup.particles)) + '): ' + str(parentGroup.orderIndex) + '.' + str(parentGroup.id) + \
                ' with childGroup (' + str(len(childGroup.particles)) + '): ' + str(childGroup.orderIndex) + '.' + str(childGroup.id)
                mylog.info(string)
                #number = parentGroup.isParent(childGroup)
                number = pushListCount(parentGroup.particles,childGroup.particles)
                #string = 'number: ' + str(number)
                #mylog.info(string)
            # create an entry in the childPercents dictionary if it's not already
            # there
            if childPercents.has_key(childGroup.id) != 1:
                childPercents[childGroup.id] = 0
            
            # If the number of particles from the parentgroup goes to the child group,
            # above some threshold, add this link to the links list; calculate how much
            # of the parent group went to the child group and add that to any existing total
            # of particles leaving the parent group; add the particles going to the child
            # group to it's percentage of particles not coming from the ether; and finally
            # flag the parentgroup as part of the lineage for future checking.
            if (number/float(len(parentGroup.particles)))*100 > 50:
                links.append([parentGroup.id,childGroup.id,number])
                string = 'parentGroup: ' + str(parentGroup.orderIndex) + '.' + str(parentGroup.id) + \
                ' | childGroup: ' + str(childGroup.orderIndex) + '.' + str(childGroup.id) + \
                ' | %% transferred: %2.2f%%' % (number/float(len(parentGroup.particles)) * 100)
                mylog.info(string)
                parentPercent += (number/float(len(parentGroup.particles))*100)
                childPercents[childGroup.id] += (number/float(len(childGroup.particles))*100)
                parentGroup.flag = 1
        # End loop over childGroups.
        # if any of the particles from the parent group went to a child group,
        # change the parent groups toEther to reflect this.
        if parentPercent > 0:
            parentPercent = 100 - parentPercent
            parentGroup.toEther = parentPercent

    # go over the childgroups and see how much of the particles came from
    # the ether.
    for childGroup in childGroups:
        if childPercents[childGroup.id] > 0:
            childPercents[childGroup.id] = 100 - childPercents[childGroup.id]
            childGroup.fromEther = childPercents[childGroup.id]

    return (links,parentGroups,childGroups)

# The functions below just write out simple blocks that GraphViz needs.
def writeTop(fopen):

    line = 'digraph galaxy {size="40,40";\n'
    fopen.write(line)

def writeClose(fopen):

    line = '};\n'
    fopen.write(line)

def writeOpen(fopen):

    line = '{\n'
    fopen.write(line)

def writeNode(fopen):
    
    line = 'node [color=lightblue, style=bold, shape=record];\n'
    fopen.write(line)

# In order to keep the boxes for the haloes on the same level, GraphViz needs
# to be told which boxes are all on the same level.
# Mar 2007, added lines for redshift levels.
def writeLevels(fopen, Groups, redshift):

    line = '{ rank = same;\n'
    fopen.write(line)

    line = '"' + str(redshift) + '";'
    fopen.write(line)
    
    for Group in Groups:
        if (Group.toEther != 100) or (Group.fromEther != 100):
            line = '"' + str(Group.orderIndex) + '.' + str(Group.id) + '";'
            fopen.write(line)
    
    line = "\n};\n"
    
    fopen.write(line)

    
# Each box (node) has four bits of information:
# 1. For all but the top nodes, the percentage of the particles that come from
# the ether.
# 2. The number of particles in it.
# 3. The center of that group in (0,0,0)->(1,1,1) coordinates
# 4. The all but the bottom notes, the percentage of the particles that go to 
# the ether.

def writeLabels(fopen, Groups, gPositions, switch):

    # switch = 0: very top node, beginning of time
    # switch = 1: intermediate node, many
    # switch = 2: very bottom node, singular (usually), the last in time

    groupLength = len(Groups)
    
    if switch == 0:
        for Group in Groups:
            if (Group.toEther != 100) or (Group.fromEther != 100):
                color = 1 - float(Group.id)/float(groupLength)
                line = '"' + str(Group.orderIndex) + '.' + str(Group.id) + '" [label="{' + str(len(Group.particles))
                # write the position triplet for the group
                line += '\\n(%1.3f,' % float(gPositions[Group.orderIndex,Group.id][0])
                line += '%1.3f,' % float(gPositions[Group.orderIndex,Group.id][1])
                line += '%1.3f)' % float(gPositions[Group.orderIndex,Group.id][2])
                line += '| %2.2f%%}", shape="record",' % Group.toEther
                line += 'color="%1.3f ' % color
                line += '1.0 %1.3f"];\n' % color
                fopen.write(line)
        
    if switch == 1:
        for Group in Groups:
            if (Group.toEther != 100) or (Group.fromEther != 100):
                color = 1 - float(Group.id)/float(groupLength)
                line = '"' + str(Group.orderIndex) + '.' + str(Group.id)
                line += '" [label="{%2.2f%%| {' % Group.fromEther
                line += str(len(Group.particles))
                # write the position triplet for the group
                line += '\\n(%1.3f,' % float(gPositions[Group.orderIndex,Group.id][0])
                line += '%1.3f,' % float(gPositions[Group.orderIndex,Group.id][1])
                line += '%1.3f)' % float(gPositions[Group.orderIndex,Group.id][2])
                line +='}|%2.2f%%}", shape="record"' % Group.toEther
                line += ', color="%1.3f ' % color
                line += '1.0 %1.3f"];\n' % color
                fopen.write(line)
            
    if switch == 2:
        for Group in Groups:
            if (Group.toEther != 100) or (Group.fromEther != 100):
                color = 1 - float(Group.id)/float(groupLength)
                line = '"' + str(Group.orderIndex) + '.' + str(Group.id)
                line += '" [label="{%2.2f%%|' % Group.fromEther
                line += str(len(Group.particles))
                # write the position triplet for the group
                line += '\\n(%1.3f,' % float(gPositions[Group.orderIndex,Group.id][0])
                line += '%1.3f,' % float(gPositions[Group.orderIndex,Group.id][1])
                line += '%1.3f)' % float(gPositions[Group.orderIndex,Group.id][2])
                line += '}", shape="record", color="%1.3f ' % color
                line += '1.0 %1.3f"];\n' % color
                fopen.write(line)

# write out the GraphViz links that connect two groups.
def writeLinks(fopen, links, parentGroups, childGroups):

    for parentGroup in parentGroups:
        for link in links:
            if str(parentGroup.id) == str(link[0]):
                line = '"' + str(parentGroup.orderIndex) + '.' + str(link[0]) + '"->"' + str(parentGroup.orderIndex+1) + '.' + str(link[1])+\
                '"[label="%2.2f%%",color="blue", fontsize=10];\n' % (link[2]/float(len(parentGroup.particles)) * 100)
                fopen.write(line)
