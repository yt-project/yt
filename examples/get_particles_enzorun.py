import yt.lagos as lagos

myRun = lagos.EnzoRun("my_awesome_run")

a = []
for i in range(1,5):
    a.append("galaxy%04i.dir/galaxy%04i" % (i,i))

myRun.addOutputByFilename(a)

# We're putting the arguments in the order they are expected by the function

myRun.runFunction(lagos.EnzoHierarchy.exportParticlesPB, [1,None,100.0], fmt_string="dm%04i.pb")
myRun.runFunction(lagos.EnzoHierarchy.exportParticlesPB, [2,["metallicity_fraction","particle mass"],100.0], fmt_string="stars%04i.pb")
