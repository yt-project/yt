# This is an example of a manner of dealing with an OutputCollection
#
# We first prompt the user for which of the fido'd runs he wants, and then we
# run the particle exporter on that.  (This is tuned to star particles.)

import yt.fido as fido
import yt.lagos as lagos

oc = fido.selectCollection()

#
# Here is a quick little way to *add* an output collection to your centralized
# db.  Note that you only have to do this once!  (Also, you should comment out
# the above line if you want to use this method.
#
#Giles = fido.Watcher()
#for i in range(42):
    #Giles.dealWithOutput("DataDump%04i.dir/DataDump%04i" % (i,i))
#oc = Giles.oc
#

oc.runFunction(lagos.EnzoHierarchy.exportParticlesPB, prop='h',
   kwargs = {'filename':'%(fn)s_stars.pb',
             'filter':2,
             'fields':["creation_time"],
             'scale':100.0})
