"""
Test that we can make plots
"""

import unittest, glob, os.path, os

from yt import ytcfg
ytcfg["yt","LogLevel"] = '20'

import yt.raven

# The dataset used is located at:
# http://yt.spacepope.org/DD0018.zip
fn = "/Users/rpwagner/Data/yt_test/DD0018/moving7_0018"

class TestRaven(unittest.TestCase):
    def setUp(self):
        self.OutputFile = yt.lagos.EnzoStaticOutput(fn)
        self.hierarchy = self.OutputFile.hierarchy
        self.pc = yt.raven.PlotCollection(self.OutputFile)
        v,self.c = self.hierarchy.findMax("Density")

    def tearDown(self):
        self.hierarchy.dataFile.close()
        del self.OutputFile, self.hierarchy, self.pc

    def DoSave(self):
        fns=self.pc.save("test")
        for fn in fns:
            print "Unlinking",fn

    def testSlice(self):
        self.pc.addSlice("Density",0)
        self.pc.plots[-1].switch_z("Temperature")
        self.pc.set_width(0.5,'1')
        self.pc.set_zlim(1,1000)
        self.pc.set_cmap("hot")
        self.DoSave()

    def testProjection(self):
        self.pc.addProjection("Temperature",1,weightField="Density")
        self.pc.set_width(0.5,'1')
        self.pc.set_zlim(1,1000)
        self.pc.set_cmap("hot")
        self.DoSave()

    def testThreePhaseSphere(self):
        print "Testing ThreePhase"
        self.pc.addThreePhaseSphere(1.0,'1',["Density","Temperature","Density"],center=self.c)
        self.DoSave()

    def testTwoPhaseSphere(self):
        print "Testing TwoPhase"
        self.pc.addTwoPhaseSphere(1.0,'1',["Temperature","Density"],center=self.c)
        self.DoSave()
        
if __name__ == "__main__":
    unittest.main()
