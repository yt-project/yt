"""
Test that we can make plots
"""

import unittest, glob, os.path, os

import yt.raven
from yt.config import ytcfg
print "Setting loglevel"
ytcfg["yt","LogLevel"] = '50'

import yt.raven

# The dataset used is located at:
# http://yt.spacepope.org/DD0018.zip
fn = "DD0018/moving7_0018"
fn = os.path.join(os.path.dirname(__file__),fn)

class TestRaven(unittest.TestCase):
    def setUp(self):
        self.OutputFile = yt.lagos.EnzoStaticOutput(fn)
        self.hierarchy = self.OutputFile.hierarchy
        self.pc = yt.raven.PlotCollection(self.OutputFile)
        v,self.c = self.hierarchy.findMax("Density")

    def tearDown(self):
        #self.hierarchy.data_file.close()
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

    def testCallbacksOnSlices(self):
        # We test a couple things here
        # Add callbacks, then remove one, then do the plot-saving
        # Both types of callback should be called here
        for ax in range(3):
            self.pc.addSlice("Density", ax)
            x,y = yt.raven.axis_labels[ax]
            v1 = "%s-velocity" % (x)
            v2 = "%s-velocity" % (y)
            qi = self.pc.plots[-1].addCallback(yt.raven.be.quiverCallback(v1,v2,ax,32))
            ti = self.pc.plots[-1].addCallback(yt.raven.be.contourCallback("Temperature",
                                               ax, ncont=3, factor=10))
            gi = self.pc.plots[-1].addCallback(yt.raven.be.contourCallback("Gas_Energy",
                                               ax, ncont=3, factor=10))
            self.pc.plots[-1].removeCallback(gi)
        self.DoSave()

    def testParticlesOnSlice(self):
        # We test a couple things here
        # Add callbacks, then remove one, then do the plot-saving
        # Both types of callback should be called here
        for ax in range(3):
            self.pc.addSlice("Density", ax)
            gi = self.pc.plots[-1].addCallback(yt.raven.be.particleCallback(
                                               ax, 0.01))
        self.DoSave()

if __name__ == "__main__":
    unittest.main()
