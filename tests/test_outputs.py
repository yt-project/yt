"""
Test that we can get outputs, and interact with them in some primitive ways.
"""

import unittest, glob, os.path, os

from yt import ytcfg
ytcfg["yt","LogLevel"] = '20'

import yt.lagos

# The dataset used is located at:
# http://yt.spacepope.org/DD0018.zip
fn = "DD0018/moving7_0018"

class TestDataContainers(unittest.TestCase):
    def setUp(self):
        self.OutputFile = yt.lagos.EnzoStaticOutput(fn)
        self.hierarchy = self.OutputFile.hierarchy
        self.v, self.c = self.hierarchy.findMax("Density")
        gp = os.path.join(os.path.dirname(fn),"*.yt")
        ytFiles = glob.glob(gp)
        for i in ytFiles: 
            print "Removing %s" % (i)
            os.unlink(i)

    def tearDown(self):
        self.hierarchy.dataFile.close()
        del self.OutputFile, self.hierarchy

    def testGetHierarchy(self):
        self.assert_(self.OutputFile.hierarchy != None)
        
    def testGetUnits(self):
        self.assert_(self.OutputFile["cm"] != 1.0)

    def testGetSmallestDx(self):
        self.assertAlmostEqual(self.hierarchy.getSmallestDx(),
                               0.0009765625, 5)
        
    def testGetNumberOfGrids(self):
        self.assertEqual(self.hierarchy.numGrids, len(self.hierarchy.grids))
        self.assertEqual(self.hierarchy.numGrids, 211)
        
    def testChildrenOfRootGrid(self):
        for child in self.hierarchy.grids[0].Children:
            self.assert_(child.Parent.id == self.hierarchy.grids[0].id)
            
    def testGetSelectLevels(self):
        grids = self.hierarchy.grids
        for level in range(self.hierarchy.maxLevel+1):
            for grid in grids[self.hierarchy.selectLevel(level)]:
                self.assert_(grid.Level == level)
                
    def testPrintStats(self):
        try:
            self.hierarchy.printStats()
        except:
            self.fail()
        self.assert_(True)

    def testRegions(self):
        # This tests a couple things
        r=self.hierarchy.region(
                     [0.5,0.5,0.5],[0.0, 0.0, 0.0],
                     [1.0, 1.0, 1.0],
                     ["Density","Temperature"])
            # Testing multiple fields fed in
        s=self.hierarchy.sphere(
                     [0.5,0.5,0.5],2.0,
                     ["CellMass","Temperature"])
            # Testing multiple fields fed in
        ms = s["CellMass"].sum() # Testing adding new field transparently
        mr = r["CellMass"].sum() # Testing adding new field transparently
        self.assertEqual(ms,mr)  # Asserting equality between the two

    def testProjection(self):
        p = self.hierarchy.proj(0,"Density") # Unweighted
        p = self.hierarchy.proj(1,"Temperature","Density") # Weighted
        p = self.hierarchy.proj(2,"Entropy") # Derived field

    def testSlice(self):
        s = self.hierarchy.slice(0,0.5,"Density") # Non-derived
        s = self.hierarchy.slice(1,0.5,"RadialVelocity") # Derived
        s = self.hierarchy.slice(2,0.5,["RadialVelocity","Temperature"])
                         # Multiple, derived

if __name__ == "__main__":
    unittest.main()
