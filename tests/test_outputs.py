"""
Test that we can get outputs, and interact with them in some primitive ways.
"""

import unittest
import yt.lagos

fn = "DD0018/moving7_0018"

class baseTest(unittest.TestCase):
    def setUp(self):
        self.OutputFile = yt.lagos.EnzoStaticOutput(fn)
        self.hierarchy = self.OutputFile.hierarchy

    def tearDown(self):
        del self.OutputFile

class TestEnzoOutputs(baseTest):
    def testGetHierarchy(self):
        self.assert_(self.OutputFile.hierarchy != None)
        
    def testGetUnits(self):
        self.assert_(self.OutputFile["cm"] != 1.0)

class TestEnzoHierarchy(baseTest):
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

class TestRegionData(baseTest):
    def setUp(self):
        baseTest.setUp(self)
        self.r=self.hierarchy.region(
                          [0.5,0.5,0.5],[0.0, 0.0, 0.0],
                          [1.0, 1.0, 1.0],
                          ["CellMass"])
        self.s=self.hierarchy.sphere(
                          [0.5,0.5,0.5],2.0,["CellMass"])

    def testMassSum(self):
        ms = self.s["CellMass"].sum()
        mr = self.r["CellMass"].sum()
        self.assertEqual(ms,mr)

if __name__ == "__main__":
    unittest.main()
