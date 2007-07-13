"""
Test that we can get outputs, and interact with them in some primitive ways.
"""

import unittest
import yt.lagos

fn = "/Users/matthewturk/Research/data/mornkr/galaxy0398.dir/galaxy0398"

class TestEnzoOutputs(unittest.TestCase):
    def setUp(self):
        self.OutputFile = yt.lagos.EnzoStaticOutput(fn)
    
    def testGetHierarchy(self):
        self.assert_(self.OutputFile.hierarchy != None)
        
    def testGetUnits(self):
        self.assert_(self.OutputFile["cm"] != 1.0)
        
    def tearDown(self):
        del self.OutputFile
        
class TestEnzoHierarchy(unittest.TestCase):
    def setUp(self):
        self.OutputFile = yt.lagos.EnzoStaticOutput(fn)
        self.hierarchy = self.OutputFile.hierarchy
        
    def testGetSmallestDx(self):
        self.assertAlmostEqual(self.hierarchy.getSmallestDx(),
                               7.62939453125e-06, 5)
        
    def testGetNumberOfGrids(self):
        self.assertEqual(self.hierarchy.numGrids, len(self.hierarchy.grids))
        self.assertEqual(self.hierarchy.numGrids, 949)
        
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
        
if __name__ == "__main__":
    unittest.main()