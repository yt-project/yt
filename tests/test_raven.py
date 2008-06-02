"""
Test that we can make plots
"""

import unittest, glob, os.path, os, sys, StringIO

print "Reporting from %s" % (os.getcwd())
sys.path = ['.'] + sys.path

from yt.config import ytcfg
ytcfg["yt","LogLevel"] = '50'
ytcfg["yt","logFile"] = "False"
ytcfg["yt","suppressStreamLogging"] = "True"
ytcfg["lagos","serialize"] = "False"

import numpy as na
import yt.lagos
import yt.raven

# The dataset used is located at:
# http://yt.spacepope.org/DD0018.zip
fn = "DD0010/moving7_0010"
fn = os.path.join(os.path.dirname(__file__),fn)

class RavenTestingBase:
    def setUp(self):
        self.OutputFile = yt.lagos.EnzoStaticOutput(fn)
        self.hierarchy = self.OutputFile.hierarchy
        self.pc = yt.raven.PlotCollection(self.OutputFile)
        self.v, self.c = self.hierarchy.find_max("Density")
        gp = os.path.join(os.path.dirname(fn),"*.yt")
        ytFiles = glob.glob(gp)
        for i in ytFiles:
            os.unlink(i)
        self.localSetup()

    def tearDown(self):
        if hasattr(self,'data'): del self.data
        if hasattr(self,'region'): del self.region
        if hasattr(self,'ind_to_get'): del self.ind_to_get
        if hasattr(self,'pc'): del self.pc
        del self.OutputFile, self.hierarchy
        
    def DoSave(self):
        fns=self.pc.save("test")
        for fn in fns:
            os.unlink(fn)

    def _testSlice(self):
        self.pc.add_slice("Density",0)
        self.pc.plots[-1].switch_z("CellMass")
        self.pc.set_width(0.5,'1')
        self.pc.set_zlim(1,1000)
        self.pc.set_cmap("hot")
        self.DoSave()

    def _testProjection(self):
        self.pc.add_projection("Temperature",1,weight_field="Density")
        self.pc.set_width(0.5,'1')
        self.pc.set_zlim(1,1000)
        self.pc.set_cmap("hot")
        self.DoSave()

    def _testThreePhaseSphere(self):
        print "Testing ThreePhase"
        self.pc.add_phase_sphere(1.0,'1',["Density","Temperature","Density"],center=self.c)
        self.DoSave()

    def _testCallbacksOnSlices(self):
        # We test a couple things here
        # Add callbacks, then remove one, then do the plot-saving
        # Both types of callback should be called here
        for ax in range(3):
            self.pc.add_slice("Density", 0)
            x,y = yt.raven.axis_labels[ax]
            v1 = "%s-velocity" % (x)
            v2 = "%s-velocity" % (y)
            qi = self.pc.plots[-1].add_callback(yt.raven.QuiverCallback(v1,v2,ax,32))
            ti = self.pc.plots[-1].add_callback(yt.raven.ContourCallback("Temperature",
                                               ncont=3, factor=10))
            gi = self.pc.plots[-1].add_callback(yt.raven.ContourCallback("Gas_Energy",
                                               ncont=3, factor=10))
            self.pc.plots[-1].remove_callback(gi)
        self.DoSave()

class PlotTestingBase(RavenTestingBase):
    def test_set_xlim(self):
        self.pc.set_xlim(0.25,0.75)
        self.DoSave()

    def test_set_ylim(self):
        self.pc.set_ylim(0.25,0.75)
        self.DoSave()

    def test_autoscale(self):
        # verify autoscale changed
        self.pc.autoscale()
        self.DoSave()

    def test_set_zlim(self):
        self.pc.set_zlim(0.5, 1.0)
        self.DoSave()

    def test_set_lim(self):
        self.pc.set_lim((0.25,0.75,0.25,0.75))
        self.DoSave()

    def test_set_width(self):
        self.pc.set_width(0.25,'1')
        self.DoSave()

    def test_set_cmap(self):
        self.pc.set_cmap("kamae")
        self.DoSave()
        self.pc.set_cmap("jet")
        self.DoSave()

    def test_switch_field(self):
        for field in ["Temperature","x-velocity"]:
            self.pc.switch_field(field)
            # Check if the logging is set correctly
            self.DoSave()

    def test_clear_plots(self):
        self.pc.clear_plots()
        self.assertTrue(len(self.pc.plots) == 0)

    def test_set_label(self):
        for p in self.pc.plots: p.set_label(r"$\rm{Hi}$")
        self.DoSave()
        for p in self.pc.plots: p.set_label("Hi!")
        self.DoSave()

    def test_set_logfield(self):
        for p in self.pc.plots: p.set_log_field(False)
        self.DoSave()
        for p in self.pc.plots: p.set_log_field(False)
        self.DoSave()

    def test_save(self):
        self.DoSave()

class TestSlices(PlotTestingBase, unittest.TestCase):
    def localSetup(self):
        self.pc.add_slice("Density",0)
        self.pc.add_slice("Density",1)
        self.pc.add_slice("Density",2)

class TestSphere(PlotTestingBase, unittest.TestCase):
    def localSetup(self):
        self.pc.add_phase_sphere(1.0,'1',
                ["Density","TotalEnergy","y-velocity"])

class TestPhaseObject(PlotTestingBase, unittest.TestCase):
    def localSetup(self):
        obj = self.hierarchy.region([0.5]*3, [0.0]*3, [1.0]*3)
        self.pc.add_phase_object(obj, ["Density","TotalEnergy","y-velocity"])

class TestProjection(PlotTestingBase, unittest.TestCase):
    def localSetup(self):
        self.pc.add_projection("Density", 0)
        self.pc.add_projection("Temperature", 1)
        self.pc.add_projection("x-velocity", 2, weight_field="Density")

class TestMixProjectionSlice(PlotTestingBase, unittest.TestCase):
    def localSetup(self):
        self.pc.add_projection("Density",0)
        self.pc.add_slice("Density",0)

class TestCuttingPlane(PlotTestingBase, unittest.TestCase):
    def localSetup(self):
        self.pc.add_cutting_plane("Density", [0.1,0.2,0.3])

if __name__ == "__main__":
    unittest.main()
