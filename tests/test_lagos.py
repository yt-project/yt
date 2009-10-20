"""
Test that we can get outputs, and interact with them in some primitive ways.
"""

# @TODO: Add unit test for deleting field from FieldInfo

import unittest, glob, os.path, os, sys, StringIO

print "Reporting from %s" % (os.getcwd())
sys.path = ['.'] + sys.path

from yt.config import ytcfg
ytcfg["yt","LogLevel"] = '50'
ytcfg["yt","logFile"] = "False"
ytcfg["yt","suppressStreamLogging"] = "True"
ytcfg["lagos","serialize"] = "False"

import cPickle
import yt.lagos
import yt.lagos.OutputTypes
import numpy as na
from yt.fido import ParameterFileStore

# The dataset used is located at:
# http://yt.spacepope.org/DD0018.zip
fn = "DD0010/moving7_0010"
fn = os.path.join(os.path.dirname(__file__), fn)

class LagosTestingBase:
    def setUp(self):
        self.OutputFile = yt.lagos.EnzoStaticOutput(fn)
        self.hierarchy = self.OutputFile.hierarchy
        self.v, self.c = self.hierarchy.find_max("Density")
        gp = os.path.join(os.path.dirname(fn),"*.yt")
        ytFiles = glob.glob(gp)
        for i in ytFiles:
            #print "Removing %s" % (i)
            os.unlink(i)

    def tearDown(self):
        if hasattr(self,'data'): del self.data
        if hasattr(self,'region'): del self.region
        if hasattr(self,'ind_to_get'): del self.ind_to_get
        del self.OutputFile, self.hierarchy
        
class TestParameterFileStore(unittest.TestCase):
    def setUp(self):
        self.original = (yt.config.ytcfg.get("yt","ParameterFileStore"),
                         yt.config.ytcfg.get("lagos","serialize"))
        ytcfg['yt','ParameterFileStore'] = "testing.csv"
        pfs = ParameterFileStore()
        os.unlink(pfs._get_db_name())
        self.pfs = ParameterFileStore() # __init__ gets called again
        ytcfg['lagos', 'serialize'] = "True"

    def testCacheFile(self):
        pf1 = yt.lagos.EnzoStaticOutput(fn)
        pf2 = self.pfs.get_pf_hash(pf1._hash())
        self.assertTrue(pf1 is pf2)

    def testGrabFile(self):
        pf1 = yt.lagos.EnzoStaticOutput(fn)
        hash = pf1._hash()
        del pf1
        pf2 = self.pfs.get_pf_hash(hash)
        self.assertTrue(hash == pf2._hash())

    def testGetCurrentTimeID(self):
        pf1 = yt.lagos.EnzoStaticOutput(fn)
        hash = pf1._hash()
        ctid = pf1["CurrentTimeIdentifier"]
        del pf1
        pf2 = self.pfs.get_pf_ctid(ctid)
        self.assertTrue(hash == pf2._hash())

    def tearDown(self):
        os.unlink(self.pfs._get_db_name())
        ytcfg['yt', 'ParameterFileStore'] = self.original[0]
        ytcfg['lagos', 'serialize'] = self.original[1]
        self.pfs.__init__()

class TestHierarchy(LagosTestingBase, unittest.TestCase):
    def testGetHierarchy(self):
        self.assert_(self.OutputFile.hierarchy != None)

    def testGetUnits(self):
        self.assert_(self.OutputFile["cm"] != 1.0)

    def testGetSmallestDx(self):
        self.assertAlmostEqual(self.hierarchy.get_smallest_dx(),
                               0.00048828125, 7)

    def testGetNumberOfGrids(self):
        self.assertEqual(self.hierarchy.num_grids, len(self.hierarchy.grids))
        self.assertEqual(self.hierarchy.num_grids, 10)

    def testChildrenOfRootGrid(self):
        for child in self.hierarchy.grids[0].Children:
            self.assert_(child.Parent.id == self.hierarchy.grids[0].id)

    def testGetSelectLevels(self):
        for level in range(self.hierarchy.maxLevel+1):
            for grid in self.hierarchy.select_grids(level):
                self.assert_(grid.Level == level)

    def testPrintStats(self):
        a = sys.stdout
        sys.stdout = StringIO.StringIO()
        try:
            self.hierarchy.print_stats()
            worked = True
        except:
            worked = False
        sys.stdout = a
        self.assert_(worked)

    def testDataTypes(self):
        r=self.hierarchy.region(
                     [0.5,0.5,0.5],[0.0, 0.0, 0.0],
                     [1.0, 1.0, 1.0],
                     ["CellMass","Temperature"])
            # Testing multiple fields fed in
        s=self.hierarchy.sphere(
                     [0.5,0.5,0.5],2.0,
                     ["CellMass","Temperature"])
        ms = s["CellMass"].sum() # Testing adding new field transparently
        mr = r["CellMass"].sum() # Testing adding new field transparently
        self.assertEqual(ms,mr)  # Asserting equality between the two

    def testProjectionCorrectnessMultipleFields(self):
        p = self.hierarchy.proj(0,["Density","Ones"], weight=None) # Unweighted
        self.assertTrue(na.all(p["Ones"] == 1.0))

    def testProjectionMakingMultipleFields(self):
        p = self.hierarchy.proj(0,["Density","Temperature","Ones"],weight_field="Ones") # Unweighted
        # One for each field, pdx, pdy, px, py, and one for the weight
        self.assertEqual(len(p.data.keys()), 8)

    def testProjectionSuccess(self):
        p = self.hierarchy.proj(0,"Density") # Unweighted
        p = self.hierarchy.proj(1,"Temperature","Density") # Weighted
        p = self.hierarchy.proj(2,"Entropy") # Derived field

    def testUnweightedProjectionCorrectness(self):
        # Now we test that we get good answers
        for axis in range(3):
            p = self.hierarchy.proj(axis, "Ones") # Derived field
            self.assertTrue(na.all(p["Ones"] == 1.0))
            # Regardless of weighting, we want ones back

    def testWeightedProjectionCorrectness(self):
        # Now we test that we get good answers
        for axis in range(3):
            # Regardless of weighting, we want ones back
            p = self.hierarchy.proj(axis, "Ones", "Density")
            self.assertTrue(na.all(p["Ones"] == 1.0))

# Now we test each datatype in turn

def _returnFieldFunction(field):
    def field_function(self):
        try:
            self.data[field.name]
            if not field.particle_type and not field.vector_field and \
                self.data[field.name].size > 1:
                self.assertEqual(na.product(self.data["Density"].shape),
                                 na.product(self.data[field.name].shape))
            del self.data[field.name]
        except yt.lagos.ValidationException:
            pass
    return field_function

def _returnProfile1DFunction(field, weight, accumulation, lazy):
    def add_field_function(self):
        self.data.set_field_parameter("center",[.5,.5,.5])
        profile = yt.lagos.BinnedProfile1D(
            self.data, 8, "RadiusCode", 0, 1.0, False, lazy)
        profile.add_fields(field, weight=weight, accumulation=accumulation)
    return add_field_function

def _returnProfile2DFunction(field, weight, accumulation, lazy):
    def add_field_function(self):
        self.data.set_field_parameter("center",[.5,.5,.5])
        cv_min = self.hierarchy.gridDxs.min()**3.0
        cv_max = self.hierarchy.gridDxs.max()**3.0
        profile = yt.lagos.BinnedProfile2D(self.data,
                    8, "RadiusCode", 1e-3, 1.0, True,
                    8, "CellVolumeCode", cv_min, cv_max, True, lazy)
        profile.add_fields(field, weight=weight, accumulation=accumulation)
    return add_field_function

class DataTypeTestingBase:
    def setUp(self):
        LagosTestingBase.setUp(self)

    def testRepr(self):
        self.assertTrue(
            ("%s" % self.data).startswith(self.data.__class__.__name__))

class Data3DBase:
    def testProfileAccumulateMass(self):
        self.data.set_field_parameter("center",[0.5]*3)
        profile = yt.lagos.BinnedProfile1D(self.data, 8, "RadiusCode", 0, 1.0,
                                           False, True)
        profile.add_fields("CellMassMsun", weight=None, accumulation=True)
        v1 = profile["CellMassMsun"].max()
        v2 = self.data["CellMassMsun"].sum()
        v2 = na.abs(1.0 - v2/v1)
        self.assertAlmostEqual(v2, 0.0, 7)

    def testExtractConnectedSetsNoCache(self):
        mi = self.data["Density"].min() * 2.0
        ma = self.data["Density"].max() * 0.99
        cons, contours = self.data.extract_connected_sets(
            "Density", 2, mi, ma)
        self.assertEqual(len(contours), 2) # number of contour levels
        self.assertEqual(len(contours[0]), 2)
        self.assertEqual(len(contours[1]), 1)

    def testExtractConnectedSetsCache(self):
        mi = self.data["Density"].min() * 2.0
        ma = self.data["Density"].max() * 0.99
        cons, contours = self.data.extract_connected_sets(
            "Density", 2, mi, ma, cache=True)
        self.assertEqual(len(contours), 2) # number of contour levels
        self.assertEqual(len(contours[0]), 2)
        self.assertEqual(len(contours[1]), 1)

    def testContoursCache(self):
        cid = yt.lagos.identify_contours(self.data, "Density",
                self.data["Density"].min()*2.00,
                self.data["Density"].max()*1.01)
        self.assertEqual(len(cid), 2)

    def testContoursObtain(self):
        cid = yt.lagos.identify_contours(self.data, "Density",
                self.data["Density"].min()*2.00, self.data["Density"].max()*1.01)
        self.assertEqual(len(cid), 2)

    def testContoursValidityMax(self):
        v1 = self.data["Density"].max()*0.99
        v2 = self.data["Density"].max()*1.01
        cid = yt.lagos.identify_contours(self.data, "Density", v1, v2)
        self.assertTrue(na.all(v1 < self.data["Density"][cid[0]])
                    and na.all(v2 > self.data["Density"][cid[0]]))
        self.assertEqual(len(cid), 1)

    def testContoursValidityMin(self):
        v1 = self.data["Density"].min()*0.99
        v2 = self.data["Density"].min()*1.01
        cid = yt.lagos.identify_contours(self.data, "Density", v1, v2)
        self.assertTrue(na.all(v1 < self.data["Density"][cid[0]])
                    and na.all(v2 > self.data["Density"][cid[0]]))
        self.assertEqual(len(cid), 3)

    def testPickle(self):
        ps = cPickle.dumps(self.data)
        pf, obj = cPickle.loads(ps)
        self.assertEqual(obj["CellMassMsun"].sum(), self.data["CellMassMsun"].sum())

for field_name in yt.lagos.FieldInfo:
    field = yt.lagos.FieldInfo[field_name]
    setattr(DataTypeTestingBase, "test%s" % field.name, _returnFieldFunction(field))

field = "Temperature"
for weight in [None, "CellMassMsun"]:
    for lazy in [True, False]:
        for accumulation in [True, False]:
            func = _returnProfile1DFunction(field, weight, accumulation, lazy)
            name = "test%sProfile1D_w%s_l%s_a%s" % (field,
                                                weight, lazy,
                                                accumulation)
            setattr(Data3DBase, name, func)

for weight in [None, "CellMassMsun"]:
    for lazy in [True, False]:
        for accumulation_x in [True, False]:
            for accumulation_y in [True, False]:
                acc = (accumulation_x, accumulation_y)
                func = _returnProfile2DFunction(field, weight, acc, lazy)
                name = "test%sProfile2D_w%s_l%s_a%s_a%s" % (field,
                                                        weight, lazy,
                                                        accumulation_x,
                                                        accumulation_y)
                setattr(Data3DBase, name, func)

class TestSmoothedCoveringGrid(LagosTestingBase, unittest.TestCase):
    def setUp(self):
        LagosTestingBase.setUp(self)

    def testAllCover(self):
        DIMS = 32
        for i in range(self.hierarchy.max_level+1):
            dx = (DIMS*2**i)**-1
            LE = na.array([0.5,0.5,0.5])-(dx*DIMS/2.0)
            RE = na.array([0.5,0.5,0.5])+(dx*DIMS/2.0)
            cg = self.hierarchy.smoothed_covering_grid(
                    level=i, left_edge=LE, right_edge=RE,
                    dims=[DIMS]*3, fields=["Density"])
            self.assertFalse(na.any(na.isnan(cg["Density"])))
            self.assertFalse(na.any(cg["Density"]==-999))

    def testAddField(self):
        DIMS = 64
        i = 5
        dx = (DIMS*2**i)**-1
        LE = na.array([0.5,0.5,0.5])-(dx*DIMS/2.0)
        RE = na.array([0.5,0.5,0.5])+(dx*DIMS/2.0)
        cg = self.hierarchy.smoothed_covering_grid(
                level=i, left_edge=LE, right_edge=RE,
                dims=[DIMS]*3, fields=["Density"])
        self.assertFalse(na.any(na.isnan(cg["Temperature"])))
        self.assertFalse(na.any(cg["Temperature"]==-999))

        

class TestDataCube(LagosTestingBase, unittest.TestCase):
    def setUp(self):
        LagosTestingBase.setUp(self)

    def testNoGhost(self):
        DW = self.OutputFile["DomainRightEdge"] \
           - self.OutputFile["DomainLeftEdge"]
        for g in self.hierarchy.grids:
            cube = g.retrieve_ghost_zones(0, "Density")
            self.assertTrue(na.all(cube["Density"] == g["Density"]))
            cube["Density"] = na.arange(cube["Density"].size).reshape(cube["Density"].shape)
            cube.flush_data(field="Density")
            self.assertTrue(na.all(g["Density"] == cube["Density"]))

    def testOffsetDomain(self):
        DW = self.OutputFile["DomainRightEdge"] \
           - self.OutputFile["DomainLeftEdge"]
        for g in self.hierarchy.grids:
            cube = self.hierarchy.covering_grid(g.Level,
                g.LeftEdge+DW, g.ActiveDimensions)
            self.assertTrue(na.all(g["Density"] == cube["Density"]))

    def testTwoGhost(self):
        for g in self.hierarchy.grids:
            cube = g.retrieve_ghost_zones(2, "Density")

    def testMultipleFields(self):
        for g in self.hierarchy.grids:
            cube1 = g.retrieve_ghost_zones(0, ["Density","Temperature"])
            self.assertTrue(na.all(cube1["Density"] == g["Density"]))
            self.assertTrue(na.all(cube1["Temperature"] == g["Temperature"]))
            cube2a = g.retrieve_ghost_zones(0, "Density")
            cube2b = g.retrieve_ghost_zones(0, "Temperature")
            self.assertTrue(na.all(cube1["Density"] == cube2a["Density"]))
            self.assertTrue(na.all(cube1["Temperature"] == cube2b["Temperature"]))
    
    def testFlushBackToGrids(self):
        ml = self.hierarchy.max_level
        cg = self.hierarchy.covering_grid(2, [0.0]*3, [64,64,64])
        cg["Ones"] *= 2.0
        cg.flush_data(field="Ones")
        for g in na.concatenate([self.hierarchy.select_grids(i) for i in range(3)]):
            self.assertEqual(g["Ones"].max(), 2.0)
            self.assertEqual(g["Ones"][g["Ones"]*g.child_mask>0].min(), 2.0)

    def testFlushBackToNewCover(self):
        ml = self.hierarchy.max_level
        cg = self.hierarchy.covering_grid(2, [0.0]*3, [64,64,64])
        cg["tempContours"] = cg["Ones"] * 2.0
        cg.flush_data(field="tempContours")
        cg2 = self.hierarchy.covering_grid(2, [0.0]*3, [64,64,64])
        self.assertTrue(na.all(cg["tempContours"] == cg2["tempContours"]))

    def testRawFlushBack(self):
        ml = self.hierarchy.max_level
        cg = self.hierarchy.covering_grid(2, [0.0]*3, [64,64,64])
        cg["DensityNew"] = cg["Density"] * 2.111
        cg.flush_data(field="DensityNew")
        for g in na.concatenate([self.hierarchy.select_grids(i) for i in range(3)]):
            ni = g["DensityNew"] > 0
            min_diff = (g["DensityNew"][ni]/g["Density"][ni]).max()
            max_diff = (g["DensityNew"][ni]/g["Density"][ni]).min()
            min_diff_i = na.argmin(g["DensityNew"][ni]/g["Density"][ni])
            max_diff_i = na.argmax(g["DensityNew"][ni]/g["Density"][ni])
            self.assertAlmostEqual(min_diff, 2.111, 5)
            self.assertAlmostEqual(max_diff, 2.111, 5)

    def testAllCover(self):
        cg = self.hierarchy.covering_grid(1, [0.0]*3, [32,32,32])
        mi, ma = 1e30, -1e30
        for g in na.concatenate([self.hierarchy.select_grids(i) for i in range(2)]):
            ma = max(ma, g["Density"].max())
            mi = min(mi, g["Density"].min())
        self.assertEqual(cg["Density"].max(), ma)
        self.assertEqual(cg["Density"].min(), mi)

    def testCellVolume(self):
        cg = self.hierarchy.covering_grid(2, [0.0]*3, [64,64,64])
        self.assertEqual(na.unique(cg["CellVolume"]).size, 1)


class TestDiskDataType(Data3DBase, DataTypeTestingBase, LagosTestingBase, unittest.TestCase):
    def setUp(self):
        DataTypeTestingBase.setUp(self)
        self.data=self.hierarchy.disk(
                     [0.5,0.5,0.5],[0.2, 0.1, 0.5],1.0,1.0)

class TestRegionDataType(Data3DBase, DataTypeTestingBase, LagosTestingBase, unittest.TestCase):
    def setUp(self):
        DataTypeTestingBase.setUp(self)
        self.data=self.hierarchy.region(
                     [0.5,0.5,0.5],[0.0, 0.0, 0.0],
                     [1.0, 1.0, 1.0])
    def testVolume(self):
        vol = self.data["CellVolume"].sum() / self.data.convert("cm")**3.0
        self.assertAlmostEqual(vol,1.0,7)

class TestRegionStrictDataType(Data3DBase, DataTypeTestingBase, LagosTestingBase, unittest.TestCase):
    def setUp(self):
        DataTypeTestingBase.setUp(self)
        self.data=self.hierarchy.region_strict(
                     [0.5,0.5,0.5],[0.0, 0.0, 0.0],
                     [1.0, 1.0, 1.0])
    def testVolume(self):
        vol = self.data["CellVolume"].sum() / self.data.convert("cm")**3.0
        self.assertAlmostEqual(vol,1.0,7)

class TestPeriodicRegionDataType(Data3DBase, DataTypeTestingBase, LagosTestingBase, unittest.TestCase):
    def setUp(self):
        DataTypeTestingBase.setUp(self)
        self.data=self.hierarchy.periodic_region(
                     [0.5,0.5,0.5],[0.5, 0.5, 0.5],
                     [1.5,1.5,1.5])
    def testVolume(self):
        vol = self.data["CellVolume"].sum() / self.data.convert("cm")**3.0
        self.assertAlmostEqual(vol,1.0,7)

class TestPeriodicRegionStrictDataType(Data3DBase,
            DataTypeTestingBase, LagosTestingBase, unittest.TestCase):
    def setUp(self):
        DataTypeTestingBase.setUp(self)
        self.data=self.hierarchy.periodic_region_strict(
                     [0.5,0.5,0.5],[0.5, 0.5, 0.5],
                     [1.5,1.5,1.5])
    def testVolume(self):
        vol = self.data["CellVolume"].sum() / self.data.convert("cm")**3.0
        self.assertAlmostEqual(vol,1.0,7)

class TestSphereDataType(Data3DBase, DataTypeTestingBase, LagosTestingBase, unittest.TestCase):
    def setUp(self):
        DataTypeTestingBase.setUp(self)
        self.data=self.hierarchy.sphere([0.5,0.5,0.5],1.0)
    def testVolume(self):
        vol = self.data["CellVolume"].sum() / self.data.convert("cm")**3.0
        self.assertAlmostEqual(vol,1.0,7)

class TestSliceDataType(DataTypeTestingBase, LagosTestingBase, unittest.TestCase):
    def setUp(self):
        DataTypeTestingBase.setUp(self)
        self.data = self.hierarchy.slice(0,0.5)

class TestCuttingPlane(DataTypeTestingBase, LagosTestingBase, unittest.TestCase):
    def setUp(self):
        DataTypeTestingBase.setUp(self)
        self.data = self.hierarchy.cutting([0.1,0.3,0.4], [0.5,0.5,0.5], ["Density"])
    def testAxisVectors(self):
        x_v = self.data._x_vec
        y_v = self.data._y_vec
        z_v = self.data._norm_vec
        self.assertAlmostEqual(na.dot(x_v, y_v), 0.0, 7)
        self.assertAlmostEqual(na.dot(x_v, z_v), 0.0, 7)
        self.assertAlmostEqual(na.dot(y_v, z_v), 0.0, 7)
    def testZHeight(self):
        self.assertTrue(na.all(self.data['pz'] < self.data['dx']))

class TestGridDataType(DataTypeTestingBase, LagosTestingBase, unittest.TestCase):
    def setUp(self):
        DataTypeTestingBase.setUp(self)
        self.data = self.hierarchy.grids[0]

class TestExtractFromSphere(TestSphereDataType):
    def setUp(self):
        TestSphereDataType.setUp(self)
        self.region = self.data
        self.ind_to_get = na.where(self.region["Temperature"]>500)
        self.data = self.region.extract_region(self.ind_to_get)

    def testNumberOfEntries(self):
        self.assertEqual(self.ind_to_get[0].shape,
                        self.data["Density"].shape)
    def testVolume(self):
        self.ind_to_get = na.where(self.region["CellVolume"]>0.0)
        vol = self.region.extract_region(self.ind_to_get)["CellVolume"].sum() \
            / self.data.convert("cm")**3.0
        self.assertAlmostEqual(vol,1.0,7)

    def testJoin(self):
        new_region = self.region.extract_region(
                self.region["Temperature"]<=500)
        joined_region = self.data.join(new_region)
        self.assertEqual(joined_region["CellMassMsun"].sum(),
                         self.region["CellMassMsun"].sum())

class TestExtractFromRegion(TestRegionDataType):
    def setUp(self):
        TestRegionDataType.setUp(self)
        self.region = self.data
        self.ind_to_get = na.where(self.region["Temperature"]>500)
        self.data = self.region.extract_region(self.ind_to_get)
    def testNumberOfEntries(self):
        self.assertEqual(self.ind_to_get[0].shape,
                        self.data["Density"].shape)
    def testVolume(self):
        ind_to_get = na.where(self.region["CellVolume"]>0.0)
        vol = self.region.extract_region(ind_to_get)["CellVolume"].sum() \
            / self.data.convert("cm")**3.0
        self.assertAlmostEqual(vol,1.0,7)

    def testJoin(self):
        new_region = self.region.extract_region(
                self.region["Temperature"]<=500)
        joined_region = self.data.join(new_region)
        self.assertEqual(joined_region["CellMassMsun"].sum(),
                         self.region["CellMassMsun"].sum())


class TestUnilinearInterpolator(unittest.TestCase):
    def setUp(self):
        x0, x1 = na.random.uniform(-100,100,2)
        nstep_x = na.random.randint(10,200)
        nvals = na.random.randint(100,1000)

        table = na.mgrid[x0:x1:nstep_x*1j]

        self.ufi_x = yt.lagos.UnilinearFieldInterpolator(table,
                      (x0,x1),'x')
        self.my_dict = {}
        self.my_dict['x'] = na.random.uniform(x0,x1,nvals)

    def testXInt(self):
        nv = self.ufi_x(self.my_dict)
        for i,v in enumerate(nv):
            self.assertAlmostEqual(v, self.my_dict['x'][i], 5)

class TestBilinearInterpolator(unittest.TestCase):
    def setUp(self):
        x0, x1 = na.random.uniform(-100,100,2)
        y0, y1 = na.random.uniform(-100,100,2)
        nstep_x = na.random.randint(10,200)
        nstep_y = na.random.randint(10,200)
        nvals = na.random.randint(100,1000)

        table = na.mgrid[x0:x1:nstep_x*1j,
                         y0:y1:nstep_y*1j]

        self.bfi_x = yt.lagos.BilinearFieldInterpolator(table[0,...],
                      (x0,x1,y0,y1),['x','y'])
        self.bfi_y = yt.lagos.BilinearFieldInterpolator(table[1,...],
                      (x0,x1,y0,y1),['x','y'])
        self.my_dict = {}
        self.my_dict['x'] = na.random.uniform(x0,x1,nvals)
        self.my_dict['y'] = na.random.uniform(y0,y1,nvals)

    def testXInt(self):
        nv = self.bfi_x(self.my_dict)
        for i,v in enumerate(nv):
            self.assertAlmostEqual(v, self.my_dict['x'][i], 5)

    def testYInt(self):
        nv = self.bfi_y(self.my_dict)
        for i,v in enumerate(nv):
            self.assertAlmostEqual(v, self.my_dict['y'][i], 5)

class TestTrilinearInterpolator(unittest.TestCase):
    def setUp(self):
        x0, x1 = na.random.uniform(-100,100,2)
        y0, y1 = na.random.uniform(-100,100,2)
        z0, z1 = na.random.uniform(-100,100,2)
        nstep_x = na.random.randint(10,200)
        nstep_y = na.random.randint(10,200)
        nstep_z = na.random.randint(10,200)
        nvals = na.random.randint(100,1000)

        table = na.mgrid[x0:x1:nstep_x*1j,
                         y0:y1:nstep_y*1j,
                         z0:z1:nstep_z*1j]

        self.tfi_x = yt.lagos.TrilinearFieldInterpolator(table[0,...],
                      (x0,x1,y0,y1,z0,z1),['x','y','z'])
        self.tfi_y = yt.lagos.TrilinearFieldInterpolator(table[1,...],
                      (x0,x1,y0,y1,z0,z1),['x','y','z'])
        self.tfi_z = yt.lagos.TrilinearFieldInterpolator(table[2,...],
                      (x0,x1,y0,y1,z0,z1),['x','y','z'])
        self.my_dict = {}
        self.my_dict['x'] = na.random.uniform(x0,x1,nvals)
        self.my_dict['y'] = na.random.uniform(y0,y1,nvals)
        self.my_dict['z'] = na.random.uniform(z0,z1,nvals)

    def testXInt(self):
        nv = self.tfi_x(self.my_dict)
        for i,v in enumerate(nv):
            self.assertAlmostEqual(v, self.my_dict['x'][i], 5)

    def testYInt(self):
        nv = self.tfi_y(self.my_dict)
        for i,v in enumerate(nv):
            self.assertAlmostEqual(v, self.my_dict['y'][i], 5)

    def testZInt(self):
        nv = self.tfi_z(self.my_dict)
        for i,v in enumerate(nv):
            self.assertAlmostEqual(v, self.my_dict['z'][i], 5)



if __name__ == "__main__":
    unittest.main()
