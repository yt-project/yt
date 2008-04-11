import unittest, numpy, tables, sys, os, os.path
sys.path.insert(0,".")

from new import classobj

from yt.lagos.HDF5LightReader import ReadData, ReadingError

my_dtypes = ['short','int','longlong','ushort','uint','ulonglong',
             'float','double']

class HDF5LightTestIOBase(object):
    def setUp(self):
        self.rand_array = numpy.random.random(3000).reshape((30,10,10)).astype(self.dtype)
    def test_check_io(self):
        my_table = tables.openFile("testing_h5lt_io.h5","w")
        my_table.createArray("/","%s" % (self.dtype),self.rand_array)
        my_table.close()
        recv_array = ReadData("testing_h5lt_io.h5", "/%s" % (self.dtype))
        self.assert_(numpy.all(recv_array == self.rand_array))
        self.assert_(recv_array.shape == self.rand_array.shape)
    def tearDown(self):
        os.unlink("testing_h5lt_io.h5")

for dtype in my_dtypes:
    temp = classobj("TestingIO_%s" % (dtype),
            (HDF5LightTestIOBase,unittest.TestCase), {'dtype':dtype})
    exec('TestingIO_%s = temp' % dtype)

class HDF5LightTestError(unittest.TestCase):
    def test_no_file(self):
        fn = "%s.h5" % int(numpy.random.random(1) * 1e6)
        self.assertRaises(ReadingError, ReadData,fn,"/Nothing")
    def test_no_dataset(self):
        fn = "%s.h5" % int(numpy.random.random(1) * 1e6)
        my_table = tables.openFile("testing_h5lt_io.h5","w")
        my_table.close()
        self.assertRaises(ReadingError, ReadData,fn,"/Nothing")
    def tearDown(self):
        if os.path.exists("testing_h5lt_io.h5"): os.unlink("testing_h5lt_io.h5")

if __name__ == "__main__":
    unittest.main()
