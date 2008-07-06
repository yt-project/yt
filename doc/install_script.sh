# Modify this line to change where we install
MY_DIR=$HOME/local
#MY_DIR=$HOME/ki02/testing_eric
MY_DIR=/data/mturk/testing_eric


mkdir -p $MY_DIR/src
cd $MY_DIR/src

# Individual processes
[ ! -e hdf5-1.8.1.tar.gz ] && wget ftp://ftp.hdfgroup.org/HDF5/current/src/hdf5-1.8.1.tar.gz
[ ! -e Python-2.5.2.tgz ] && wget http://python.org/ftp/python/2.5.2/Python-2.5.2.tgz
[ ! -e pytables-2.0.4.tar.gz ] && wget http://www.pytables.org/download/stable/pytables-2.0.4.tar.gz
[ ! -e matplotlib-0.91.2.tar.gz ] && wget "http://downloads.sourceforge.net/matplotlib/matplotlib-0.91.2.tar.gz?modtime=1199627250&big_mirror=0"
[ ! -e numpy-1.0.4.tar.gz ] && wget "http://downloads.sourceforge.net/numpy/numpy-1.0.4.tar.gz?modtime=1194536674&big_mirror=0"
[ ! -e wxPython-src-2.8.7.1.tar.bz2 ] && wget http://downloads.sourceforge.net/wxpython/wxPython-src-2.8.7.1.tar.bz2

[ ! -e hdf5-1.8.1/done ] && tar xvfz hdf5-1.8.1.tar.gz
echo "Doing HDF5"
cd hdf5-1.8.1
./configure --prefix=$MY_DIR/ || exit 1
make install || exit 1
touch done
cd ..

[ ! -e Python-2.5.2/done ] && tar xvfz Python-2.5.2.tgz
echo "Doing Python"
cd Python-2.5.2
./configure --prefix=$MY_DIR/ || exit 1
make || exit 1
make install || exit 1
touch done
cd ..

[ ! -e numpy-1.0.4/done ] && tar xvfz numpy-1.0.4.tar.gz
echo "Doing NumPy"
cd numpy-1.0.4
$MY_DIR/bin/python2.5 setup.py install || exit 1
touch done
cd ..

[ ! -e pytables-2.0.4/done ] && tar xvfz pytables-2.0.4.tar.gz
echo "Doing PyTables"
cd pytables-2.0.4
$MY_DIR/bin/python2.5 setup.py install --hdf5=$MY_DIR/ || exit 1
touch done
cd ..

[ ! -e matplotlib-0.91.2/done ] && tar xvfz matplotlib-0.91.2.tar.gz
echo "Doing Matplotlib"
cd matplotlib-0.91.2
$MY_DIR/bin/python2.5 setup.py install || exit 1
touch done
cd ..

[ ! -e wxPython-src-2.8.7.1/done ] && tar xvfj wxPython-src-2.8.7.1.tar.bz2
echo "Doing wxPython"
cd wxPython-src-2.8.7.1
./configure --prefix=$MY_DIR/ || exit 1
make
cd contrib
make
cd ..
make install || exit 1
cd wxPython
$MY_DIR/bin/python2.5 setup.py install || exit 1
touch done
cd ../..


