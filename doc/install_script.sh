# Modify this line and uncomment to set where we install
#MY_DIR=$HOME/yt-`uname -p`/

if [ -z "$MY_DIR" ]
then
    echo "Edit this script, set the MY_DIR parameter and re-run."
    exit 1
fi

mkdir -p $MY_DIR/src
cd $MY_DIR/src

# Individual processes
[ ! -e hdf5-1.6.7.tar.gz ] && wget ftp://ftp.hdfgroup.org/HDF5/current16/src/hdf5-1.6.7.tar.gz
[ ! -e Python-2.5.2.tgz ] && wget http://python.org/ftp/python/2.5.2/Python-2.5.2.tgz
[ ! -e pytables-2.0.4.tar.gz ] && wget http://www.pytables.org/download/stable/pytables-2.0.4.tar.gz
[ ! -e matplotlib-0.91.4.tar.gz ] && wget "http://downloads.sourceforge.net/matplotlib/matplotlib-0.91.4.tar.gz"
[ ! -e numpy-1.0.4.tar.gz ] && wget "http://downloads.sourceforge.net/numpy/numpy-1.0.4.tar.gz?modtime=1194536674&big_mirror=0"
[ ! -e wxPython-src-2.8.7.1.tar.bz2 ] && wget http://downloads.sourceforge.net/wxpython/wxPython-src-2.8.7.1.tar.bz2
[ ! -e yt-0.3.tar.gz ] && wget http://yt.enzotools.org/files/yt-0.3.tar.gz

if [ ! -e hdf5-1.6.7/done ]
then
    [ ! -e hdf5-1.6.7 ] && tar xvfz hdf5-1.6.7.tar.gz
    echo "Doing HDF5"
    cd hdf5-1.6.7
    ./configure --prefix=$MY_DIR/ || exit 1
    make install || exit 1
    touch done
    cd ..
fi

if [ ! -e Python-2.5.2/done ]
then
    echo "Doing Python"
    [ ! -e Python-2.5.2 ] && tar xvfz Python-2.5.2.tgz
    cd Python-2.5.2
    ./configure --prefix=$MY_DIR/ || exit 1
    make || exit 1
    make install || exit 1
    touch done
    cd ..
fi

export PYTHONPATH=$MY_DIR/lib/python2.5/site-packages/

if [ ! -e numpy-1.0.4/done ]
then
    echo "Doing NumPy"
    [ ! -e numpy-1.0.4 ] && tar xvfz numpy-1.0.4.tar.gz
    cd numpy-1.0.4
    $MY_DIR/bin/python2.5 setup.py install || exit 1
    touch done
    cd ..
fi

if [ ! -e pytables-2.0.4/done ]
then
    echo "Doing PyTables"
    [ ! -e pytables-2.0.4 ] && tar xvfz pytables-2.0.4.tar.gz
    cd pytables-2.0.4
    $MY_DIR/bin/python2.5 setup.py install --hdf5=$MY_DIR/ || exit 1
    touch done
    cd ..
fi

if [ ! -e matplotlib-0.91.4/done ]
then
    echo "Doing Matplotlib"
    [ ! -e matplotlib-0.91.4 ] && tar xvfz matplotlib-0.91.4.tar.gz
    cd matplotlib-0.91.4
    $MY_DIR/bin/python2.5 setup.py install || exit 1
    touch done
    cd ..
fi

if [ ! -e wxPython-src-2.8.7.1/done ]
then
    echo "Doing wxPython"
    [ ! -e wxPython-src-2.8.7.1 ] && tar xvfj wxPython-src-2.8.7.1.tar.bz2
    cd wxPython-src-2.8.7.1
    ./configure --prefix=$MY_DIR/ --with-opengl || exit 1
    make install || exit 1
    cd contrib
    make install || exit 1
    cd ../wxPython/
    $MY_DIR/bin/python2.5 setup.py WX_CONFIG=$MY_DIR/bin/wx-config install || exit 1
    touch ../done
    cd ../..
fi

if [ ! -e yt-0.3/done ]
then
    echo "Doing yt installation"
    [ ! -e yt-0.3 ] && tar xvfz yt-0.3.tar.gz
    cd yt-0.3
    echo $MY_DIR > hdf5.cfg
    $MY_DIR/bin/python2.5 setup.py install || exit 1
    touch done
    cd ..
fi

echo "yt is now installed in $MY_DIR ."
echo "To run from this new installation, the a few variables need to be"
echo "prepended with the following information:"
echo
echo "PATH            => $MY_DIR/bin/"
echo "PYTHONPATH      => $MY_DIR/lib/python2.5/site-packages/"
echo "LD_LIBRARY_PATH => $MY_DIR/lib/"
echo
echo "You can get a fully-loaded yt prompt by running:"
echo "$MY_DIR/bin/yt"
echo
