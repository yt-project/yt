#
# Hi there!  Welcome to the yt installation script.
#
# This script is designed to create a fully isolated Python installation
# with the dependencies you need to run yt.
#
# There are a few options, but you only need to set *one* of them.  And
# that's the next one, DEST_DIR.  But, if you want to use an existing HDF5
# installation you can set HDF5_DIR, or if you want to use some other
# subversion checkout of YT, you can set YT_DIR, too.  (It'll already
# check the current directory and one up.
#
# NOTE: If you have trouble with wxPython, set INST_WXPYTHON=0 .
#
# And, feel free to drop me a line: matthewturk@gmail.com
#

#DEST_DIR="$HOME/local/"

# Here's where you put the HDF5 path if you like
#HDF5_DIR=
# If you don't want to install wxPython, turn this to zero
INST_WXPYTHON=1
INST_IPYTHON=1
INST_ZLIB=0
# If you've got YT some other place, set this to point to it.
#YT_DIR=

#
# Okay, the script starts here.  Feel free to play with it, but hopefully it'll
# work as is.
#

ORIG_PWD=`pwd`

if [ -z "${DEST_DIR}" ]
then
    echo "Edit this script, set the DEST_DIR parameter and re-run."
    exit 1
fi

echo "Installing into ${DEST_DIR}"
echo "INST_WXPYTHON=${INST_WXPYTHON}"
echo "INST_IPYTHON=${INST_IPYTHON}"

mkdir -p ${DEST_DIR}/src
cd ${DEST_DIR}/src

# Individual processes
if [ -z "$HDF5_DIR" ]
then
    [ ! -e hdf5-1.6.7.tar.gz ] && wget ftp://ftp.hdfgroup.org/HDF5/current16/src/hdf5-1.6.7.tar.gz
fi

[ $INST_ZLIB -eq 1 ] && [ ! -e zlib-1.2.3.tar.bz2 ] && wget http://www.zlib.net/zlib-1.2.3.tar.bz2
[ ! -e Python-2.5.2.tgz ] && wget http://python.org/ftp/python/2.5.2/Python-2.5.2.tgz
[ ! -e pytables-2.0.4.tar.gz ] && wget http://www.pytables.org/download/stable/pytables-2.0.4.tar.gz
[ ! -e matplotlib-0.91.4.tar.gz ] && wget "http://downloads.sourceforge.net/matplotlib/matplotlib-0.91.4.tar.gz"
[ ! -e numpy-1.0.4.tar.gz ] && wget "http://downloads.sourceforge.net/numpy/numpy-1.0.4.tar.gz?modtime=1194536674&big_mirror=0"
[ $INST_WXPYTHON -eq 1 ] && [ ! -e wxPython-src-2.8.7.1.tar.bz2 ] && wget http://downloads.sourceforge.net/wxpython/wxPython-src-2.8.7.1.tar.bz2
[ $INST_IPYTHON -eq 1 ] && [ ! -e ipython-0.9.1.tar.gz ] && wget http://ipython.scipy.org/dist/ipython-0.9.1.tar.gz

if [ -z "$YT_DIR" ]
then
    if [ -e $ORIG_PWD/yt/mods.py ]
    then
        YT_DIR="$ORIG_PWD"
    elif [ -e $ORIG_PWD/../yt/mods.py ]
    then
        YT_DIR=`dirname $ORIG_PWD`
    elif [ ! -e yt-trunk-svn ] 
    then
        svn co http://svn.enzotools.org/yt/trunk/ ./yt-trunk-svn
        YT_DIR="$PWD/yt-trunk-svn/"
    fi
    echo Setting YT_DIR=${YT_DIR}
fi

if [ $INST_ZLIB -eq 1 ]
then
    if [ ! -e zlib-1.2.3/done ]
    then
        [ ! -e zlib-1.2.3 ] && tar xvfj zlib-1.2.3.tar.bz2
        echo "Doing ZLIB"
        cd zlib-1.2.3
        ./configure --shared --prefix=${DEST_DIR}/ || exit 1
        make install || exit 1
        touch done
        cd ..
    fi
    ZLIB_DIR=${DEST_DIR}
    LDFLAGS="${LDFLAGS} -L${ZLIB_DIR}/lib/"
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${ZLIB_DIR}/lib/"
fi

if [ -z "$HDF5_DIR" ]
then
    if [ ! -e hdf5-1.6.7/done ]
    then
        [ ! -e hdf5-1.6.7 ] && tar xvfz hdf5-1.6.7.tar.gz
        echo "Doing HDF5"
        cd hdf5-1.6.7
        ./configure --prefix=${DEST_DIR}/ || exit 1
        make install || exit 1
        touch done
        cd ..
    fi
    HDF5_DIR=${DEST_DIR}
fi

if [ ! -e Python-2.5.2/done ]
then
    echo "Doing Python"
    [ ! -e Python-2.5.2 ] && tar xvfz Python-2.5.2.tgz
    cd Python-2.5.2
    ./configure --prefix=${DEST_DIR}/ || exit 1

    make || exit 1
    make install || exit 1
    touch done
    cd ..
fi

export PYTHONPATH=${DEST_DIR}/lib/python2.5/site-packages/

if [ ! -e numpy-1.0.4/done ]
then
    echo "Doing NumPy"

    [ ! -e numpy-1.0.4 ] && tar xvfz numpy-1.0.4.tar.gz
    cd numpy-1.0.4
    ${DEST_DIR}/bin/python2.5 setup.py install || exit 1
    touch done
    cd ..
fi

if [ ! -e pytables-2.0.4/done ]

then
    echo "Doing PyTables"
    [ ! -e pytables-2.0.4 ] && tar xvfz pytables-2.0.4.tar.gz
    cd pytables-2.0.4
    ${DEST_DIR}/bin/python2.5 setup.py install --hdf5=$HDF5_DIR/ || exit 1
    touch done

    cd ..
fi

if [ $INST_WXPYTHON -eq 1 ] && [ ! -e wxPython-src-2.8.7.1/done ]
then
    echo "Doing wxPython"
    [ ! -e wxPython-src-2.8.7.1 ] && tar xvfj wxPython-src-2.8.7.1.tar.bz2
    cd wxPython-src-2.8.7.1

    ./configure --prefix=${DEST_DIR}/ --with-opengl || exit 1
    make install || exit 1
    cd contrib
    make install || exit 1
    cd ../wxPython/
    ${DEST_DIR}/bin/python2.5 setup.py WX_CONFIG=${DEST_DIR}/bin/wx-config install || exit 1

    touch ../done
    cd ../..
fi

if [ ! -e matplotlib-0.91.4/done ]
then
    echo "Doing Matplotlib"
    [ ! -e matplotlib-0.91.4 ] && tar xvfz matplotlib-0.91.4.tar.gz
    cd matplotlib-0.91.4
    ${DEST_DIR}/bin/python2.5 setup.py install || exit 1

    touch done
    cd ..
fi

if [ $INST_IPYTHON -eq 1 ] && [ ! -e ipython-0.9.1/done ]
then
    echo "Doing IPython"
    [ ! -e ipython-0.9.1 ] && tar xvfz ipython-0.9.1.tar.gz
    cd ipython-0.9.1
    ${DEST_DIR}/bin/python2.5 setup.py install || exit 1
    touch done
    cd ..
fi

echo "Doing yt update"
MY_PWD=`pwd`
cd $YT_DIR
svn up
echo $HDF5_DIR > hdf5.cfg
${DEST_DIR}/bin/python2.5 ez_setup.py
${DEST_DIR}/bin/python2.5 setup.py install || exit 1
touch done
cd $MY_PWD

echo "yt is now installed in $DEST_DIR ."
echo "To run from this new installation, the a few variables need to be"

echo "prepended with the following information:"
echo
echo "PATH            => $DEST_DIR/bin/"
echo "PYTHONPATH      => $DEST_DIR/lib/python2.5/site-packages/"
echo "LD_LIBRARY_PATH => $DEST_DIR/lib/"

echo
echo "You can get a fully-loaded yt prompt by running:"
echo "$DEST_DIR/bin/yt"
echo
if [ $INST_IPYTHON -eq 1 ]
then
    echo "For interactive data analysis and visualization, we now recommend running"
    echo "the IPython interface, which will become more fully featured with time!"
    echo
    echo "$DEST_DIR/bin/iyt"
    echo
fi
