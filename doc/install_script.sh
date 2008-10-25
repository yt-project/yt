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

#DEST_DIR="`pwd`/yt-`uname -p`"   # Installation location

# Here's where you put the HDF5 path if you like
#HDF5_DIR=


INST_WXPYTHON=0 # If you 't want to install wxPython, set this to 1
INST_ZLIB=0     # On some systems (Kraken) matplotlib has issues with 
                # the system zlib, which is compiled statically

# If you've got YT some other place, set this to point to it.
YT_DIR=""

#------------------------------------------------------------------------------#
#                                                                              #
# Okay, the script starts here.  Feel free to play with it, but hopefully      #
# it'll work as is.                                                            #
#                                                                              #
#------------------------------------------------------------------------------#

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
[ $INST_WXPYTHON -eq 1 ] && [ ! -e wxPython-src-2.8.7.1.tar.bz2 ] && wget http://downloads.sourceforge.net/wxpython/wxPython-src-2.8.7.1.tar.bz2

if [ -z "$YT_DIR" ]
then
    if [ -e $ORIG_PWD/yt/mods.py ]
    then
        YT_DIR="$ORIG_PWD"
    elif [ -e $ORIG_PWD/../yt/mods.py ]
    then
        YT_DIR=`dirname $ORIG_PWD`
    elif [ ! -e yt-1.0-svn ] 
    then
        svn co http://svn.enzotools.org/yt/branches/yt-1.0/ ./yt-1.0-svn
        YT_DIR="$PWD/yt-1.0-svn/"
    elif [ -e yt-1.0-svn ]
    then
        YT_DIR="$PWD/yt-1.0-svn/"
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

echo "Installing numpy, matplotlib, ipython"
${DEST_DIR}/bin/python2.5 ${YT_DIR}/ez_setup.py
${DEST_DIR}/bin/easy_install numpy      || exit 1
${DEST_DIR}/bin/easy_install matplotlib || exit 1
${DEST_DIR}/bin/easy_install ipython    || exit 1

echo "Doing yt update"
MY_PWD=`pwd`
cd $YT_DIR
svn up
echo $HDF5_DIR > hdf5.cfg
${DEST_DIR}/bin/python2.5 setup.py develop || exit 1
touch done
cd $MY_PWD

echo
echo
echo "========================================================================"
echo
echo "yt is now installed in $DEST_DIR ."
echo "To run from this new installation, the a few variables need to be"

echo "prepended with the following information:"
echo
echo "PATH            => $DEST_DIR/bin/"
echo "PYTHONPATH      => $DEST_DIR/lib/python2.5/site-packages/"
echo "LD_LIBRARY_PATH => $DEST_DIR/lib/"

echo
echo "For interactive data analysis and visualization, we recommend running"
echo "the IPython interface, which will become more fully featured with time:"
echo
echo "$DEST_DIR/bin/iyt"
echo
echo "Note of interest: this installation will use the directory"
echo "$YT_DIR"
echo "as the source for all the YT code.  This means you probably shouldn't"
echo "delete it, but on the plus side, any changes you make there are"
echo "automatically propagated."
echo
echo "========================================================================"
