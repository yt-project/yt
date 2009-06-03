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

DEST_SUFFIX="yt-`uname -p`"
DEST_DIR="`pwd`/${DEST_SUFFIX/ /}"   # Installation location

# Here's where you put the HDF5 path if you like; otherwise it'll download it
# and install it on its own
#HDF5_DIR=

# If you need to supply arguments to the NumPy build, supply them here
# This one turns on gfortran manually:
#NUMPY_ARGS="--fcompiler=gnu95"
# If you absolutely can't get the fortran to work, try this:
#NUMPY_ARGS="--fcompiler=fake"

INST_WXPYTHON=0 # If you 't want to install wxPython, set this to 1
INST_ZLIB=1     # On some systems (Kraken) matplotlib has issues with 
                # the system zlib, which is compiled statically.
                # If need be, you can turn this off.
INST_HG=0       # Install Mercurial or not?

# If you've got YT some other place, set this to point to it.
YT_DIR=""

#------------------------------------------------------------------------------#
#                                                                              #
# Okay, the script starts here.  Feel free to play with it, but hopefully      #
# it'll work as is.                                                            #
#                                                                              #
#------------------------------------------------------------------------------#

LOG_FILE="${DEST_DIR}/yt_install.log"

function get_willwont
{
    if [ $1 -eq 1 ]
    then
        echo -n "will  "
    else
        echo -n "won't "
    fi
}

function host_specific
{
    MYHOST=`hostname -s`  # just give the short one, not FQDN
    MYHOSTLONG=`hostname` # FQDN, for Ranger
    if [ "${MYHOST##kraken}" != "${MYHOST}" ]
    then
        echo "Looks like you're on Kraken."
        echo
        echo "NOTE: YOU MUST BE IN THE GNU PROGRAMMING ENVIRONMENT"
        echo "   $ module swap PrgEnv-pgi PrgEnv-gnu"
        echo
        return
    fi
    if [ "${MYHOST##verne}" != "${MYHOST}" ]
    then
        echo "Looks like you're on Verne."
        echo
        echo "NOTE: YOU MUST BE IN THE GNU PROGRAMMING ENVIRONMENT"
        echo "This command will take care of that for you:"
        echo
        echo "   $ module swap PE-pgi PE-gnu"
        echo
    fi
    if [ "${MYHOSTLONG%%ranger.tacc.utexas.edu}" != "${MYHOSTLONG}" ]
    then
        echo "Looks like you're on Ranger."
        echo
        echo "NOTE: YOU MUST BE IN THE GNU PROGRAMMING ENVIRONMENT"
        echo "These commands should take care of that for you:"
        echo
        echo "   $ module unload mvapich-devel"
        echo "   $ module swap pgi gcc"
        echo "   $ module load mvapich-devel"
        echo
    fi
}


echo
echo
echo "========================================================================"
echo
echo "Hi there!  This is the YT installation script.  We're going to download"
echo "some stuff and install it to create a self-contained, isolated"
echo "environment for YT to run within."
echo
echo "Inside the installation script you can set a few variables.  Here's what"
echo "they're currently set to -- you can hit Ctrl-C and edit the values in "
echo "the script if you aren't such a fan."
echo
printf "%-15s = %s so I " "INST_WXPYTHON" "${INST_WXPYTHON}"
get_willwont $INST_WXPYTHON
echo "be installing wxPython"

printf "%-15s = %s so I " "INST_ZLIB" "${INST_ZLIB}"
get_willwont ${INST_ZLIB}
echo "be installing zlib"

printf "%-15s = %s so I " "INST_HG" "${INST_HG}"
get_willwont ${INST_HG}
echo "be installing Mercurial"

echo

if [ -z "$HDF5_DIR" ]
then
    echo "HDF5_DIR is not set, so I will be installing HDF5"
else
    echo "HDF5_DIR=${HDF5_DIR} , so I will not be installing HDF5"
fi

if [ -z "$YT_DIR" ]
then
    echo "YT_DIR is not set, so I will be checking out a fresh copy"
else
    echo "YT_DIR=${YT_DIR} , so I will use that for YT"
fi

echo
echo "Installation will be to"
echo "  ${DEST_DIR}"
echo
echo "and I'll be logging the installation in"
echo "  ${LOG_FILE}"
echo
echo "I think that about wraps it up.  If you want to continue, hit enter.  "
echo "If you'd rather stop, maybe think things over, even grab a sandwich, "
echo "hit Ctrl-C."
echo
host_specific
echo "========================================================================"
echo
read -p "[hit enter] "
echo
echo "Awesome!  Here we go."
echo

function do_exit
{
    echo "Failure.  Check ${LOG_FILE}."
    exit 1
}

function do_setup_py
{
    [ -e $1/done ] && return
    echo "Installing $1 (arguments: '$*')"
    [ ! -e $1 ] && tar xfz $1.tar.gz
    cd $1
    shift
    ( ${DEST_DIR}/bin/python2.6 setup.py build   $* 2>&1 ) 1>> ${LOG_FILE} || do_exit
    ( ${DEST_DIR}/bin/python2.6 setup.py install    2>&1 ) 1>> ${LOG_FILE} || do_exit
    touch done
    cd ..
}

function get_enzotools
{
    echo "Downloading $1 from yt.enzotools.org"
    [ -e $1 ] && return
    wget -nv "http://yt.enzotools.org/dependencies/$1" || do_exit
    wget -nv "http://yt.enzotools.org/dependencies/$1.md5" || do_exit
    ( which md5sum &> /dev/null ) || return # return if we don't have md5sum
    ( md5sum -c $1.md5 2>&1 ) 1>> ${LOG_FILE} || do_exit
}

ORIG_PWD=`pwd`

if [ -z "${DEST_DIR}" ]
then
    echo "Edit this script, set the DEST_DIR parameter and re-run."
    exit 1
fi

mkdir -p ${DEST_DIR}/src
cd ${DEST_DIR}/src

# Individual processes
if [ -z "$HDF5_DIR" ]
then
    echo "Downloading HDF5"
    get_enzotools hdf5-1.6.8.tar.gz
fi

[ $INST_ZLIB -eq 1 ] && get_enzotools zlib-1.2.3.tar.bz2 
[ $INST_WXPYTHON -eq 1 ] && get_enzotools wxPython-src-2.8.7.1.tar.bz2
get_enzotools Python-2.6.1.tgz
get_enzotools numpy-1.2.1.tar.gz
get_enzotools matplotlib-0.98.5.2.tar.gz
get_enzotools ipython-0.9.1.tar.gz
get_enzotools h5py-1.1.0.tar.gz

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
        ( svn co http://svn.enzotools.org/yt/trunk/ ./yt-trunk-svn 2>&1 ) 1>> ${LOG_FILE}
        YT_DIR="$PWD/yt-trunk-svn/"
    elif [ -e yt-trunk-svn ] 
    then
        YT_DIR="$PWD/yt-trunk-svn/"
    fi
    echo Setting YT_DIR=${YT_DIR}
fi

if [ $INST_ZLIB -eq 1 ]
then
    if [ ! -e zlib-1.2.3/done ]
    then
        [ ! -e zlib-1.2.3 ] && tar xfj zlib-1.2.3.tar.bz2
        echo "Installing ZLIB"
        cd zlib-1.2.3
        ( ./configure --shared --prefix=${DEST_DIR}/ 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
    ZLIB_DIR=${DEST_DIR}
    export LDFLAGS="${LDFLAGS} -L${ZLIB_DIR}/lib/ -L${ZLIB_DIR}/lib64/"
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${ZLIB_DIR}/lib/"
fi

if [ -z "$HDF5_DIR" ]
then
    if [ ! -e hdf5-1.6.8/done ]
    then
        [ ! -e hdf5-1.6.8 ] && tar xfz hdf5-1.6.8.tar.gz
        echo "Installing HDF5"
        cd hdf5-1.6.8
        ( ./configure --prefix=${DEST_DIR}/ --enable-shared 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
    export HDF5_DIR=${DEST_DIR}
fi
export HDF5_API=16

if [ ! -e Python-2.6.1/done ]
then
    echo "Installing Python.  This may take a while, but don't worry.  YT loves you."
    [ ! -e Python-2.6.1 ] && tar xfz Python-2.6.1.tgz
    cd Python-2.6.1
    ( ./configure --prefix=${DEST_DIR}/ 2>&1 ) 1>> ${LOG_FILE} || do_exit

    ( make 2>&1 ) 1>> ${LOG_FILE} || do_exit
    ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
    touch done
    cd ..
fi

export PYTHONPATH=${DEST_DIR}/lib/python2.6/site-packages/

if [ $INST_WXPYTHON -eq 1 ] && [ ! -e wxPython-src-2.8.7.1/done ]
then
    echo "Installing wxPython.  This may take a while, but don't worry.  YT loves you."
    [ ! -e wxPython-src-2.8.7.1 ] && tar xfj wxPython-src-2.8.7.1.tar.bz2
    cd wxPython-src-2.8.7.1

    ( ./configure --prefix=${DEST_DIR}/ --with-opengl 2>&1 ) 1>> ${LOG_FILE} || do_exit
    ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
    cd contrib
    ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
    cd ../wxPython/
    ( ${DEST_DIR}/bin/python2.6 setup.py WX_CONFIG=${DEST_DIR}/bin/wx-config install 2>&1 ) 1>> ${LOG_FILE} || do_exit

    touch ../done
    cd ../..
fi

# This fixes problems with gfortran linking.
unset LDFLAGS 

echo "Installing setuptools"
( ${DEST_DIR}/bin/python2.6 ${YT_DIR}/ez_setup.py 2>&1 ) 1>> ${LOG_FILE} || do_exit

do_setup_py numpy-1.2.1 ${NUMPY_ARGS}
do_setup_py matplotlib-0.98.5.2
do_setup_py ipython-0.9.1
do_setup_py h5py-1.1.0

echo "Doing yt update"
MY_PWD=`pwd`
cd $YT_DIR
( svn up 2>&1 ) 1>> ${LOG_FILE}

echo "Installing yt"
echo $HDF5_DIR > hdf5.cfg
( ${DEST_DIR}/bin/python2.6 setup.py develop 2>&1 ) 1>> ${LOG_FILE} || do_exit
touch done
cd $MY_PWD

if [ $INST_HG -eq 1 ]
then
    echo "Installing Mercurial."
    ( ${DEST_DIR}/bin/easy_install-2.6 mercurial 2>&1 ) 1>> ${LOG_FILE} || do_exit
fi

echo
echo
echo "========================================================================"
echo
echo "yt is now installed in $DEST_DIR ."
echo "To run from this new installation, the a few variables need to be"
echo "prepended with the following information:"
echo
echo "PATH            => $DEST_DIR/bin/"
echo "PYTHONPATH      => $DEST_DIR/lib/python2.6/site-packages/"
echo "LD_LIBRARY_PATH => $DEST_DIR/lib/"
echo
echo "For interactive data analysis and visualization, we recommend running"
echo "the IPython interface, which will become more fully featured with time:"
echo
echo "$DEST_DIR/bin/iyt"
echo
echo "For command line analysis run:"
echo
echo "$DEST_DIR/bin/yt"
echo
echo "Note of interest: this installation will use the directory"
echo "$YT_DIR"
echo "as the source for all the YT code.  This means you probably shouldn't"
echo "delete it, but on the plus side, any changes you make there are"
echo "automatically propagated."
if [ $INST_HG -eq 1 ]
then
  echo "Mercurial has also been installed:"
  echo
  echo "$DEST_DIR/bin/hg"
  echo
fi
echo
echo "For support, see one of the following websites:"
echo
echo "    http://yt.enzotools.org/wiki/"
echo "    http://yt.enzotools.org/doc/"
echo
echo "Or join the mailing list:"
echo 
echo "    http://lists.spacepope.org/listinfo.cgi/yt-users-spacepope.org"
echo
echo "========================================================================"
