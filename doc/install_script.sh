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
BRANCH="yt" # This is the branch to which we will forcibly update.

# Here's where you put the HDF5 path if you like; otherwise it'll download it
# and install it on its own
#HDF5_DIR=

# If you need to supply arguments to the NumPy build, supply them here
# This one turns on gfortran manually:
#NUMPY_ARGS="--fcompiler=gnu95"
# If you absolutely can't get the fortran to work, try this:
#NUMPY_ARGS="--fcompiler=fake"

INST_HG=1       # Install Mercurial or not?  If hg is not already
                # installed, yt cannot be installed.
INST_ZLIB=1     # On some systems (Kraken) matplotlib has issues with 
                # the system zlib, which is compiled statically.
                # If need be, you can turn this off.
INST_BZLIB=1    # On some systems, libbzip2 is missing.  This can
                # lead to broken mercurial installations.
INST_PNG=1      # Install a local libpng?  Same things apply as with zlib.
INST_ENZO=0     # Clone a copy of Enzo?

# If you've got YT some other place, set this to point to it.
YT_DIR=""

# If you need to pass anything to matplotlib, do so here:
#MPL_SUPP_LDFLAGS=""
#MPL_SUPP_CCFLAGS=""
#MPL_SUPP_CXXFLAGS=""

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
    if [ "${MYHOST##honest}" != "${MYHOST}" ]
    then
        echo "Looks like you're on Abe."
        echo "We're going to have to set some supplemental environment"
		echo "variables to get this to work..."
		MPL_SUPP_LDFLAGS="${MPL_SUPP_LDFLAGS} -L${DEST_DIR}/lib -L${DEST_DIR}/lib64 -L/usr/local/lib64 -L/usr/local/lib"
    fi
    if [ "${MYHOST##steele}" != "${MYHOST}" ]
    then
        echo "Looks like you're on Steele."
        echo
        echo "NOTE: YOU MUST BE IN THE GNU PROGRAMMING ENVIRONMENT"
        echo "These commands should take care of that for you:"
        echo
        echo "   $ module purge"
        echo "   $ module load gcc"
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
printf "%-15s = %s so I " "INST_ZLIB" "${INST_ZLIB}"
get_willwont ${INST_ZLIB}
echo "be installing zlib"

printf "%-15s = %s so I " "INST_BZLIB" "${INST_BZLIB}"
get_willwont ${INST_BZLIB}
echo "be installing bzlib"

printf "%-15s = %s so I " "INST_HG" "${INST_HG}"
get_willwont ${INST_HG}
echo "be installing Mercurial"

printf "%-15s = %s so I " "INST_ENZO" "${INST_ENZO}"
get_willwont ${INST_ENZO}
echo "be checking out Enzo"

echo

if [ -z "$HDF5_DIR" ]
then
    echo "HDF5_DIR is not set, so I will be installing HDF5"
else
    echo "HDF5_DIR=${HDF5_DIR} , so I will not be installing HDF5"
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
    if [ ! -z `echo $1 | grep h5py` ]
    then
	    echo "${DEST_DIR}/bin/python2.6 setup.py configure --hdf5=${HDF5_DIR}"
	    ( ${DEST_DIR}/bin/python2.6 setup.py configure --hdf5=${HDF5_DIR} 2>&1 ) 1>> ${LOG_FILE} || do_exit
    fi
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
    get_enzotools hdf5-1.6.9.tar.gz
fi

[ $INST_ZLIB -eq 1 ] && get_enzotools zlib-1.2.3.tar.bz2 
[ $INST_BZLIB -eq 1 ] && get_enzotools bzip2-1.0.5.tar.gz
[ $INST_PNG -eq 1 ] && get_enzotools libpng-1.2.43.tar.gz
get_enzotools Python-2.6.3.tgz
get_enzotools numpy-1.5.0.tar.gz
get_enzotools matplotlib-1.0.0.tar.gz
get_enzotools mercurial-1.6.3.tar.gz
get_enzotools ipython-0.10.tar.gz
get_enzotools h5py-1.2.0.tar.gz

if [ $INST_BZLIB -eq 1 ]
then
    if [ ! -e bzip2-1.0.5/done ]
    then
        [ ! -e bzip2-1.0.5 ] && tar xfz bzip2-1.0.5.tar.gz
        echo "Installing BZLIB"
        cd bzip2-1.0.5
        ( make install CFLAGS=-fPIC LDFLAGS=-fPIC PREFIX=${DEST_DIR} 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make -f Makefile-libbz2_so CFLAGS=-fPIC LDFLAGS=-fPIC PREFIX=${DEST_DIR} 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( cp -v libbz2.so.1.0.4 ${DEST_DIR}/lib 2>&1 ) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
    BZLIB_DIR=${DEST_DIR}
    export LDFLAGS="${LDFLAGS} -L${BZLIB_DIR}/lib/ -L${BZLIB_DIR}/lib64/"
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${BZLIB_DIR}/lib/"
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

if [ $INST_PNG -eq 1 ]
then
    if [ ! -e libpng-1.2.43/done ]
    then
        [ ! -e libpng-1.2.43 ] && tar xfz libpng-1.2.43.tar.gz
        echo "Installing PNG"
        cd libpng-1.2.43
        ( ./configure CFLAGS=-I${DEST_DIR}/include --prefix=${DEST_DIR}/ 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
    PNG_DIR=${DEST_DIR}
    export LDFLAGS="${LDFLAGS} -L${PNG_DIR}/lib/ -L${PNG_DIR}/lib64/"
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${PNG_DIR}/lib/"
fi

if [ -z "$HDF5_DIR" ]
then
    if [ ! -e hdf5-1.6.9/done ]
    then
        [ ! -e hdf5-1.6.9 ] && tar xfz hdf5-1.6.9.tar.gz
        echo "Installing HDF5"
        cd hdf5-1.6.9
        ( ./configure --prefix=${DEST_DIR}/ --enable-shared 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
    export HDF5_DIR=${DEST_DIR}
else
    export HDF5_DIR=${HDF5_DIR}
fi
export HDF5_API=16

if [ ! -e Python-2.6.3/done ]
then
    echo "Installing Python.  This may take a while, but don't worry.  YT loves you."
    [ ! -e Python-2.6.3 ] && tar xfz Python-2.6.3.tgz
    cd Python-2.6.3
    ( ./configure --prefix=${DEST_DIR}/ 2>&1 ) 1>> ${LOG_FILE} || do_exit

    ( make 2>&1 ) 1>> ${LOG_FILE} || do_exit
    ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
    touch done
    cd ..
fi

export PYTHONPATH=${DEST_DIR}/lib/python2.6/site-packages/

if [ $INST_HG -eq 1 ]
then
    echo "Installing Mercurial."
    do_setup_py mercurial-1.6.3
    export HG_EXEC=${DEST_DIR}/bin/hg
else
    # We assume that hg can be found in the path.
    if type -P hg &>/dev/null 
    then
        export HG_EXEC=hg
    else
        echo "Cannot find mercurial.  Please set INST_HG=1."
        do_exit
    fi
fi

if [ -z "$YT_DIR" ]
then
    if [ -e $ORIG_PWD/yt/mods.py ]
    then
        YT_DIR="$ORIG_PWD"
    elif [ -e $ORIG_PWD/../yt/mods.py ]
    then
        YT_DIR=`dirname $ORIG_PWD`
    elif [ ! -e yt-hg ] 
    then
        # Note that we clone the entire repository, not just the branch in
        # question.  We update to the correct branch momentarily...
        ( ${HG_EXEC} --debug clone --pull http://hg.enzotools.org/yt/ ./yt-hg 2>&1 ) 1>> ${LOG_FILE}
        YT_DIR="$PWD/yt-hg/"
        ( ${HG_EXEC} up -R ${YT_DIR} -C ${BRANCH} 2>&1 ) 1>> ${LOG_FILE}
    elif [ -e yt-hg ] 
    then
        YT_DIR="$PWD/yt-hg/"
    fi
    echo Setting YT_DIR=${YT_DIR}
fi

# This fixes problems with gfortran linking.
unset LDFLAGS 

echo "Installing distribute"
( ${DEST_DIR}/bin/python2.6 ${YT_DIR}/distribute_setup.py 2>&1 ) 1>> ${LOG_FILE} || do_exit

echo "Installing pip"
( ${DEST_DIR}/bin/easy_install-2.6 pip 2>&1 ) 1>> ${LOG_FILE} || do_exit

do_setup_py numpy-1.5.0 ${NUMPY_ARGS}

if [ -n "${MPL_SUPP_LDFLAGS}" ]
then
    OLD_LDFLAGS=${LDFLAGS}
    export LDFLAGS="${MPL_SUPP_LDFLAGS}"
    echo "Setting LDFLAGS ${LDFLAGS}"
fi
if [ -n "${MPL_SUPP_CXXFLAGS}" ]
then
    OLD_CXXFLAGS=${CXXFLAGS}
    export CXXFLAGS="${MPL_SUPP_CXXFLAGS}"
    echo "Setting CXXFLAGS ${CXXFLAGS}"
fi
if [ -n "${MPL_SUPP_CCFLAGS}" ]
then
    OLD_CCFLAGS=${CCFLAGS}
    export CCFLAGS="${MPL_SUPP_CCFLAGS}"
    echo "Setting CCFLAGS ${CCFLAGS}"
fi
do_setup_py matplotlib-1.0.0
if [ -n "${OLD_LDFLAGS}" ]
then
    export LDFLAG=${OLD_LDFLAGS}
fi
[ -n "${OLD_LDFLAGS}" ] && export LDFLAGS=${OLD_LDFLAGS}
[ -n "${OLD_CXXFLAGS}" ] && export CXXFLAGS=${OLD_CXXFLAGS}
[ -n "${OLD_CCFLAGS}" ] && export CCFLAGS=${OLD_CCFLAGS}
do_setup_py ipython-0.10
do_setup_py h5py-1.2.0

echo "Doing yt update, wiping local changes and updating to branch ${BRANCH}"
MY_PWD=`pwd`
cd $YT_DIR
( ${HG_EXEC} pull && ${HG_EXEC} up -C ${BRANCH} 2>&1 ) 1>> ${LOG_FILE}

echo "Installing yt"
echo $HDF5_DIR > hdf5.cfg
[ $INST_PNG -eq 1 ] && echo $PNG_DIR > png.cfg
( ${DEST_DIR}/bin/python2.6 setup.py develop 2>&1 ) 1>> ${LOG_FILE} || do_exit
touch done
cd $MY_PWD

if [ $INST_ENZO -eq 1 ]
then
    echo "Cloning a copy of Enzo."
    cd ${DEST_DIR}/src/
    ${HG_EXEC} clone https://enzo.googlecode.com/hg/ ./enzo-hg-stable
    cd $MY_PWD
fi

echo
echo
echo "========================================================================"
echo
echo "yt is now installed in $DEST_DIR ."
echo "To run from this new installation, the a few variables need to be"
echo "prepended with the following information:"
echo
echo "YT_DEST         => $DEST_DIR"
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
if [ $INST_ENZO -eq 1 ]
then
  echo "Enzo has also been checked out, but not built."
  echo
  echo "$DEST_DIR/src/enzo-hg-stable"
  echo
  echo "The value of YT_DEST can be used as an HDF5 installation location."
  echo "Questions about Enzo should be directed to the Enzo User List."
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
echo
echo "Oh, look at me, still talking when there's science to do!"
echo "Good luck, and email the user list if you run into any problems."
