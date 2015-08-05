#
# Hi there!  Welcome to the yt installation script.
#
# First things first, if you experience problems, please visit the Help 
# section at http://yt-project.org.
#
# This script is designed to create a fully isolated Python installation
# with the dependencies you need to run yt.
#
# There are a few options, but you only need to set *one* of them, which is 
# the next one, DEST_DIR:

DEST_SUFFIX="yt-`uname -m`"
DEST_DIR="`pwd`/${DEST_SUFFIX/ /}"   # Installation location
BRANCH="yt" # This is the branch to which we will forcibly update.

if [ ${REINST_YT} ] && [ ${REINST_YT} -eq 1 ] && [ -n ${YT_DEST} ]
then
    DEST_DIR=${YT_DEST}
fi

# What follows are some other options that you may or may not need to change.

# Here's where you put the HDF5 path if you like; otherwise it'll download it
# and install it on its own
#HDF5_DIR=

# If you've got yt some other place, set this to point to it. The script will
# already check the current directory and the one above it in the tree.
YT_DIR=""

# If you need to supply arguments to the NumPy or SciPy build, supply them here
# This one turns on gfortran manually:
#NUMPY_ARGS="--fcompiler=gnu95"
# If you absolutely can't get the fortran to work, try this:
#NUMPY_ARGS="--fcompiler=fake"

INST_PY3=0      # Install Python 3 along with Python 2. If this is turned
                # on, all Python packages (including yt) will be installed
                # in Python 3 (except Mercurial, which requires Python 2).
INST_HG=1       # Install Mercurial or not?  If hg is not already
                # installed, yt cannot be installed.
INST_ZLIB=1     # On some systems (Kraken) matplotlib has issues with
                # the system zlib, which is compiled statically.
                # If need be, you can turn this off.
INST_BZLIB=1    # On some systems, libbzip2 is missing.  This can
                # lead to broken mercurial installations.
INST_PNG=1      # Install a local libpng?  Same things apply as with zlib.
INST_FTYPE=1    # Install FreeType2 locally?
INST_SQLITE3=1  # Install a local version of SQLite3?
INST_PYX=0      # Install PyX?  Sometimes PyX can be problematic without a
                # working TeX installation.
INST_0MQ=1      # Install 0mq (for IPython) and affiliated bindings?
INST_ROCKSTAR=0 # Install the Rockstar halo finder?
INST_SCIPY=0    # Install scipy?

# If you need to pass anything to matplotlib, do so here.
MPL_SUPP_LDFLAGS=""
MPL_SUPP_CFLAGS=""
MPL_SUPP_CXXFLAGS=""

# If you want to spawn multiple Make jobs, here's the place to set the
# arguments.  For instance, "-j4"
MAKE_PROCS=""

# Make sure we are NOT being run as root
if [[ $EUID -eq 0 ]]
then
   echo "******************************************************"
   echo "*                                                    *"
   echo "*                                                    *"
   echo "*  IT IS A BAD IDEA TO RUN THIS SCRIPT AS ROOT!!!!   *"
   echo "*                                                    *"
   echo "*                                                    *"
   echo "******************************************************"
   echo
   echo "If you really want to do this, you must manually edit"
   echo "the script to re-enable root-level installation.  Sorry!"
   exit 1
fi
if [[ ${DEST_DIR%/} == /usr/local ]]
then
   echo "******************************************************"
   echo "*                                                    *"
   echo "*                                                    *"
   echo "*  THIS SCRIPT WILL NOT INSTALL TO /usr/local !!!!   *"
   echo "*                                                    *"
   echo "*                                                    *"
   echo "******************************************************"
   exit 1
fi

#------------------------------------------------------------------------------#
#                                                                              #
# Okay, the script starts here.  Feel free to play with it, but hopefully      #
# it'll work as is.                                                            #
#                                                                              #
#------------------------------------------------------------------------------#

LOG_FILE="${DEST_DIR}/yt_install.log"

function write_config
{
    CONFIG_FILE=${DEST_DIR}/.yt_config

    echo INST_HG=${INST_HG} > ${CONFIG_FILE}
    echo INST_ZLIB=${INST_ZLIB} >> ${CONFIG_FILE}
    echo INST_BZLIB=${INST_BZLIB} >> ${CONFIG_FILE}
    echo INST_PNG=${INST_PNG} >> ${CONFIG_FILE}
    echo INST_FTYPE=${INST_FTYPE} >> ${CONFIG_FILE}
    echo INST_SQLITE3=${INST_SQLITE3} >> ${CONFIG_FILE}
    echo INST_PYX=${INST_PYX} >> ${CONFIG_FILE}
    echo INST_0MQ=${INST_0MQ} >> ${CONFIG_FILE}
    echo INST_PY3=${INST_PY3} >> ${CONFIG_FILE}
    echo INST_ROCKSTAR=${INST_ROCKSTAR} >> ${CONFIG_FILE}
    echo INST_SCIPY=${INST_SCIPY} >> ${CONFIG_FILE}
    echo YT_DIR=${YT_DIR} >> ${CONFIG_FILE}
    echo MPL_SUPP_LDFLAGS=${MPL_SUPP_LDFLAGS} >> ${CONFIG_FILE}
    echo MPL_SUPP_CFLAGS=${MPL_SUPP_CFLAGS} >> ${CONFIG_FILE}
    echo MPL_SUPP_CXXFLAGS=${MPL_SUPP_CXXFLAGS} >> ${CONFIG_FILE}
    echo MAKE_PROCS=${MAKE_PROCS} >> ${CONFIG_FILE}
    if [ ${HDF5_DIR} ]
    then
        echo ${HDF5_DIR} >> ${CONFIG_FILE}
    fi
    if [ ${NUMPY_ARGS} ]
    then
        echo ${NUMPY_ARGS} >> ${CONFIG_FILE}
    fi
}

# Write config settings to file.
CONFIG_FILE=${DEST_DIR}/.yt_config
mkdir -p ${DEST_DIR}
if [ -z ${REINST_YT} ] || [ ${REINST_YT} -neq 1 ]
then
    write_config
elif [ ${REINST_YT} ] && [ ${REINST_YT} -eq 1 ] && [ -f ${CONFIG_FILE} ]
then
    USED_CONFIG=1
    source ${CONFIG_FILE}
fi

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
    MYOS=`uname -s`       # A guess at the OS
    if [ "${MYHOST##kraken}" != "${MYHOST}" ]
    then
        echo "Looks like you're on Kraken."
        echo
        echo " ******************************************"
        echo " * It may be better to use the yt module! *"
        echo " *                                        *"
        echo " *   $ module load yt                     *"
        echo " *                                        *"
        echo " ******************************************"
        echo
        echo "IF YOU CHOOSE TO PROCEED:"
        echo "YOU MUST BE IN THE GNU PROGRAMMING ENVIRONMENT"
        echo "   $ module swap PrgEnv-pgi PrgEnv-gnu"
        echo
        return
    fi
    if [ "${MYHOST##nautilus}" != "${MYHOST}" ]
    then
        echo "Looks like you're on Nautilus."
        echo
        echo " ******************************************"
        echo " * It may be better to use the yt module! *"
        echo " *                                        *"
        echo " *   $ module load yt                     *"
        echo " *                                        *"
        echo " ******************************************"
        echo
        echo "NOTE: YOU MUST BE IN THE GNU PROGRAMMING ENVIRONMENT"
        echo "   $ module swap PE-intel PE-gnu"
        echo
        echo "Additionally, note that by default, yt will OVERWRITE"
        echo "any existing installations from Kraken!  You might want"
        echo "to adjust the variable DEST_SUFFIX in the install script."
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
    if [ "${MYHOST##midway}" != "${MYHOST}" ]
    then
        echo "Looks like you're on Midway."
        echo
        echo " ******************************************"
        echo " * It may be better to use the yt module! *"
        echo " *                                        *"
        echo " *   $ module load yt                     *"
        echo " *                                        *"
        echo " ******************************************"
        echo
        return
    fi
    if [ "${MYOS##Darwin}" != "${MYOS}" ]
    then
        echo "Looks like you're running on Mac OSX."
        echo
        echo "NOTE: you must have the Xcode command line tools installed."
        echo
	echo "The instructions for obtaining these tools varies according"
	echo "to your exact OS version.  On older versions of OS X, you"
	echo "must register for an account on the apple developer tools"
	echo "website: https://developer.apple.com/downloads to obtain the"
	echo "download link."
	echo
	echo "We have gathered some additional instructions for each"
	echo "version of OS X below. If you have trouble installing yt"
	echo "after following these instructions, don't hesitate to contact"
	echo "the yt user's e-mail list."
	echo
	echo "You can see which version of OSX you are running by clicking"
	echo "'About This Mac' in the apple menu on the left hand side of"
	echo "menu bar.  We're assuming that you've installed all operating"
	echo "system updates; if you have an older version, we suggest"
	echo "running software update and installing all available updates."
	echo
        echo "OS X 10.5.8: search for and download Xcode 3.1.4 from the"
	echo "Apple developer tools website."
        echo
        echo "OS X 10.6.8: search for and download Xcode 3.2 from the Apple"
	echo "developer tools website.  You can either download the"
	echo "Xcode 3.2.2 Developer Tools package (744 MB) and then use"
	echo "Software Update to update to XCode 3.2.6 or"
	echo "alternatively, you can download the Xcode 3.2.6/iOS SDK"
	echo "bundle (4.1 GB)."
        echo
        echo "OS X 10.7.5: download Xcode 4.2 from the mac app store"
	echo "(search for Xcode)."
        echo "Alternatively, download the Xcode command line tools from"
        echo "the Apple developer tools website."
        echo
	echo "OS X 10.8.4, 10.9, and 10.10: download the appropriate version of"
	echo "Xcode from the mac app store (search for Xcode)."
    echo
	echo "Additionally, you will have to manually install the Xcode"
	echo "command line tools."
    echo
    echo "For OS X 10.8, see:"
   	echo "http://stackoverflow.com/questions/9353444"
	echo
    echo "For OS X 10.9 and 10.10, the command line tools can be installed"
    echo "with the following command:"
    echo "    xcode-select --install"
    echo
    OSX_VERSION=`sw_vers -productVersion`
    if [ "${OSX_VERSION##10.8}" != "${OSX_VERSION}" ]
        then
            MPL_SUPP_CFLAGS="${MPL_SUPP_CFLAGS} -mmacosx-version-min=10.7"
            MPL_SUPP_CXXFLAGS="${MPL_SUPP_CXXFLAGS} -mmacosx-version-min=10.7"
        fi
    fi
    if [ -f /etc/redhat-release ]
    then
        echo "Looks like you're on an Redhat-compatible machine."
        echo
        echo "You need to have these packages installed:"
        echo
        echo "  * openssl-devel"
        echo "  * uuid-devel"
        echo "  * readline-devel"
        echo "  * ncurses-devel"
        echo "  * zip"
        echo "  * gcc-{,c++,gfortran}"
        echo "  * make"
        echo "  * patch"
        echo 
        echo "You can accomplish this by executing:"
        echo "$ sudo yum install gcc gcc-c++ gcc-gfortran make patch zip"
        echo "$ sudo yum install ncurses-devel uuid-devel openssl-devel readline-devel"
    fi
    if [ -f /etc/SuSE-release ] && [ `grep --count SUSE /etc/SuSE-release` -gt 0 ]
    then
        echo "Looks like you're on an OpenSUSE-compatible machine."
        echo
        echo "You need to have these packages installed:"
        echo
        echo "  * devel_C_C++"
        echo "  * libopenssl-devel"
        echo "  * libuuid-devel"
        echo "  * zip"
        echo "  * gcc-c++"
        echo
        echo "You can accomplish this by executing:"
        echo
        echo "$ sudo zypper install -t pattern devel_C_C++"
        echo "$ sudo zypper install gcc-c++ libopenssl-devel libuuid-devel zip"
        echo
        echo "I am also setting special configure arguments to Python to"
        echo "specify control lib/lib64 issues."
        PYCONF_ARGS="--libdir=${DEST_DIR}/lib"
    fi
    if [ -f /etc/lsb-release ] && [ `grep --count buntu /etc/lsb-release` -gt 0 ]
    then
        echo "Looks like you're on an Ubuntu-compatible machine."
        echo
        echo "You need to have these packages installed:"
        echo
        echo "  * libssl-dev"
        echo "  * build-essential"
        echo "  * libncurses5"
        echo "  * libncurses5-dev"
        echo "  * zip"
        echo "  * uuid-dev"
        echo "  * libfreetype6-dev"
        echo "  * tk-dev"
        echo
        echo "You can accomplish this by executing:"
        echo
        echo "$ sudo apt-get install libssl-dev build-essential libncurses5 libncurses5-dev zip uuid-dev libfreetype6-dev tk-dev"
        echo
        echo
        echo " Additionally, if you want to put yt's lib dir in your LD_LIBRARY_PATH"
        echo " so you can use yt without the activate script, you might "
        echo " want to consider turning off LIBZ and FREETYPE in this"
        echo " install script by editing this file and setting"
        echo
        echo " INST_ZLIB=0"
        echo " INST_FTYPE=0"
        echo
        echo " to avoid conflicts with other command-line programs "
        echo " (like eog and evince, for example)."
    fi
    if [ $INST_SCIPY -eq 1 ]
    then
	echo
	echo "Looks like you've requested that the install script build SciPy."
	echo
	echo "If the SciPy build fails, please uncomment one of the the lines"
	echo "at the top of the install script that sets NUMPY_ARGS, delete"
	echo "any broken installation tree, and re-run the install script"
	echo "verbatim."
	echo
	echo "If that doesn't work, don't hesitate to ask for help on the yt"
	echo "user's mailing list."
	echo
    fi
    if [ ! -z "${CFLAGS}" ]
    then
        echo "******************************************"
        echo "******************************************"
        echo "**                                      **"
        echo "**    Your CFLAGS is not empty.         **"
        echo "**    This can break h5py compilation.  **"
        echo "**                                      **"
        echo "******************************************"
        echo "******************************************"
    fi
}


echo
echo
echo "========================================================================"
echo
echo "Hi there!  This is the yt installation script.  We're going to download"
echo "some stuff and install it to create a self-contained, isolated"
echo "environment for yt to run within."
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

printf "%-15s = %s so I " "INST_PNG" "${INST_PNG}"
get_willwont ${INST_PNG}
echo "be installing libpng"

printf "%-15s = %s so I " "INST_FTYPE" "${INST_FTYPE}"
get_willwont ${INST_FTYPE}
echo "be installing freetype2"

printf "%-15s = %s so I " "INST_SQLITE3" "${INST_SQLITE3}"
get_willwont ${INST_SQLITE3}
echo "be installing SQLite3"

printf "%-15s = %s so I " "INST_PY3" "${INST_PY3}"
get_willwont ${INST_PY3}
echo "be installing Python 3"

printf "%-15s = %s so I " "INST_HG" "${INST_HG}"
get_willwont ${INST_HG}
echo "be installing Mercurial"

printf "%-15s = %s so I " "INST_PYX" "${INST_PYX}"
get_willwont ${INST_PYX}
echo "be installing PyX"

printf "%-15s = %s so I " "INST_SCIPY" "${INST_SCIPY}"
get_willwont ${INST_SCIPY}
echo "be installing scipy"

printf "%-15s = %s so I " "INST_0MQ" "${INST_0MQ}"
get_willwont ${INST_0MQ}
echo "be installing ZeroMQ"

printf "%-15s = %s so I " "INST_ROCKSTAR" "${INST_ROCKSTAR}"
get_willwont ${INST_ROCKSTAR}
echo "be installing Rockstar"

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
if [ ${USED_CONFIG} ]
then
    echo "Settings were loaded from ${CONFIG_FILE}."
    echo "Remove this file if you wish to return to the default settings."
    echo
fi
echo "========================================================================"
echo
read -p "[hit enter] "
echo
echo "Awesome!  Here we go."
echo

function do_exit
{
    echo "********************************************"
    echo "        FAILURE REPORT:"
    echo "********************************************"
    echo
    tail -n 10 ${LOG_FILE}
    echo
    echo "********************************************"
    echo "********************************************"
    echo "Failure.  Check ${LOG_FILE}.  The last 10 lines are above."
    exit 1
}

if [ $INST_PY3 -eq 1 ]
then
	 PYTHON_EXEC='python3.4'
else 
	 PYTHON_EXEC='python2.7'
fi

function do_setup_py
{
    [ -e $1/done ] && return
    LIB=$1
    shift
    if [ -z "$@" ]
    then
        echo "Installing $LIB"
    else
        echo "Installing $LIB (arguments: '$@')"
    fi
    [ ! -e $LIB/extracted ] && tar xfz $LIB.tar.gz
    touch $LIB/extracted
    BUILD_ARGS=""
    if [[ $LIB =~ .*mercurial.* ]] 
    then
        PYEXE="python2.7"
    else
        PYEXE=${PYTHON_EXEC}
    fi
    case $LIB in
        *h5py*)
            pushd $LIB &> /dev/null
            ( ${DEST_DIR}/bin/${PYTHON_EXEC} setup.py configure --hdf5=${HDF5_DIR} 2>&1 ) 1>> ${LOG_FILE} || do_exit
            popd &> /dev/null
            ;;
        *numpy*)
            if [ -e ${DEST_DIR}/lib/${PYTHON_EXEC}/site-packages/numpy/__init__.py ]
            then
                VER=$(${DEST_DIR}/bin/${PYTHON_EXEC} -c 'from distutils.version import StrictVersion as SV; \
                                                 import numpy; print SV(numpy.__version__) < SV("1.8.0")')
                if [ $VER == "True" ]
                then
                    echo "Removing previous NumPy instance (see issue #889)"
                    rm -rf ${DEST_DIR}/lib/${PYTHON_EXEC}/site-packages/{numpy*,*.pth}
                fi
            fi
            ;;
        *)
            ;;
    esac
    cd $LIB
    ( ${DEST_DIR}/bin/${PYEXE} setup.py build ${BUILD_ARGS} $* 2>&1 ) 1>> ${LOG_FILE} || do_exit
    ( ${DEST_DIR}/bin/${PYEXE} setup.py install    2>&1 ) 1>> ${LOG_FILE} || do_exit
    touch done
    cd ..
}

if type -P wget &>/dev/null
then
    echo "Using wget"
    export GETFILE="wget -nv"
else
    echo "Using curl"
    export GETFILE="curl -sSO"
fi

if type -P sha512sum &> /dev/null
then
    echo "Using sha512sum"
    export SHASUM="sha512sum"
elif type -P shasum &> /dev/null
then
    echo "Using shasum -a 512"
    export SHASUM="shasum -a 512"
else
    echo
    echo "I am unable to locate any shasum-like utility."
    echo "ALL FILE INTEGRITY IS NOT VERIFIABLE."
    echo "THIS IS PROBABLY A BIG DEAL."
    echo
    echo "(I'll hang out for a minute for you to consider this.)"
    sleep 60
fi

function get_ytproject
{
    [ -e $1 ] && return
    echo "Downloading $1 from yt-project.org"
    ${GETFILE} "http://yt-project.org/dependencies/$1" || do_exit
    ( ${SHASUM} -c $1.sha512 2>&1 ) 1>> ${LOG_FILE} || do_exit
}

function get_ytdata
{
    echo "Downloading $1 from yt-project.org"
    [ -e $1 ] && return
    ${GETFILE} "http://yt-project.org/data/$1" || do_exit
    ( ${SHASUM} -c $1.sha512 2>&1 ) 1>> ${LOG_FILE} || do_exit
}

ORIG_PWD=`pwd`

if [ -z "${DEST_DIR}" ]
then
    echo "Edit this script, set the DEST_DIR parameter and re-run."
    exit 1
fi

# Get supplemental data.

mkdir -p ${DEST_DIR}/data
cd ${DEST_DIR}/data
echo 'de6d8c6ea849f0206d219303329a0276b3cce7c051eec34377d42aacbe0a4f47ac5145eb08966a338ecddd2b83c8f787ca9956508ad5c39ee2088ad875166410  cloudy_emissivity.h5' > cloudy_emissivity.h5.sha512
[ ! -e cloudy_emissivity.h5 ] && get_ytdata cloudy_emissivity.h5
echo '0f714ae2eace0141b1381abf1160dc8f8a521335e886f99919caf3beb31df1fe271d67c7b2a804b1467949eb16b0ef87a3d53abad0e8160fccac1e90d8d9e85f  apec_emissivity.h5' > apec_emissivity.h5.sha512
[ ! -e apec_emissivity.h5 ] && get_ytdata apec_emissivity.h5

# Set paths to what they should be when yt is activated.
export PATH=${DEST_DIR}/bin:$PATH
export LD_LIBRARY_PATH=${DEST_DIR}/lib:$LD_LIBRARY_PATH
export PYTHONPATH=${DEST_DIR}/lib/${PYTHON_EXEC}/site-packages

mkdir -p ${DEST_DIR}/src
cd ${DEST_DIR}/src

PYTHON2='Python-2.7.9'
PYTHON3='Python-3.4.3'
CYTHON='Cython-0.22'
PYX='PyX-0.12.1'
BZLIB='bzip2-1.0.6'
FREETYPE_VER='freetype-2.4.12' 
H5PY='h5py-2.5.0'
HDF5='hdf5-1.8.14' 
IPYTHON='ipython-2.4.1'
LAPACK='lapack-3.4.2'
PNG=libpng-1.6.3
MATPLOTLIB='matplotlib-1.4.3'
MERCURIAL='mercurial-3.4'
NOSE='nose-1.3.6'
NUMPY='numpy-1.9.2'
PYTHON_HGLIB='python-hglib-1.6'
PYZMQ='pyzmq-14.5.0'
ROCKSTAR='rockstar-0.99.6'
SCIPY='scipy-0.15.1'
SQLITE='sqlite-autoconf-3071700'
SYMPY='sympy-0.7.6'
TORNADO='tornado-4.0.2'
ZEROMQ='zeromq-4.0.5'
ZLIB='zlib-1.2.8'
SETUPTOOLS='setuptools-18.0.1'

# Now we dump all our SHA512 files out.
echo '856220fa579e272ac38dcef091760f527431ff3b98df9af6e68416fcf77d9659ac5abe5c7dee41331f359614637a4ff452033085335ee499830ed126ab584267  Cython-0.22.tar.gz' > Cython-0.22.tar.gz.sha512
echo '4941f5aa21aff3743546495fb073c10d2657ff42b2aff401903498638093d0e31e344cce778980f28a7170c6d29eab72ac074277b9d4088376e8692dc71e55c1  PyX-0.12.1.tar.gz' > PyX-0.12.1.tar.gz.sha512
echo 'a42f28ed8e49f04cf89e2ea7434c5ecbc264e7188dcb79ab97f745adf664dd9ab57f9a913543731635f90859536244ac37dca9adf0fc2aa1b215ba884839d160  Python-2.7.9.tgz' > Python-2.7.9.tgz.sha512
echo '609cc82586fabecb25f25ecb410f2938e01d21cde85dd3f8824fe55c6edde9ecf3b7609195473d3fa05a16b9b121464f5414db1a0187103b78ea6edfa71684a7  Python-3.4.3.tgz' > Python-3.4.3.tgz.sha512
echo '276bd9c061ec9a27d478b33078a86f93164ee2da72210e12e2c9da71dcffeb64767e4460b93f257302b09328eda8655e93c4b9ae85e74472869afbeae35ca71e  blas.tar.gz' > blas.tar.gz.sha512
echo '00ace5438cfa0c577e5f578d8a808613187eff5217c35164ffe044fbafdfec9e98f4192c02a7d67e01e5a5ccced630583ad1003c37697219b0f147343a3fdd12  bzip2-1.0.6.tar.gz' > bzip2-1.0.6.tar.gz.sha512
echo 'a296dfcaef7e853e58eed4e24b37c4fa29cfc6ac688def048480f4bb384b9e37ca447faf96eec7b378fd764ba291713f03ac464581d62275e28eb2ec99110ab6  reason-js-20120623.zip' > reason-js-20120623.zip.sha512
echo '609a68a3675087e0cc95268574f31e104549daa48efe15a25a33b8e269a93b4bd160f4c3e8178dca9c950ef5ca514b039d6fd1b45db6af57f25342464d0429ce  freetype-2.4.12.tar.gz' > freetype-2.4.12.tar.gz.sha512
echo '4a83f9ae1855a7fad90133b327d426201c8ccfd2e7fbe9f39b2d61a2eee2f3ebe2ea02cf80f3d4e1ad659f8e790c173df8cc99b87d0b7ce63d34aa88cfdc7939  h5py-2.5.0.tar.gz' > h5py-2.5.0.tar.gz.sha512
echo '4073fba510ccadaba41db0939f909613c9cb52ba8fb6c1062fc9118edc601394c75e102310be1af4077d07c9b327e6bbb1a6359939a7268dc140382d0c1e0199  hdf5-1.8.14.tar.gz' > hdf5-1.8.14.tar.gz.sha512
echo 'a9cffc08ba10c47b0371b05664e55eee0562a30ef0d4bbafae79e52e5b9727906c45840c0918122c06c5672ac65e6eb381399f103e1a836aca003eda81b2acde  ipython-2.4.1.tar.gz' > ipython-2.4.1.tar.gz.sha512
echo '8770214491e31f0a7a3efaade90eee7b0eb20a8a6ab635c5f854d78263f59a1849133c14ef5123d01023f0110cbb9fc6f818da053c01277914ae81473430a952  lapack-3.4.2.tar.gz' > lapack-3.4.2.tar.gz.sha512
echo '887582e5a22e4cde338aa8fec7a89f6dd31f2f02b8842735f00f970f64582333fa03401cea6d01704083403c7e8b7ebc26655468ce930165673b33efa4bcd586  libpng-1.6.3.tar.gz' > libpng-1.6.3.tar.gz.sha512
echo '51b0f58b2618b47b653e17e4f6b6a1215d3a3b0f1331ce3555cc7435e365d9c75693f289ce12fe3bf8f69fd57b663e545f0f1c2c94e81eaa661cac0689e125f5  matplotlib-1.4.3.tar.gz' > matplotlib-1.4.3.tar.gz.sha512
echo 'a61b0d4cf528136991243bb23ac972c11c50ab5681d09f8b2d12cf7d37d3a9d76262f7fe6e7a1834bf6d03e8dc0ebbd9231da982e049e09830341dabefe5d064  mercurial-3.4.tar.gz' > mercurial-3.4.tar.gz.sha512
echo 'd0cede08dc33a8ac0af0f18063e57f31b615f06e911edb5ca264575174d8f4adb4338448968c403811d9dcc60f38ade3164662d6c7b69b499f56f0984bb6283c  nose-1.3.6.tar.gz' > nose-1.3.6.tar.gz.sha512
echo '70470ebb9afef5dfd0c83ceb7a9d5f1b7a072b1a9b54b04f04f5ed50fbaedd5b4906bd500472268d478f94df9e749a88698b1ff30f2d80258e7f3fec040617d9  numpy-1.9.2.tar.gz' > numpy-1.9.2.tar.gz.sha512
echo 'bfd10455e74e30df568c4c4827140fb6cc29893b0e062ce1764bd52852ec7487a70a0f5ea53c3fca7886f5d36365c9f4db52b8c93cad35fb67beeb44a2d56f2d  python-hglib-1.6.tar.gz' > python-hglib-1.6.tar.gz.sha512
echo '20164f7b05c308e0f089c07fc46b1c522094f3ac136f2e0bba84f19cb63dfd36152a2465df723dd4d93c6fbd2de4f0d94c160e2bbc353a92cfd680eb03cbdc87  pyzmq-14.5.0.tar.gz' > pyzmq-14.5.0.tar.gz.sha512
echo 'fff4412d850c431a1b4e6ee3b17958ee5ab3beb81e6cb8a8e7d56d368751eaa8781d7c3e69d932dc002d718fddc66a72098acfe74cfe29ec80b24e6736317275  scipy-0.15.1.tar.gz' > scipy-0.15.1.tar.gz.sha512
echo '96f3e51b46741450bc6b63779c10ebb4a7066860fe544385d64d1eda52592e376a589ef282ace2e1df73df61c10eab1a0d793abbdaf770e60289494d4bf3bcb4  sqlite-autoconf-3071700.tar.gz' > sqlite-autoconf-3071700.tar.gz.sha512
echo 'ce0f1a17ac01eb48aec31fc0ad431d9d7ed9907f0e8584a6d79d0ffe6864fe62e203fe3f2a3c3e4e3d485809750ce07507a6488e776a388a7a9a713110882fcf  sympy-0.7.6.tar.gz' > sympy-0.7.6.tar.gz.sha512
echo '93591068dc63af8d50a7925d528bc0cccdd705232c529b6162619fe28dddaf115e8a460b1842877d35160bd7ed480c1bd0bdbec57d1f359085bd1814e0c1c242  tornado-4.0.2.tar.gz' > tornado-4.0.2.tar.gz.sha512
echo '0d928ed688ed940d460fa8f8d574a9819dccc4e030d735a8c7db71b59287ee50fa741a08249e356c78356b03c2174f2f2699f05aa7dc3d380ed47d8d7bab5408  zeromq-4.0.5.tar.gz' > zeromq-4.0.5.tar.gz.sha512
echo 'ece209d4c7ec0cb58ede791444dc754e0d10811cbbdebe3df61c0fd9f9f9867c1c3ccd5f1827f847c005e24eef34fb5bf87b5d3f894d75da04f1797538290e4a  zlib-1.2.8.tar.gz' > zlib-1.2.8.tar.gz.sha512
echo '9b318ce2ee2cf787929dcb886d76c492b433e71024fda9452d8b4927652a298d6bd1bdb7a4c73883a98e100024f89b46ea8aa14b250f896e549e6dd7e10a6b41  setuptools-18.0.1.tar.gz' > setuptools-18.0.1.tar.gz.sha512
# Individual processes
[ -z "$HDF5_DIR" ] && get_ytproject $HDF5.tar.gz
[ $INST_ZLIB -eq 1 ] && get_ytproject $ZLIB.tar.gz
[ $INST_BZLIB -eq 1 ] && get_ytproject $BZLIB.tar.gz
[ $INST_PNG -eq 1 ] && get_ytproject $PNG.tar.gz
[ $INST_FTYPE -eq 1 ] && get_ytproject $FREETYPE_VER.tar.gz
[ $INST_SQLITE3 -eq 1 ] && get_ytproject $SQLITE.tar.gz
[ $INST_PYX -eq 1 ] && get_ytproject $PYX.tar.gz
[ $INST_0MQ -eq 1 ] && get_ytproject $ZEROMQ.tar.gz
[ $INST_0MQ -eq 1 ] && get_ytproject $PYZMQ.tar.gz
[ $INST_0MQ -eq 1 ] && get_ytproject $TORNADO.tar.gz
[ $INST_SCIPY -eq 1 ] && get_ytproject $SCIPY.tar.gz
[ $INST_SCIPY -eq 1 ] && get_ytproject blas.tar.gz
[ $INST_SCIPY -eq 1 ] && get_ytproject $LAPACK.tar.gz
[ $INST_HG -eq 1 ] && get_ytproject $MERCURIAL.tar.gz
[ $INST_PY3 -eq 1 ] && get_ytproject $PYTHON3.tgz
get_ytproject $PYTHON2.tgz
get_ytproject $NUMPY.tar.gz
get_ytproject $MATPLOTLIB.tar.gz
get_ytproject $IPYTHON.tar.gz
get_ytproject $H5PY.tar.gz
get_ytproject $CYTHON.tar.gz
get_ytproject reason-js-20120623.zip
get_ytproject $NOSE.tar.gz
get_ytproject $PYTHON_HGLIB.tar.gz
get_ytproject $SYMPY.tar.gz
get_ytproject $SETUPTOOLS.tar.gz
if [ $INST_BZLIB -eq 1 ]
then
    if [ ! -e $BZLIB/done ]
    then
        [ ! -e $BZLIB ] && tar xfz $BZLIB.tar.gz
        echo "Installing BZLIB"
        cd $BZLIB
        if [ `uname` = "Darwin" ]
        then
            if [ -z "${CC}" ]
            then
                sed -i.bak 's/soname/install_name/' Makefile-libbz2_so
            else
                sed -i.bak -e 's/soname/install_name/' -e "s|CC=gcc|CC=${CC}|" Makefile-libbz2_so
            fi
        fi
        ( make install CFLAGS=-fPIC LDFLAGS=-fPIC PREFIX=${DEST_DIR} 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make -f Makefile-libbz2_so CFLAGS=-fPIC LDFLAGS=-fPIC PREFIX=${DEST_DIR} 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( cp -v libbz2.so.1.0.6 ${DEST_DIR}/lib 2>&1 ) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
    BZLIB_DIR=${DEST_DIR}
    export LDFLAGS="${LDFLAGS} -L${BZLIB_DIR}/lib/ -L${BZLIB_DIR}/lib64/"
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${BZLIB_DIR}/lib/"
fi

if [ $INST_ZLIB -eq 1 ]
then
    if [ ! -e $ZLIB/done ]
    then
        [ ! -e $ZLIB ] && tar xfz $ZLIB.tar.gz
        echo "Installing ZLIB"
        cd $ZLIB
        ( ./configure --shared --prefix=${DEST_DIR}/ 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make clean 2>&1) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
    ZLIB_DIR=${DEST_DIR}
    export LDFLAGS="${LDFLAGS} -L${ZLIB_DIR}/lib/ -L${ZLIB_DIR}/lib64/"
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${ZLIB_DIR}/lib/"
fi

if [ $INST_PNG -eq 1 ]
then
    if [ ! -e $PNG/done ]
    then
        [ ! -e $PNG ] && tar xfz $PNG.tar.gz
        echo "Installing PNG"
        cd $PNG
        ( ./configure CPPFLAGS=-I${DEST_DIR}/include CFLAGS=-I${DEST_DIR}/include --prefix=${DEST_DIR}/ 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make clean 2>&1) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
    PNG_DIR=${DEST_DIR}
    export LDFLAGS="${LDFLAGS} -L${PNG_DIR}/lib/ -L${PNG_DIR}/lib64/"
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${PNG_DIR}/lib/"
fi

if [ $INST_FTYPE -eq 1 ]
then
    if [ ! -e $FREETYPE_VER/done ]
    then
        [ ! -e $FREETYPE_VER ] && tar xfz $FREETYPE_VER.tar.gz
        echo "Installing FreeType2"
        cd $FREETYPE_VER
        ( ./configure CFLAGS=-I${DEST_DIR}/include --prefix=${DEST_DIR}/ 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make clean 2>&1) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
    FTYPE_DIR=${DEST_DIR}
    export LDFLAGS="${LDFLAGS} -L${FTYPE_DIR}/lib/ -L${FTYPE_DIR}/lib64/"
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${FTYPE_DIR}/lib/"
fi

if [ -z "$HDF5_DIR" ]
then
    if [ ! -e $HDF5/done ]
    then
        [ ! -e $HDF5 ] && tar xfz $HDF5.tar.gz
        echo "Installing HDF5"
        cd $HDF5
        ( ./configure --prefix=${DEST_DIR}/ --enable-shared 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make ${MAKE_PROCS} install 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make clean 2>&1) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
    export HDF5_DIR=${DEST_DIR}
else
    export HDF5_DIR=${HDF5_DIR}
fi
export HDF5_API=16

if [ $INST_SQLITE3 -eq 1 ]
then
    if [ ! -e $SQLITE/done ]
    then
        [ ! -e $SQLITE ] && tar xfz $SQLITE.tar.gz
        echo "Installing SQLite3"
        cd $SQLITE
        ( ./configure --prefix=${DEST_DIR}/ 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make ${MAKE_PROCS} install 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make clean 2>&1) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
fi

if [ ! -e $PYTHON2/done ]
then
    echo "Installing Python 2. This may take a while, but don't worry. yt loves you."
    [ ! -e $PYTHON2 ] && tar xfz $PYTHON2.tgz
    cd $PYTHON2
    ( ./configure --prefix=${DEST_DIR}/ ${PYCONF_ARGS} 2>&1 ) 1>> ${LOG_FILE} || do_exit

    ( make ${MAKE_PROCS} 2>&1 ) 1>> ${LOG_FILE} || do_exit
    ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
    ( ln -sf ${DEST_DIR}/bin/python2.7 ${DEST_DIR}/bin/pyyt 2>&1 ) 1>> ${LOG_FILE}
    ( make clean 2>&1) 1>> ${LOG_FILE} || do_exit
    touch done
    cd ..
fi

if [ $INST_PY3 -eq 1 ]
then
    if [ ! -e $PYTHON3/done ]
    then
        echo "Installing Python 3. Because two Pythons are better than one."
        [ ! -e $PYTHON3 ] && tar xfz $PYTHON3.tgz
        cd $PYTHON3
        ( ./configure --prefix=${DEST_DIR}/ ${PYCONF_ARGS} 2>&1 ) 1>> ${LOG_FILE} || do_exit

        ( make ${MAKE_PROCS} 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( ln -sf ${DEST_DIR}/bin/python3.4 ${DEST_DIR}/bin/pyyt 2>&1 ) 1>> ${LOG_FILE}
        ( ln -sf ${DEST_DIR}/bin/python3.4 ${DEST_DIR}/bin/python 2>&1 ) 1>> ${LOG_FILE}
        ( ln -sf ${DEST_DIR}/bin/python3-config ${DEST_DIR}/bin/python-config 2>&1 ) 1>> ${LOG_FILE}
        ( make clean 2>&1) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
fi

export PYTHONPATH=${DEST_DIR}/lib/${PYTHON_EXEC}/site-packages/

# Install setuptools
do_setup_py $SETUPTOOLS

if [ $INST_HG -eq 1 ]
then
    do_setup_py $MERCURIAL
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
        echo "Cloning yt"
        YT_DIR="$PWD/yt-hg/"
        ( ${HG_EXEC} --debug clone https://bitbucket.org/yt_analysis/yt-supplemental/ 2>&1 ) 1>> ${LOG_FILE}
        # Recently the hg server has had some issues with timeouts.  In lieu of
        # a new webserver, we are now moving to a three-stage process.
        # First we clone the repo, but only up to r0.
        ( ${HG_EXEC} --debug clone https://bitbucket.org/yt_analysis/yt/ ./yt-hg 2>&1 ) 1>> ${LOG_FILE}
        # Now we update to the branch we're interested in.
        ( ${HG_EXEC} -R ${YT_DIR} up -C ${BRANCH} 2>&1 ) 1>> ${LOG_FILE}
    elif [ -e yt-hg ]
    then
        YT_DIR="$PWD/yt-hg/"
    fi
    echo Setting YT_DIR=${YT_DIR}
fi

# This fixes problems with gfortran linking.
unset LDFLAGS
 
echo "Installing pip"
( ${GETFILE} https://bootstrap.pypa.io/get-pip.py 2>&1 ) 1>> ${LOG_FILE} || do_exit
( ${DEST_DIR}/bin/${PYTHON_EXEC} get-pip.py 2>&1 ) 1>> ${LOG_FILE} || do_exit

if [ $INST_SCIPY -eq 0 ]
then
    do_setup_py $NUMPY ${NUMPY_ARGS}
else
    if [ ! -e $SCIPY/done ]
    then
	if [ ! -e BLAS/done ]
	then
	    tar xfz blas.tar.gz
	    echo "Building BLAS"
	    cd BLAS
	    gfortran -O2 -fPIC -fno-second-underscore -c *.f
	    ( ar r libfblas.a *.o 2>&1 ) 1>> ${LOG_FILE}
	    ( ranlib libfblas.a 2>&1 ) 1>> ${LOG_FILE}
	    rm -rf *.o
	    touch done
	    cd ..
	fi
	if [ ! -e $LAPACK/done ]
	then
	    tar xfz $LAPACK.tar.gz
	    echo "Building LAPACK"
	    cd $LAPACK/
	    cp INSTALL/make.inc.gfortran make.inc
	    ( make lapacklib OPTS="-fPIC -O2" NOOPT="-fPIC -O0" CFLAGS=-fPIC LDFLAGS=-fPIC 2>&1 ) 1>> ${LOG_FILE} || do_exit
	    touch done
	    cd ..
	fi
    fi
    export BLAS=$PWD/BLAS/libfblas.a
    export LAPACK=$PWD/$LAPACK/liblapack.a
    do_setup_py $NUMPY ${NUMPY_ARGS}
    do_setup_py $SCIPY ${NUMPY_ARGS}
fi

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
if [ -n "${MPL_SUPP_CFLAGS}" ]
then
    OLD_CFLAGS=${CFLAGS}
    export CFLAGS="${MPL_SUPP_CFLAGS}"
    echo "Setting CFLAGS ${CFLAGS}"
fi
# Now we set up the basedir for matplotlib:
mkdir -p ${DEST_DIR}/src/$MATPLOTLIB
echo "[directories]" >> ${DEST_DIR}/src/$MATPLOTLIB/setup.cfg
echo "basedirlist = ${DEST_DIR}" >> ${DEST_DIR}/src/$MATPLOTLIB/setup.cfg
if [ `uname` = "Darwin" ]
then
   echo "[gui_support]" >> ${DEST_DIR}/src/$MATPLOTLIB/setup.cfg
   echo "macosx = False" >> ${DEST_DIR}/src/$MATPLOTLIB/setup.cfg
fi

_user_DISPLAY=$DISPLAY
unset DISPLAY   # see (yt-user link missing: "Installation failure" 01/29/15)
do_setup_py $MATPLOTLIB
export DISPLAY=${_user_DISPLAY}
if [ -n "${OLD_LDFLAGS}" ]
then
    export LDFLAG=${OLD_LDFLAGS}
fi
[ -n "${OLD_LDFLAGS}" ] && export LDFLAGS=${OLD_LDFLAGS}
[ -n "${OLD_CXXFLAGS}" ] && export CXXFLAGS=${OLD_CXXFLAGS}
[ -n "${OLD_CFLAGS}" ] && export CFLAGS=${OLD_CFLAGS}

# Now we do our IPython installation, which has two optional dependencies.
if [ $INST_0MQ -eq 1 ]
then
    if [ ! -e $ZEROMQ/done ]
    then
        [ ! -e $ZEROMQ ] && tar xfz $ZEROMQ.tar.gz
        echo "Installing ZeroMQ"
        cd $ZEROMQ
        ( ./configure --prefix=${DEST_DIR}/ 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make clean 2>&1) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
    do_setup_py $PYZMQ --zmq=${DEST_DIR}
    do_setup_py $TORNADO
fi

do_setup_py $IPYTHON
do_setup_py $CYTHON
do_setup_py $H5PY
do_setup_py $NOSE
do_setup_py $PYTHON_HGLIB
do_setup_py $SYMPY
[ $INST_PYX -eq 1 ] && do_setup_py $PYX

( ${DEST_DIR}/bin/pip install jinja2 2>&1 ) 1>> ${LOG_FILE}

# Now we build Rockstar and set its environment variable.
if [ $INST_ROCKSTAR -eq 1 ]
then
    if [ ! -e rockstar/done ]
    then
        echo "Building Rockstar"
        if [ ! -e rockstar ]
        then
            ( hg clone http://bitbucket.org/MatthewTurk/rockstar 2>&1 ) 1>> ${LOG_FILE}
        fi
        cd rockstar
        ( hg pull 2>&1 ) 1>> ${LOG_FILE}
        ( hg up -C tip 2>&1 ) 1>> ${LOG_FILE}
        ( make lib 2>&1 ) 1>> ${LOG_FILE} || do_exit
        cp librockstar.so ${DEST_DIR}/lib
        ROCKSTAR_DIR=${DEST_DIR}/src/rockstar
        echo $ROCKSTAR_DIR > ${YT_DIR}/rockstar.cfg
        touch done
        cd ..
    fi
fi

echo "Doing yt update, wiping local changes and updating to branch ${BRANCH}"
MY_PWD=`pwd`
cd $YT_DIR
( ${HG_EXEC} pull 2>1 && ${HG_EXEC} up -C 2>1 ${BRANCH} 2>&1 ) 1>> ${LOG_FILE}

echo "Installing yt"
[ $INST_PNG -eq 1 ] && echo $PNG_DIR > png.cfg
( export PATH=$DEST_DIR/bin:$PATH ; ${DEST_DIR}/bin/${PYTHON_EXEC} setup.py develop 2>&1 ) 1>> ${LOG_FILE} || do_exit
touch done
cd $MY_PWD

if !( ( ${DEST_DIR}/bin/${PYTHON_EXEC} -c "import readline" 2>&1 )>> ${LOG_FILE}) || \
	[[ "${MYOS##Darwin}" != "${MYOS}" && $INST_PY3 -eq 1 ]] 
then
    if !( ( ${DEST_DIR}/bin/${PYTHON_EXEC} -c "import gnureadline" 2>&1 )>> ${LOG_FILE})
    then
        echo "Installing pure-python readline"
        ( ${DEST_DIR}/bin/pip install gnureadline 2>&1 ) 1>> ${LOG_FILE}
    fi
fi

if [ -e $HOME/.matplotlib/fontList.cache ] && \
   ( grep -q python2.6 $HOME/.matplotlib/fontList.cache )
then
    echo "WARNING WARNING WARNING WARNING WARNING WARNING WARNING"
    echo "*******************************************************"
    echo
    echo "  You likely need to remove your old fontList.cache!"
    echo "  You can do this with this command:"
    echo ""
    echo "  rm $HOME/.matplotlib/fontList.cache"
    echo
    echo "*******************************************************"
fi

# Add the environment scripts
( cp ${YT_DIR}/doc/activate ${DEST_DIR}/bin/activate 2>&1 ) 1>> ${LOG_FILE}
sed -i.bak -e "s,__YT_DIR__,${DEST_DIR}," ${DEST_DIR}/bin/activate
( cp ${YT_DIR}/doc/activate.csh ${DEST_DIR}/bin/activate.csh 2>&1 ) 1>> ${LOG_FILE}
sed -i.bak -e "s,__YT_DIR__,${DEST_DIR}," ${DEST_DIR}/bin/activate.csh

function print_afterword
{
    echo
    echo
    echo "========================================================================"
    echo
    echo "yt is now installed in $DEST_DIR ."
    echo
    echo "To run from this new installation, use the activate script for this "
    echo "environment."
    echo
    echo "    $ source $DEST_DIR/bin/activate"
    echo
    echo "This modifies the environment variables YT_DEST, PATH, PYTHONPATH, and"
    echo "LD_LIBRARY_PATH to match your new yt install.  If you use csh, just"
    echo "append .csh to the above."
    echo
    echo "To get started with yt, check out the orientation:"
    echo
    echo "    http://yt-project.org/doc/quickstart/"
    echo
    echo "The source for yt is located at:"
    echo "    $YT_DIR"
    if [ $INST_HG -eq 1 ]
    then
      echo
      echo "Mercurial has also been installed:"
      echo
      echo "$DEST_DIR/bin/hg"
      echo
    fi
    echo
    echo "For support, see the website and join the mailing list:"
    echo
    echo "    http://yt-project.org/"
    echo "    http://yt-project.org/data/      (Sample data)"
    echo "    http://yt-project.org/doc/       (Docs)"
    echo
    echo "    http://lists.spacepope.org/listinfo.cgi/yt-users-spacepope.org"
    echo
    echo "========================================================================"
    echo
    echo "Oh, look at me, still talking when there's science to do!"
    echo "Good luck, and email the user list if you run into any problems."
}

print_afterword
print_afterword >> ${LOG_FILE}

echo "yt dependencies were last updated on" > ${DEST_DIR}/.yt_update
date >> ${DEST_DIR}/.yt_update
