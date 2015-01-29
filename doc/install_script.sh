#
# Hi there!  Welcome to the yt installation script.
#
# This script is designed to create a fully isolated Python installation
# with the dependencies you need to run yt.
#
# There are a few options, but you only need to set *one* of them.  And
# that's the next one, DEST_DIR.  But, if you want to use an existing HDF5
# installation you can set HDF5_DIR, or if you want to use some other
# subversion checkout of yt, you can set YT_DIR, too.  (It'll already
# check the current directory and one up.
#
# If you experience problems, please visit the Help section at 
# http://yt-project.org.
#

DEST_SUFFIX="yt-`uname -m`"
DEST_DIR="`pwd`/${DEST_SUFFIX/ /}"   # Installation location
BRANCH="yt" # This is the branch to which we will forcibly update.

if [ ${REINST_YT} ] && [ ${REINST_YT} -eq 1 ] && [ -n ${YT_DEST} ]
then
    DEST_DIR=${YT_DEST}
fi

# Here's where you put the HDF5 path if you like; otherwise it'll download it
# and install it on its own
#HDF5_DIR=

# If you need to supply arguments to the NumPy or SciPy build, supply them here
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
INST_FTYPE=1    # Install FreeType2 locally?
INST_ENZO=0     # Clone a copy of Enzo?
INST_SQLITE3=1  # Install a local version of SQLite3?
INST_PYX=0      # Install PyX?  Sometimes PyX can be problematic without a
                # working TeX installation.
INST_0MQ=1      # Install 0mq (for IPython) and affiliated bindings?
INST_ROCKSTAR=0 # Install the Rockstar halo finder?
INST_SCIPY=0    # Install scipy?

# If you've got yt some other place, set this to point to it.
YT_DIR=""

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
    echo INST_ENZO=${INST_ENZO} >> ${CONFIG_FILE}
    echo INST_SQLITE3=${INST_SQLITE3} >> ${CONFIG_FILE}
    echo INST_PYX=${INST_PYX} >> ${CONFIG_FILE}
    echo INST_0MQ=${INST_0MQ} >> ${CONFIG_FILE}
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
        echo "$ sudo yum install gcc gcc-g++ gcc-gfortran make patch zip"
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

printf "%-15s = %s so I " "INST_HG" "${INST_HG}"
get_willwont ${INST_HG}
echo "be installing Mercurial"

printf "%-15s = %s so I " "INST_ENZO" "${INST_ENZO}"
get_willwont ${INST_ENZO}
echo "be checking out Enzo"

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
    case $LIB in
        *h5py*)
            BUILD_ARGS="--hdf5=${HDF5_DIR}"
            ;;
        *numpy*)
            if [ -e ${DEST_DIR}/lib/python2.7/site-packages/numpy/__init__.py ]
            then
                VER=$(${DEST_DIR}/bin/python -c 'from distutils.version import StrictVersion as SV; \
                                                 import numpy; print SV(numpy.__version__) < SV("1.8.0")')
                if [ $VER == "True" ]
                then
                    echo "Removing previous NumPy instance (see issue #889)"
                    rm -rf ${DEST_DIR}/lib/python2.7/site-packages/{numpy*,*.pth}
                fi
            fi
            ;;
        *)
            ;;
    esac
    cd $LIB
    ( ${DEST_DIR}/bin/python2.7 setup.py build ${BUILD_ARGS} $* 2>&1 ) 1>> ${LOG_FILE} || do_exit
    ( ${DEST_DIR}/bin/python2.7 setup.py install    2>&1 ) 1>> ${LOG_FILE} || do_exit
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
export PYTHONPATH=${DEST_DIR}/lib/python2.7/site-packages

mkdir -p ${DEST_DIR}/src
cd ${DEST_DIR}/src

CYTHON='Cython-0.20.2'
PYX='PyX-0.12.1'
PYTHON='Python-2.7.8'
BZLIB='bzip2-1.0.6'
FREETYPE_VER='freetype-2.4.12'
H5PY='h5py-2.3.1'
HDF5='hdf5-1.8.14'
IPYTHON='ipython-2.2.0'
LAPACK='lapack-3.4.2'
PNG=libpng-1.6.3
MATPLOTLIB='matplotlib-1.4.0'
MERCURIAL='mercurial-3.1'
NOSE='nose-1.3.4'
NUMPY='numpy-1.8.2'
PYTHON_HGLIB='python-hglib-1.0'
PYZMQ='pyzmq-14.3.1'
ROCKSTAR='rockstar-0.99.6'
SCIPY='scipy-0.14.0'
SQLITE='sqlite-autoconf-3071700'
SYMPY='sympy-0.7.5'
TORNADO='tornado-4.0.1'
ZEROMQ='zeromq-4.0.4'
ZLIB='zlib-1.2.8'

# Now we dump all our SHA512 files out.
echo '118e3ebd76f50bda8187b76654e65caab2c2c403df9b89da525c2c963dedc7b38d898ae0b92d44b278731d969a891eb3f7b5bcc138cfe3e037f175d4c87c29ec  Cython-0.20.2.tar.gz' > Cython-0.20.2.tar.gz.sha512
echo '4941f5aa21aff3743546495fb073c10d2657ff42b2aff401903498638093d0e31e344cce778980f28a7170c6d29eab72ac074277b9d4088376e8692dc71e55c1  PyX-0.12.1.tar.gz' > PyX-0.12.1.tar.gz.sha512
echo '4b05f0a490ddee37e8fc7970403bb8b72c38e5d173703db40310e78140d9d5c5732789d69c68dbd5605a623e4582f5b9671f82b8239ecdb34ad4261019dace6a  Python-2.7.8.tgz' > Python-2.7.8.tgz.sha512
echo '276bd9c061ec9a27d478b33078a86f93164ee2da72210e12e2c9da71dcffeb64767e4460b93f257302b09328eda8655e93c4b9ae85e74472869afbeae35ca71e  blas.tar.gz' > blas.tar.gz.sha512
echo '00ace5438cfa0c577e5f578d8a808613187eff5217c35164ffe044fbafdfec9e98f4192c02a7d67e01e5a5ccced630583ad1003c37697219b0f147343a3fdd12  bzip2-1.0.6.tar.gz' > bzip2-1.0.6.tar.gz.sha512
echo 'a296dfcaef7e853e58eed4e24b37c4fa29cfc6ac688def048480f4bb384b9e37ca447faf96eec7b378fd764ba291713f03ac464581d62275e28eb2ec99110ab6  reason-js-20120623.zip' > reason-js-20120623.zip.sha512
echo '609a68a3675087e0cc95268574f31e104549daa48efe15a25a33b8e269a93b4bd160f4c3e8178dca9c950ef5ca514b039d6fd1b45db6af57f25342464d0429ce  freetype-2.4.12.tar.gz' > freetype-2.4.12.tar.gz.sha512
echo 'f0da1d2ac855c02fb828444d719a1b23a580adb049335f3e732ace67558a125ac8cd3b3a68ac6bf9d10aa3ab19e4672b814eb28cc8c66910750c62efb655d744  h5py-2.3.1.tar.gz' > h5py-2.3.1.tar.gz.sha512
echo '4073fba510ccadaba41db0939f909613c9cb52ba8fb6c1062fc9118edc601394c75e102310be1af4077d07c9b327e6bbb1a6359939a7268dc140382d0c1e0199  hdf5-1.8.14.tar.gz' > hdf5-1.8.14.tar.gz.sha512
echo '4953bf5e9d6d5c6ad538d07d62b5b100fd86a37f6b861238501581c0059bd4655345ca05cf395e79709c38ce4cb9c6293f5d11ac0252a618ad8272b161140d13  ipython-2.2.0.tar.gz' > ipython-2.2.0.tar.gz.sha512
echo '8770214491e31f0a7a3efaade90eee7b0eb20a8a6ab635c5f854d78263f59a1849133c14ef5123d01023f0110cbb9fc6f818da053c01277914ae81473430a952  lapack-3.4.2.tar.gz' > lapack-3.4.2.tar.gz.sha512
echo '887582e5a22e4cde338aa8fec7a89f6dd31f2f02b8842735f00f970f64582333fa03401cea6d01704083403c7e8b7ebc26655468ce930165673b33efa4bcd586  libpng-1.6.3.tar.gz' > libpng-1.6.3.tar.gz.sha512
echo '60aa386639dec17b4f579955df60f2aa7c8ccd589b3490bb9afeb2929ea418d5d1a36a0b02b8d4a6734293076e9069429956c56cf8bd099b756136f2657cf9d4  matplotlib-1.4.0.tar.gz' > matplotlib-1.4.0.tar.gz.sha512
echo '1ee2fe7a241bf81087e55d9e4ee8fa986f41bb0655d4828d244322c18f3958a1f3111506e2df15aefcf86100b4fe530fcab2d4c041b5945599ed3b3a889d50f5  mercurial-3.1.tar.gz' > mercurial-3.1.tar.gz.sha512
echo '19499ab08018229ea5195cdac739d6c7c247c5aa5b2c91b801cbd99bad12584ed84c5cfaaa6fa8b4893a46324571a2f8a1988a1381f4ddd58390e597bd7bdc24  nose-1.3.4.tar.gz' > nose-1.3.4.tar.gz.sha512
echo '996e6b8e2d42f223e44660f56bf73eb8ab124f400d89218f8f5e4d7c9860ada44a4d7c54526137b0695c7a10f36e8834fbf0d42b7cb20bcdb5d5c245d673385c  numpy-1.8.2.tar.gz' > numpy-1.8.2.tar.gz.sha512
echo '9c0a61299779aff613131aaabbc255c8648f0fa7ab1806af53f19fbdcece0c8a68ddca7880d25b926d67ff1b9201954b207919fb09f6a290acb078e8bbed7b68  python-hglib-1.0.tar.gz' > python-hglib-1.0.tar.gz.sha512
echo '3d93a8fbd94fc3f1f90df68257cda548ba1adf3d7a819e7a17edc8681894003ac7ae6abd319473054340c11443a6a3817b931366fd7dae78e3807d549c544f8b  pyzmq-14.3.1.tar.gz' > pyzmq-14.3.1.tar.gz.sha512
echo 'ad1278740c1dc44c5e1b15335d61c4552b66c0439325ed6eeebc5872a1c0ba3fce1dd8509116b318d01e2d41da2ee49ec168da330a7fafd22511138b29f7235d  scipy-0.14.0.tar.gz' > scipy-0.14.0.tar.gz.sha512
echo '96f3e51b46741450bc6b63779c10ebb4a7066860fe544385d64d1eda52592e376a589ef282ace2e1df73df61c10eab1a0d793abbdaf770e60289494d4bf3bcb4  sqlite-autoconf-3071700.tar.gz' > sqlite-autoconf-3071700.tar.gz.sha512
echo '8a46e75abc3ed2388b5da9cb0e5874ae87580cf3612e2920b662d8f8eee8047efce5aa998eee96661d3565070b1a6b916c8bed74138b821f4e09115f14b6677d  sympy-0.7.5.tar.gz' > sympy-0.7.5.tar.gz.sha512
echo 'a4e0231e77ebbc2885bab648b292b842cb15c84d66a1972de18cb00fcc611eae2794b872f070ab7d5af32dd0c6c1773527fe1332bd382c1821e1f2d5d76808fb  tornado-4.0.1.tar.gz' > tornado-4.0.1.tar.gz.sha512
echo '7d70855d0537971841810a66b7a943a88304f6991ce445df19eea034aadc53dbce9d13be92bf44cfef1f3e19511a754eb01006a3968edc1ec3d1766ea4730cda  zeromq-4.0.4.tar.gz' > zeromq-4.0.4.tar.gz.sha512
echo 'ece209d4c7ec0cb58ede791444dc754e0d10811cbbdebe3df61c0fd9f9f9867c1c3ccd5f1827f847c005e24eef34fb5bf87b5d3f894d75da04f1797538290e4a  zlib-1.2.8.tar.gz' > zlib-1.2.8.tar.gz.sha512
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
get_ytproject $PYTHON.tgz
get_ytproject $NUMPY.tar.gz
get_ytproject $MATPLOTLIB.tar.gz
get_ytproject $MERCURIAL.tar.gz
get_ytproject $IPYTHON.tar.gz
get_ytproject $H5PY.tar.gz
get_ytproject $CYTHON.tar.gz
get_ytproject reason-js-20120623.zip
get_ytproject $NOSE.tar.gz
get_ytproject $PYTHON_HGLIB.tar.gz
get_ytproject $SYMPY.tar.gz
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

if [ ! -e $PYTHON/done ]
then
    echo "Installing Python.  This may take a while, but don't worry.  yt loves you."
    [ ! -e $PYTHON ] && tar xfz $PYTHON.tgz
    cd $PYTHON
    ( ./configure --prefix=${DEST_DIR}/ ${PYCONF_ARGS} 2>&1 ) 1>> ${LOG_FILE} || do_exit

    ( make ${MAKE_PROCS} 2>&1 ) 1>> ${LOG_FILE} || do_exit
    ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
    ( ln -sf ${DEST_DIR}/bin/python2.7 ${DEST_DIR}/bin/pyyt 2>&1 ) 1>> ${LOG_FILE}
    ( make clean 2>&1) 1>> ${LOG_FILE} || do_exit
    touch done
    cd ..
fi

export PYTHONPATH=${DEST_DIR}/lib/python2.7/site-packages/

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

echo "Installing distribute"
( ${DEST_DIR}/bin/python2.7 ${YT_DIR}/distribute_setup.py 2>&1 ) 1>> ${LOG_FILE} || do_exit

echo "Installing pip"
( ${DEST_DIR}/bin/easy_install-2.7 pip 2>&1 ) 1>> ${LOG_FILE} || do_exit

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
do_setup_py $H5PY
do_setup_py $CYTHON
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
( export PATH=$DEST_DIR/bin:$PATH ; ${DEST_DIR}/bin/python2.7 setup.py develop 2>&1 ) 1>> ${LOG_FILE} || do_exit
touch done
cd $MY_PWD

if !( ( ${DEST_DIR}/bin/python2.7 -c "import readline" 2>&1 )>> ${LOG_FILE})
then
    if !( ( ${DEST_DIR}/bin/python2.7 -c "import gnureadline" 2>&1 )>> ${LOG_FILE})
    then
        echo "Installing pure-python readline"
        ( ${DEST_DIR}/bin/pip install gnureadline 2>&1 ) 1>> ${LOG_FILE}
    fi
fi

if [ $INST_ENZO -eq 1 ]
then
    echo "Cloning a copy of Enzo."
    cd ${DEST_DIR}/src/
    ${HG_EXEC} clone https://bitbucket.org/enzo/enzo-stable ./enzo-hg-stable
    cd $MY_PWD
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
