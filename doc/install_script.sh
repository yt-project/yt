#
# Hi there!  Welcome to the yt installation script.
#
# First things first, if you experience problems, please visit the Help
# section at http://yt-project.org.
#
# This script is designed to create a fully isolated Python installation
# with the dependencies you need to run yt.
#
# If you would like to customize the yt installation, then please edit
# the following options.

# If you do not have a working compiler environment, use the following
# configuration:

INST_YT_SOURCE=0   # Should yt itself be installed from source?

# What follows are some other options that you may or may not need to change.

# If you've got a clone of the yt repository some other place, set this to
# point to it. The script will already check the current directory and the one
# above it in the tree.
YT_DIR=""

# These options can be set to customize the installation.

INST_GIT=1      # Install git or not?  If git is not already installed, yt
                # cannot be installed from source.
INST_EMBREE=0   # Install dependencies needed for Embree-accelerated ray tracing
INST_PYX=0      # Install PyX?  Sometimes PyX can be problematic without a
                # working TeX installation.
INST_SCIPY=0    # Install scipy?
INST_H5PY=1     # Install h5py?
INST_ASTROPY=0  # Install astropy?
INST_CARTOPY=0  # Install cartopy?
INST_NOSE=1     # Install nose?
INST_NETCDF4=1  # Install netcdf4 and its python bindings?
INST_POOCH=1    # Install pooch?

# This is the branch we will install from for INST_YT_SOURCE=1
BRANCH="main"

# These variables control which miniconda version is used

MINICONDA_URLBASE="http://repo.continuum.io/miniconda"
MINICONDA_VERSION="latest"

if [ ! -z "${CONDA_DEFAULT_ENV}" ]
then
    echo "Aborting the yt installation because you appear to already"
    echo "have a conda environment activated. Either deactivate it with:"
    echo
    echo "    $ source deactivate"
    echo
    echo "or install yt into your current environment with:"
    echo
    echo "    $ conda install -c conda-forge yt"
    echo
    exit 1
fi
DEST_SUFFIX="yt-conda"
if [ -n "${PYTHONPATH}" ]
then
    echo "WARNING WARNING WARNING WARNING WARNING WARNING WARNING"
    echo "*******************************************************"
    echo
    echo "The PYTHONPATH environment variable is set to:"
    echo
    echo "    $PYTHONPATH"
    echo
    echo "If dependencies of yt (numpy, scipy, matplotlib) are installed"
    echo "to this path, this may cause issues. Exit the install script"
    echo "with Ctrl-C and unset PYTHONPATH if you are unsure."
    echo "Hit enter to continue."
    echo
    echo "WARNING WARNING WARNING WARNING WARNING WARNING WARNING"
    echo "*******************************************************"
    read -p "[hit enter]"
fi

if [ -z "${DEST_DIR}" ]
then
    DEST_DIR="`pwd`/${DEST_SUFFIX/ /}"   # Installation location
fi

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
   echo
   echo "If you really want to do this, you must manually edit"
   echo "the script."
   echo "Sorry!"
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
    echo INST_YT_SOURCE=${INST_YT_SOURCE} > ${CONFIG_FILE}
    echo INST_GIT=${INST_GIT} >> ${CONFIG_FILE}
    echo INST_PYX=${INST_PYX} >> ${CONFIG_FILE}
    echo INST_SCIPY=${INST_SCIPY} >> ${CONFIG_FILE}
    echo INST_EMBREE=${INST_EMBREE} >> ${CONFIG_FILE}
    echo INST_H5PY=${INST_H5PY} >> ${CONFIG_FILE}
    echo INST_ASTROPY=${INST_ASTROPY} >> ${CONFIG_FILE}
    echo INST_CARTOPY=${INST_CARTOPY} >> ${CONFIG_FILE}
    echo INST_NOSE=${INST_NOSE} >> ${CONFIG_FILE}
    echo INST_POOCH=${INST_POOCH} >> ${CONFIG_FILE}

    echo YT_DIR=${YT_DIR} >> ${CONFIG_FILE}
}

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
    if [ "${MYOS##Darwin}" != "${MYOS}" ]
    then
        echo "Looks like you're running on MacOS."
        echo
        echo "NOTE: you must have the Xcode command line tools installed."
        echo
        echo "Download the appropriate version of Xcode from the"
        echo "Mac App Store (search for Xcode)."
        echo
        echo "Additionally, you will have to manually install the Xcode"
        echo "command line tools."
        echo
        echo "For MacOS 10.10 and newer the command line tools can be installed"
        echo "with the following command:"
        echo "    xcode-select --install"
        echo
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
        echo "  * git"
        echo
        echo "You can accomplish this by executing:"
        echo "$ sudo yum install gcc gcc-c++ gcc-gfortran make patch zip git"
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
        echo "  * git"
        echo
        echo "You can accomplish this by executing:"
        echo
        echo "$ sudo zypper install -t pattern devel_C_C++"
        echo "$ sudo zypper install git-core gcc-c++ libopenssl-devel libuuid-devel zip"
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
        echo "  * git"
        echo
        echo "You can accomplish this by executing:"
        echo
        echo "$ sudo apt-get install libssl-dev build-essential libncurses5 libncurses5-dev zip uuid-dev libfreetype6-dev tk-dev git"
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
}

function log_cmd
{
    echo "EXECUTING:" >> ${LOG_FILE}
    echo "  $*" >> ${LOG_FILE}
    ( $* 2>&1 ) 1>> ${LOG_FILE} || do_exit
}

# set paths needed for Embree

if [ $INST_EMBREE -ne 0 ]
then
    if [ $INST_YT_SOURCE -eq 0 ]
    then
        echo "yt must be compiled from source to install Embree support."
        echo "Please set INST_YT_SOURCE to 1 and re-run the install script."
        exit 1
    fi
    if [ `uname` = "Darwin" ]
    then
        EMBREE="embree-2.8.0.x86_64.macosx"
        EMBREE_URL="https://github.com/embree/embree/releases/download/v2.8.0/$EMBREE.tar.gz"
    elif [ `uname` = "Linux" ]
    then
        EMBREE="embree-2.8.0.x86_64.linux"
        EMBREE_URL="https://github.com/embree/embree/releases/download/v2.8.0/$EMBREE.tar.gz"
    else
        echo "Embree is not supported on this platform."
        echo "Set INST_EMBREE=0 and re-run the install script."
        exit 1
    fi
    PYEMBREE_URL="https://github.com/scopatz/pyembree/archive/master.zip"
fi

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

printf "%-18s = %s so I " "INST_YT_SOURCE" "${INST_YT_SOURCE}"
get_willwont ${INST_YT_SOURCE}
echo "be compiling yt from source"

printf "%-18s = %s so I " "INST_GIT" "${INST_GIT}"
get_willwont ${INST_GIT}
echo "be installing git"

printf "%-18s = %s so I " "INST_EMBREE" "${INST_EMBREE}"
get_willwont ${INST_EMBREE}
echo "be installing Embree"

printf "%-18s = %s so I " "INST_PYX" "${INST_PYX}"
get_willwont ${INST_PYX}
echo "be installing PyX"

printf "%-18s = %s so I " "INST_H5PY" "${INST_H5PY}"
get_willwont ${INST_H5PY}
echo "be installing h5py"

printf "%-18s = %s so I " "INST_ASTROPY" "${INST_ASTROPY}"
get_willwont ${INST_ASTROPY}
echo "be installing astropy"

printf "%-18s = %s so I " "INST_CARTOPY" "${INST_CARTOPY}"
get_willwont ${INST_CARTOPY}
echo "be installing cartopy"

printf "%-18s = %s so I " "INST_NOSE" "${INST_NOSE}"
get_willwont ${INST_NOSE}
echo "be installing nose"

printf "%-18s = %s so I " "INST_POOCH" "${INST_POOCH}"
get_willwont ${INST_POOCH}
echo "be installing pooch"

echo

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

if [ $INST_YT_SOURCE -ne 0 ]
then
   host_specific
fi

echo

if [[ $1 != "--yes" ]]
then
    read -p "[hit enter] "
fi

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

PYTHON_EXEC='python3'


if type -P curl &>/dev/null
then
    echo "Using curl"
    export GETFILE="curl -sSOL"
else
    echo "Using wget"
    export GETFILE="wget -nv"
fi

function test_install
{
    echo "Testing that yt can be imported"
    ( ${DEST_DIR}/bin/${PYTHON_EXEC} -c "import yt" 2>&1 ) 1>> ${LOG_FILE} || do_exit
}

ORIG_PWD=`pwd`

MYARCH=`uname -m`
MYOS=`uname -s`

if [ $MYOS = "Darwin" ]
then
    MINICONDA_OS="MacOSX"
    MINICONDA_ARCH="x86_64"
elif [ $MYOS = "Linux" ]
then
    MINICONDA_OS="Linux"
    if [ $MYARCH = "i386" ]
    then
        MINICONDA_ARCH="x86"
    elif [ $MYARCH = "i686"  ]
    then
        MINICONDA_ARCH="x86"
    elif [ $MYARCH = "x86_64"  ]
    then
        MINICONDA_ARCH="x86_64"
    elif [ $MYARCH = "ppc64le"  ]
    then
        MINICONDA_ARCH="ppc64le"
    else
        echo "Not sure which architecture you are running."
        echo "Going with x86_64 architecture."
        MINICONDA_ARCH="x86_64"
    fi
else
    echo "The yt install script is not supported on the ${MYOS}"
    echo "operating system."
    exit 1
fi

PY_VERSION='3'

MINICONDA_PKG="Miniconda${PY_VERSION}-${MINICONDA_VERSION}-${MINICONDA_OS}-${MINICONDA_ARCH}.sh"

echo
echo "Downloading ${MINICONDA_URLBASE}/${MINICONDA_PKG}"
echo

if [ -f ${MINICONDA_PKG} ]
then
    rm $MINICONDA_PKG
fi

echo "Installing the Miniconda python environment."

if [ -e ${DEST_DIR} ]
then
    rm -rf $DEST_DIR/*
else
    mkdir $DEST_DIR
fi

log_cmd ${GETFILE} ${MINICONDA_URLBASE}/${MINICONDA_PKG} || do_exit

log_cmd bash ./${MINICONDA_PKG} -b -p $DEST_DIR -f

# Need to set PATH so we use miniconda's python environment
export PATH=${DEST_DIR}/bin:$PATH

echo "Installing the necessary packages for yt."
echo "This may take a while, but don't worry.  yt loves you."

declare -a YT_DEPS
YT_DEPS+=('python')
YT_DEPS+=('setuptools')
YT_DEPS+=('numpy')
YT_DEPS+=('jupyter')
YT_DEPS+=('ipython')
YT_DEPS+=('sphinx')
if [ ${INST_GIT} -eq 1 ]
then
    YT_DEPS+=('git')
    YT_DEPS+=('gitpython')
fi
if [ $INST_H5PY -ne 0 ]
then
    YT_DEPS+=('h5py')
fi
YT_DEPS+=('matplotlib')
YT_DEPS+=('cython')
if [ $INST_NOSE -ne 0 ]
then
    YT_DEPS+=('nose')
fi
if [ $INST_SCIPY -ne 0 ]
then
    YT_DEPS+=('scipy')
fi
if [ $INST_ASTROPY -ne 0 ]
then
    YT_DEPS+=('astropy')
fi
if [ $INST_CARTOPY -ne 0 ]
then
    YT_DEPS+=('cartopy')
fi
if [ $INST_POOCH -ne 0 ]
then
    YT_DEPS+=('pooch')
fi
YT_DEPS+=('conda-build')
YT_DEPS+=('sympy')

if [ $INST_NETCDF4 -eq 1 ]
then
    YT_DEPS+=('netcdf4')
fi

log_cmd ${DEST_DIR}/bin/conda update --yes conda

if [ $INST_GIT -eq 1 ]
then
    GIT_EXE=${DEST_DIR}/bin/git
else
    if type -P git &>/dev/null
    then
        GIT_EXE="git"
    else
        if [ $INST_YT_SOURCE -eq 1 ]
        then
            echo "Cannot find git. Please install git or set INST_GIT=1."
            do_exit
        fi
    fi
fi

log_cmd echo "DEPENDENCIES" ${YT_DEPS[@]}
for YT_DEP in "${YT_DEPS[@]}"; do
    echo "Installing $YT_DEP"
    log_cmd ${DEST_DIR}/bin/conda install -c conda-forge --yes ${YT_DEP}
done

if [ $INST_YT_SOURCE -eq 1 ]
then
    log_cmd ${GIT_EXE} clone https://github.com/yt-project/yt_conda ${DEST_DIR}/src/yt_conda
fi

if [ $INST_EMBREE -eq 1 ]
then
    echo "Installing Embree"
    if [ ! -d ${DEST_DIR}/src ]
    then
        mkdir ${DEST_DIR}/src
    fi
    cd ${DEST_DIR}/src
    ( ${GETFILE} "$EMBREE_URL" 2>&1 ) 1>> ${LOG_FILE} || do_exit
    log_cmd tar xfz ${EMBREE}.tar.gz
    log_cmd mv ${DEST_DIR}/src/${EMBREE}/include/embree2 ${DEST_DIR}/include
    log_cmd mv ${DEST_DIR}/src/${EMBREE}/lib/lib*.* ${DEST_DIR}/lib
    if [ `uname` = "Darwin" ]
    then
        ln -s ${DEST_DIR}/lib/libembree.2.dylib ${DEST_DIR}/lib/libembree.dylib
        install_name_tool -id ${DEST_DIR}/lib/libembree.2.dylib ${DEST_DIR}/lib/libembree.2.dylib
    else
        ln -s ${DEST_DIR}/lib/libembree.so.2 ${DEST_DIR}/lib/libembree.so
    fi

    echo "Installing pyembree from source"
    ( ${GETFILE} "$PYEMBREE_URL" 2>&1 ) 1>> ${LOG_FILE} || do_exit
    log_cmd unzip ${DEST_DIR}/src/master.zip
    pushd ${DEST_DIR}/src/pyembree-master &> /dev/null
    log_cmd ${DEST_DIR}/bin/${PYTHON_EXEC} setup.py install build_ext -I${DEST_DIR}/include -L${DEST_DIR}/lib
    popd &> /dev/null
fi

# conda doesn't package pyx, so we install manually with pip
if [ $INST_PYX -eq 1 ]
then
    log_cmd ${DEST_DIR}/bin/pip install pyx
fi

if [ $INST_YT_SOURCE -eq 0 ]
then
    echo "Installing yt"
    log_cmd ${DEST_DIR}/bin/conda install -c conda-forge --yes yt
else
    echo "Building yt from source"
    if [ -z "$YT_DIR" ]
    then
        if [ -e $ORIG_PWD/yt/mods.py ]
        then
            YT_DIR="$ORIG_PWD"
        elif [ -e $ORIG_PWD/../yt/mods.py ]
        then
            YT_DIR=$(dirname $ORIG_PWD)
        else
            YT_DIR="${DEST_DIR}/src/yt-git"
            log_cmd ${GIT_EXE} clone https://github.com/yt-project/yt ${YT_DIR}
            log_cmd ${GIT_EXE} -C ${YT_DIR} checkout ${BRANCH}
        fi
        echo Setting YT_DIR=${YT_DIR}
    else
        if [ ! -e $YT_DIR/.git ]
        then
            echo "$YT_DIR is not a clone of the yt git repository, exiting"
            do_exit
        fi
    fi
    if [ $INST_EMBREE -eq 1 ]
    then
        echo $DEST_DIR > ${YT_DIR}/embree.cfg
    fi
    pushd ${YT_DIR} &> /dev/null
    ( ${DEST_DIR}/bin/${PYTHON_EXEC} setup.py develop 2>&1) 1>> ${LOG_FILE} || do_exit
    popd &> /dev/null
fi

test_install

echo
echo
echo "========================================================================"
echo
echo "yt and the Conda system are now installed in $DEST_DIR"
echo
echo "To get started with yt, check out the orientation:"
echo
echo "    http://yt-project.org/doc/orientation/"
echo
echo "For support, see the website and join the mailing list:"
echo
echo "    http://yt-project.org/"
echo "    http://yt-project.org/data/      (Sample data)"
echo "    http://yt-project.org/doc/       (Docs)"
echo
echo "    https://mail.python.org/archives/list/yt-users@python.org/"
echo
echo "You must now prepend the following folder to your PATH environment variable:"
echo
echo "    $DEST_DIR/bin"
echo
echo "On Bash-style shells you can copy/paste the following command to "
echo "temporarily activate the yt installation:"
echo
echo "    export PATH=$DEST_DIR/bin:\$PATH"
echo
echo "and on csh-style shells:"
echo
echo "    setenv PATH $DEST_DIR/bin:\$PATH"
echo
echo "You can also update the init file appropriate for your shell"
echo "(e.g. .bashrc, .bash_profile, .cshrc, or .zshrc) to include"
echo "the same command."
echo
echo "========================================================================"
echo
echo "Oh, look at me, still talking when there's science to do!"
echo "Good luck, and email the mailing list if you run into any problems."
