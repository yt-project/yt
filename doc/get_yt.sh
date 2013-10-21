#
# Hi there!  Welcome to the yt installation script.
#
# This script is designed to create a fully isolated Python installation
# with the dependencies you need to run yt.
#
# This script is based on Conda, a distribution mechanism from Continuum
# Analytics.  The process is as follows:
#
#  1. Download the appropriate Conda installation package
#  2. Install Conda into the specified directory
#  3. Install yt-specific dependencies
#  4. Install yt
#
# There are a few options listed below, but by default, this will install
# everything.  At the end, it will tell you what to do to use yt.
#
# By default this will install yt from source.
#
# If you experience problems, please visit the Help section at 
# http://yt-project.org.
#
DEST_SUFFIX="yt-conda"
DEST_DIR="`pwd`/${DEST_SUFFIX/ /}"   # Installation location
BRANCH="yt" # This is the branch to which we will forcibly update.
INST_YT_SOURCE=1 # Do we do a source install of yt?

##################################################################
#                                                                #
# You will likely not have to modify anything below this region. #
#                                                                #
##################################################################

LOG_FILE="`pwd`/yt_install.log"

# Here is the idiom for redirecting to the log file:
# ( SOMECOMMAND 2>&1 ) 1>> ${LOG_FILE} || do_exit

MINICONDA_URLBASE="http://repo.continuum.io/miniconda"
MINICONDA_VERSION="1.9.1"
YT_RECIPE_REPO="https://bitbucket.org/yt_analysis/yt_conda/raw/default"

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

function log_cmd
{
    echo "EXECUTING:" >> ${LOG_FILE}
    echo "  $*" >> ${LOG_FILE}
    ( $* 2>&1 ) 1>> ${LOG_FILE} || do_exit
}

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

function get_ytrecipe {
    RDIR=${DEST_DIR}/src/yt-recipes/$1
    mkdir -p ${RDIR}
    pushd ${RDIR}
    log_cmd ${GETFILE} ${YT_RECIPE_REPO}/$1/meta.yaml
    log_cmd ${GETFILE} ${YT_RECIPE_REPO}/$1/build.sh
    NEW_PKG=`conda build --output ${RDIR}`
    log_cmd conda build --no-binstar-upload ${RDIR}
    log_cmd conda install ${NEW_PKG}
    popd
}


echo
echo
echo "========================================================================"
echo
echo "Hi there!  This is the yt installation script.  We're going to download"
echo "some stuff and install it to create a self-contained, isolated"
echo "environment for yt to run within."
echo
echo "This will install Miniconda from Continuum Analytics, the necessary"
echo "packages to run yt, and create a self-contained environment for you to"
echo "use yt.  Additionally, Conda itself provides the ability to install"
echo "many other packages that can be used for other purposes."
echo
MYOS=`uname -s`       # A guess at the OS
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
  echo "OS X 10.8.2: download Xcode 4.6.1 from the mac app store."
  echo "(search for Xcode)."
  echo "Additionally, you will have to manually install the Xcode"
  echo "command line tools, see:"
  echo "http://stackoverflow.com/questions/9353444"
  echo "Alternatively, download the Xcode command line tools from"
  echo "the Apple developer tools website."
  echo
  echo "NOTE: It's possible that the installation will fail, if so,"
  echo "please set the following environment variables, remove any"
  echo "broken installation tree, and re-run this script verbatim."
  echo
  echo "$ export CC=gcc"
  echo "$ export CXX=g++"
  echo
  MINICONDA_OS="MacOSX-x86_64"
fi
if [ "${MYOS##Linux}" != "${MYOS}" ]
then
  echo "Looks like you're on Linux."
  echo
  echo "Please make sure you have the developer tools for your OS installed."
  echo
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
    echo "  * chrpath"
    echo
    echo "You can accomplish this by executing:"
    echo
    echo "$ sudo zypper install -t pattern devel_C_C++"
    echo "$ sudo zypper install gcc-c++ libopenssl-devel libuuid-devel zip"
    echo "$ sudo zypper install chrpath"
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
    echo "  * chrpath"
    echo
    echo "You can accomplish this by executing:"
    echo
    echo "$ sudo apt-get install libssl-dev build-essential libncurses5 libncurses5-dev zip uuid-dev chrpath"
    echo
  fi
  echo
  echo "If you are running on a supercomputer or other module-enabled"
  echo "system, please make sure that the GNU module has been loaded."
  echo
  if [ "${MYOS##x86_64}" != "${MYOS}" ]
  then
    MINICONDA_OS="Linux-x86_64"
  elif [ "${MYOS##i386}" != "${MYOS}" ]
  then
    MINICONDA_OS="Linux-x86"
  else
    echo "Not sure which type of Linux you're on.  Going with x86_64."
    MINICONDA_OS="Linux-x86_64"
  fi
fi
echo
echo "If you'd rather not continue, hit Ctrl-C."
echo
echo "========================================================================"
echo
read -p "[hit enter] "
echo
echo "Awesome!  Here we go."
echo

MINICONDA_PKG=Miniconda-${MINICONDA_VERSION}-${MINICONDA_OS}.sh

if type -P wget &>/dev/null
then
    echo "Using wget"
    export GETFILE="wget -nv"
else
    echo "Using curl"
    export GETFILE="curl -sSO"
fi

echo
echo "Downloading ${MINICONDA_URLBASE}/${MINICONDA_PKG}"
echo "Downloading ${MINICONDA_URLBASE}/${MINICONDA_PKG}" >> ${LOG_FILE}
echo

${GETFILE} ${MINICONDA_URLBASE}/${MINICONDA_PKG} || do_exit

echo "Installing the Miniconda python environment."

log_cmd bash ./${MINICONDA_PKG} -b -p $DEST_DIR

# I don't think we need OR want this anymore:
#export LD_LIBRARY_PATH=${DEST_DIR}/lib:$LD_LIBRARY_PATH

# This we *do* need.
export PATH=${DEST_DIR}/bin:$PATH

echo "Installing the necessary packages for yt."
echo "This may take a while, but don't worry.  yt loves you."

declare -a YT_DEPS
YT_DEPS+=('python')
YT_DEPS+=('distribute')
YT_DEPS+=('libpng')
YT_DEPS+=('numpy')
YT_DEPS+=('pygments')
YT_DEPS+=('jinja2')
YT_DEPS+=('tornado')
YT_DEPS+=('pyzmq')
YT_DEPS+=('ipython')
YT_DEPS+=('sphinx')
YT_DEPS+=('h5py')
YT_DEPS+=('matplotlib')
YT_DEPS+=('cython')
YT_DEPS+=('nose')

# Here is our dependency list for yt
log_cmd conda config --system --add channels http://repo.continuum.io/pkgs/free
log_cmd conda config --system --add channels http://repo.continuum.io/pkgs/dev
log_cmd conda config --system --add channels http://repo.continuum.io/pkgs/gpl
log_cmd conda update --yes conda

echo "Current dependencies: ${YT_DEPS[@]}"
log_cmd echo "DEPENDENCIES" ${YT_DEPS[@]}
log_cmd conda install --yes ${YT_DEPS[@]}

echo "Installing mercurial."
get_ytrecipe mercurial

if [ $INST_YT_SOURCE -eq 0 ]
then
  echo "Installing yt as a package."
  get_ytrecipe yt
else
  # We do a source install.
  YT_DIR="${DEST_DIR}/src/yt-hg"
  export PNG_DIR=${DEST_DIR}
  export FTYPE_DIR=${DEST_DIR}
  export HDF5_DIR=${DEST_DIR}
  log_cmd hg clone -r ${BRANCH} https://bitbucket.org/yt_analysis/yt ${YT_DIR}
  pushd ${YT_DIR}
  log_cmd python setup.py develop
  popd
  log_cmd cp ${YT_DIR}/doc/activate ${DEST_DIR}/bin/activate 
  log_cmd sed -i.bak -e "s,__YT_DIR__,${DEST_DIR}," ${DEST_DIR}/bin/activate
  log_cmd cp ${YT_DIR}/doc/activate.csh ${DEST_DIR}/bin/activate.csh
  log_cmd sed -i.bak -e "s,__YT_DIR__,${DEST_DIR}," ${DEST_DIR}/bin/activate.csh
fi

echo
echo
echo "========================================================================"
echo
echo "yt and the Conda system are now installed in $DEST_DIR ."
echo
if [ $INST_YT_SOURCE -eq 0 ]
then
  echo "You must now modify your PATH variable by prepending:"
  echo 
  echo "   $DEST_DIR/bin"
  echo
  echo "For example, if you use bash, place something like this at the end"
  echo "of your ~/.bashrc :"
  echo
  echo "   export PATH=$DEST_DIR/bin:$PATH"
else
  echo "To run from this new installation, use the activate script for this "
  echo "environment."
  echo
  echo "    $ source $DEST_DIR/bin/activate"
  echo
  echo "This modifies the environment variables YT_DEST, PATH, PYTHONPATH, and"
  echo "LD_LIBRARY_PATH to match your new yt install.  If you use csh, just"
  echo "append .csh to the above."
fi
echo
echo "To get started with yt, check out the orientation:"
echo
echo "    http://yt-project.org/doc/orientation/"
echo
echo "or just activate your environment and run 'yt serve' to bring up the"
echo "yt GUI."
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
