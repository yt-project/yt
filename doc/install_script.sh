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
# And, feel free to drop me a line: matthewturk@gmail.com
#

DEST_SUFFIX="yt-`uname -m`"
DEST_DIR="`pwd`/${DEST_SUFFIX/ /}"   # Installation location
BRANCH="yt" # This is the branch to which we will forcibly update.

if [ ${REINST_YT} -eq 1 ] && [ -n ${YT_DEST} ]
then
    DEST_DIR=${YT_DEST}
fi

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
INST_FTYPE=1    # Install FreeType2 locally?
INST_ENZO=0     # Clone a copy of Enzo?
INST_SQLITE3=1  # Install a local version of SQLite3?
INST_PYX=0      # Install PyX?  Sometimes PyX can be problematic without a
                # working TeX installation.
INST_0MQ=1      # Install 0mq (for IPython) and affiliated bindings?
INST_ROCKSTAR=0 # Install the Rockstar halo finder?
INST_SCIPY=0    # Install scipy?

# If you've got YT some other place, set this to point to it.
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
    if [ "${MYHOSTLONG%%ranger}" != "${MYHOSTLONG}" ]
    then
        echo "Looks like you're on Ranger."
        echo
        echo "NOTE: YOU MUST BE IN THE GNU PROGRAMMING ENVIRONMENT"
        echo "These commands should take care of that for you:"
        echo
        echo "   $ module unload mvapich2"
        echo "   $ module swap pgi gcc"
        echo "   $ module load mvapich2"
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
    if [ "${MYOS##Darwin}" != "${MYOS}" ]
    then
        echo "Looks like you're running on Mac OSX."
        echo
        echo "NOTE: you must have the Xcode command line tools installed."
        echo
        echo "OS X 10.5: download Xcode 3.0 from the mac developer tools"
        echo "website"
        echo
        echo "OS X 10.6: download Xcode 3.2 from the mac developer tools"
        echo "website"
        echo
        echo "OS X 10.7: download Xcode 4.0 from the mac app store or"
        echo "alternatively download the Xcode command line tools from"
        echo "the mac developer tools website"
        echo
        echo "NOTE: You may have problems if you are running OSX 10.6 (Snow"
        echo "Leopard) or newer.  If you do, please set the following"
        echo "environment variables, remove any broken installation tree, and"
        echo "re-run this script verbatim."
        echo
        echo "$ export CC=gcc-4.2"
        echo "$ export CXX=g++-4.2"
        echo
        OSX_VERSION=`sw_vers -productVersion`
        if [ "${OSX_VERSION##10.8}" != "${OSX_VERSION}" ]
        then
            MPL_SUPP_CFLAGS="${MPL_SUPP_CFLAGS} -mmacosx-version-min=10.7"
            MPL_SUPP_CXXFLAGS="${MPL_SUPP_CXXFLAGS} -mmacosx-version-min=10.7"
        fi
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
get_willwont ${INST_PYX}
echo "be installing scipy"

printf "%-15s = %s so I " "INST_0MQ" "${INST_0MQ}"
get_willwont ${INST_0MQ}
echo "be installing ZeroMQ"

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
    echo "Installing $1 (arguments: '$*')"
    [ ! -e $1/extracted ] && tar xfz $1.tar.gz
    touch $1/extracted
    cd $1
    if [ ! -z `echo $1 | grep h5py` ]
    then
        shift
	( ${DEST_DIR}/bin/python2.7 setup.py build --hdf5=${HDF5_DIR} $* 2>&1 ) 1>> ${LOG_FILE} || do_exit
    else
        shift
        ( ${DEST_DIR}/bin/python2.7 setup.py build   $* 2>&1 ) 1>> ${LOG_FILE} || do_exit
    fi
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
    echo "Downloading $1 from yt-project.org"
    [ -e $1 ] && return
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
echo 'de6d8c6ea849f0206d219303329a0276b3cce7c051eec34377d42aacbe0a4f47ac5145eb08966a338ecddd2b83c8f787ca9956508ad5c39ee2088ad875166410  xray_emissivity.h5' > xray_emissivity.h5.sha512
get_ytdata xray_emissivity.h5

mkdir -p ${DEST_DIR}/src
cd ${DEST_DIR}/src

# Now we dump all our SHA512 files out.

echo 'eda1b8090e5e21e7e039ef4dd03de186a7b416df9d5a4e4422abeeb4d51383b9a6858e1ac4902d8e5010f661b295bbb2452c43c8738be668379b4eb4835d0f61  Cython-0.17.1.tar.gz' > Cython-0.17.1.tar.gz.sha512
echo '44eea803870a66ff0bab08d13a8b3388b5578ebc1c807d1d9dca0a93e6371e91b15d02917a00b3b20dc67abb5a21dabaf9b6e9257a561f85eeff2147ac73b478  PyX-0.11.1.tar.gz' > PyX-0.11.1.tar.gz.sha512
echo 'b981f8464575bb24c297631c87a3b9172312804a0fc14ce1fa7cb41ce2b0d2fd383cd1c816d6e10c36467d18bf9492d6faf557c81c04ff3b22debfa93f30ad0b  Python-2.7.3.tgz' > Python-2.7.3.tgz.sha512
echo 'c017d3d59dd324ac91af0edc178c76b60a5f90fbb775cf843e39062f95bd846238f2c53705f8890ed3f34bc0e6e75671a73d13875eb0287d6201cb45f0a2d338  bzip2-1.0.5.tar.gz' > bzip2-1.0.5.tar.gz.sha512
echo 'a296dfcaef7e853e58eed4e24b37c4fa29cfc6ac688def048480f4bb384b9e37ca447faf96eec7b378fd764ba291713f03ac464581d62275e28eb2ec99110ab6  reason-js-20120623.zip' > reason-js-20120623.zip.sha512
echo 'b519218f93946400326e9b656669269ecb3e5232b944e18fbc3eadc4fe2b56244d68aae56d6f69042b4c87c58c881ee2aaa279561ea0f0f48d5842155f4de9de  freetype-2.4.4.tar.gz' > freetype-2.4.4.tar.gz.sha512
echo 'b3290c498191684781ca5286ab454eb1bd045e8d894f5b86fb86beb88f174e22ac3ab008fb02d6562051d9fa6a9593920cab433223f6d5473999913223b8e183  h5py-2.1.0.tar.gz' > h5py-2.1.0.tar.gz.sha512
echo 'c68a425bacaa7441037910b9166f25b89e1387776a7749a5350793f89b1690350df5f018060c31d03686e7c3ed2aa848bd2b945c96350dc3b6322e087934783a  hdf5-1.8.9.tar.gz' > hdf5-1.8.9.tar.gz.sha512
echo 'dbefad00fa34f4f21dca0f1e92e95bd55f1f4478fa0095dcf015b4d06f0c823ff11755cd777e507efaf1c9098b74af18f613ec9000e5c3a5cc1c7554fb5aefb8  libpng-1.5.12.tar.gz' > libpng-1.5.12.tar.gz.sha512
echo '5b1a0fb52dcb21ca5f0ab71c8a49550e1e8cf633552ec6598dc43f0b32c03422bf5af65b30118c163231ecdddfd40846909336f16da318959106076e80a3fad0  matplotlib-1.2.0.tar.gz' > matplotlib-1.2.0.tar.gz.sha512
echo '52d1127de2208aaae693d16fef10ffc9b8663081bece83b7597d65706e9568af3b9e56bd211878774e1ebed92e21365ee9c49602a0ff5e48f89f12244d79c161  mercurial-2.4.tar.gz' > mercurial-2.4.tar.gz.sha512
echo 'de3dd37f753614055dcfed910e9886e03688b8078492df3da94b1ec37be796030be93291cba09e8212fffd3e0a63b086902c3c25a996cf1439e15c5b16e014d9  numpy-1.6.1.tar.gz' > numpy-1.6.1.tar.gz.sha512
echo '5ad681f99e75849a5ca6f439c7a19bb51abc73d121b50f4f8e4c0da42891950f30407f761a53f0fe51b370b1dbd4c4f5a480557cb2444c8c7c7d5412b328a474  sqlite-autoconf-3070500.tar.gz' > sqlite-autoconf-3070500.tar.gz.sha512
echo 'edae735960279d92acf58e1f4095c6392a7c2059b8f1d2c46648fc608a0fb06b392db2d073f4973f5762c034ea66596e769b95b3d26ad963a086b9b2d09825f2  zlib-1.2.3.tar.bz2' > zlib-1.2.3.tar.bz2.sha512
echo '05ac335727a2c3036f31a2506fdd2615aa436bfbe2f81799fe6c51bffe2591ad6a8427f3b25c34e7e709fb4e7607a0589dc7a22185c1f9b894e90de6711a88aa  ipython-0.13.1.tar.gz' > ipython-0.13.1.tar.gz.sha512
echo 'fb3cf421b2dc48c31956b3e3ee4ab6ebc743deec3bf626c2238a1996c8c51be87260bd6aa662793a1f0c34dcda9b3146763777bb162dfad6fec4ca7acc403b2e  zeromq-2.2.0.tar.gz' > zeromq-2.2.0.tar.gz.sha512
echo 'd761b492352841cdc125d9f0c99ee6d6c435812472ea234728b7f0fb4ad1048e1eec9b399df2081fbc926566f333f7780fedd0ce23255a6633fe5c60ed15a6af  pyzmq-2.1.11.tar.gz' > pyzmq-2.1.11.tar.gz.sha512
echo '57fa5e57dfb98154a42d2d477f29401c2260ae7ad3a8128a4098b42ee3b35c54367b1a3254bc76b9b3b14b4aab7c3e1135858f68abc5636daedf2f01f9b8a3cf  tornado-2.2.tar.gz' > tornado-2.2.tar.gz.sha512
echo '1332e3d5465ca249c357314cf15d2a4e5e83a941841021b8f6a17a107dce268a7a082838ade5e8db944ecde6bfb111211ab218aa414ee90aafbb81f1491b3b93  Forthon-0.8.10.tar.gz' > Forthon-0.8.10.tar.gz.sha512
echo 'c13116c1f0547000cc565e15774687b9e884f8b74fb62a84e578408a868a84961704839065ae4f21b662e87f2aaedf6ea424ea58dfa9d3d73c06281f806d15dd  nose-1.2.1.tar.gz' > nose-1.2.1.tar.gz.sha512
echo '73de2c99406a38f85273931597525cec4ebef55b93712adca3b0bfea8ca3fc99446e5d6495817e9ad55cf4d48feb7fb49734675c4cc8938db8d4a5225d30eca7  python-hglib-0.2.tar.gz' > python-hglib-0.2.tar.gz.sha512
echo 'ffc602eb346717286b3d0a6770c60b03b578b3cf70ebd12f9e8b1c8c39cdb12ef219ddaa041d7929351a6b02dbb8caf1821b5452d95aae95034cbf4bc9904a7a  sympy-0.7.2.tar.gz' > sympy-0.7.2.tar.gz.sha512
echo '172f2bc671145ebb0add2669c117863db35851fb3bdb192006cd710d4d038e0037497eb39a6d01091cb923f71a7e8982a77b6e80bf71d6275d5d83a363c8d7e5  rockstar-0.99.6.tar.gz' > rockstar-0.99.6.tar.gz.sha512
echo 'd4fdd62f2db5285cd133649bd1bfa5175cb9da8304323abd74e0ef1207d55e6152f0f944da1da75f73e9dafb0f3bb14efba3c0526c732c348a653e0bd223ccfa  scipy-0.11.0.tar.gz' > scipy-0.11.0.tar.gz.sha512
echo '276bd9c061ec9a27d478b33078a86f93164ee2da72210e12e2c9da71dcffeb64767e4460b93f257302b09328eda8655e93c4b9ae85e74472869afbeae35ca71e  blas.tar.gz' > blas.tar.gz.sha512
echo '8770214491e31f0a7a3efaade90eee7b0eb20a8a6ab635c5f854d78263f59a1849133c14ef5123d01023f0110cbb9fc6f818da053c01277914ae81473430a952  lapack-3.4.2.tar.gz' > lapack-3.4.2.tar.gz.sha512
# Individual processes
[ -z "$HDF5_DIR" ] && get_ytproject hdf5-1.8.9.tar.gz
[ $INST_ZLIB -eq 1 ] && get_ytproject zlib-1.2.3.tar.bz2 
[ $INST_BZLIB -eq 1 ] && get_ytproject bzip2-1.0.5.tar.gz
[ $INST_PNG -eq 1 ] && get_ytproject libpng-1.5.12.tar.gz
[ $INST_FTYPE -eq 1 ] && get_ytproject freetype-2.4.4.tar.gz
[ $INST_SQLITE3 -eq 1 ] && get_ytproject sqlite-autoconf-3070500.tar.gz
[ $INST_PYX -eq 1 ] && get_ytproject PyX-0.11.1.tar.gz
[ $INST_0MQ -eq 1 ] && get_ytproject zeromq-2.2.0.tar.gz
[ $INST_0MQ -eq 1 ] && get_ytproject pyzmq-2.1.11.tar.gz
[ $INST_0MQ -eq 1 ] && get_ytproject tornado-2.2.tar.gz
[ $INST_SCIPY -eq 1 ] && get_ytproject scipy-0.11.0.tar.gz
[ $INST_SCIPY -eq 1 ] && get_ytproject blas.tar.gz
[ $INST_SCIPY -eq 1 ] && get_ytproject lapack-3.4.2.tar.gz
get_ytproject Python-2.7.3.tgz
get_ytproject numpy-1.6.1.tar.gz
get_ytproject matplotlib-1.2.0.tar.gz
get_ytproject mercurial-2.4.tar.gz
get_ytproject ipython-0.13.1.tar.gz
get_ytproject h5py-2.1.0.tar.gz
get_ytproject Cython-0.17.1.tar.gz
get_ytproject reason-js-20120623.zip
get_ytproject Forthon-0.8.10.tar.gz
get_ytproject nose-1.2.1.tar.gz 
get_ytproject python-hglib-0.2.tar.gz
get_ytproject sympy-0.7.2.tar.gz
get_ytproject rockstar-0.99.6.tar.gz
if [ $INST_BZLIB -eq 1 ]
then
    if [ ! -e bzip2-1.0.5/done ]
    then
        [ ! -e bzip2-1.0.5 ] && tar xfz bzip2-1.0.5.tar.gz
        echo "Installing BZLIB"
        cd bzip2-1.0.5
        if [ `uname` = "Darwin" ] 
        then
            if [ -z "${CC}" ] 
            then
                sed -i.bak 's/soname/install_name/' Makefile-libbz2_so
            else
                sed -i.bak -e 's/soname/install_name/' -e "s/CC=gcc/CC=${CC}/" Makefile-libbz2_so 
            fi
        fi
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
    if [ ! -e libpng-1.5.12/done ]
    then
        [ ! -e libpng-1.5.12 ] && tar xfz libpng-1.5.12.tar.gz
        echo "Installing PNG"
        cd libpng-1.5.12
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
    if [ ! -e freetype-2.4.4/done ]
    then
        [ ! -e freetype-2.4.4 ] && tar xfz freetype-2.4.4.tar.gz
        echo "Installing FreeType2"
        cd freetype-2.4.4
        ( ./configure CFLAGS=-I${DEST_DIR}/include --prefix=${DEST_DIR}/ 2>&1 ) 1>> ${LOG_FILE} || do_exit
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
    if [ ! -e hdf5-1.8.9/done ]
    then
        [ ! -e hdf5-1.8.9 ] && tar xfz hdf5-1.8.9.tar.gz
        echo "Installing HDF5"
        cd hdf5-1.8.9
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
    if [ ! -e sqlite-autoconf-3070500/done ]
    then
        [ ! -e sqlite-autoconf-3070500 ] && tar xfz sqlite-autoconf-3070500.tar.gz
        echo "Installing SQLite3"
        cd sqlite-autoconf-3070500
        ( ./configure --prefix=${DEST_DIR}/ 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make ${MAKE_PROCS} install 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make clean 2>&1) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
fi

if [ ! -e Python-2.7.3/done ]
then
    echo "Installing Python.  This may take a while, but don't worry.  YT loves you."
    [ ! -e Python-2.7.3 ] && tar xfz Python-2.7.3.tgz
    cd Python-2.7.3
    ( ./configure --prefix=${DEST_DIR}/ 2>&1 ) 1>> ${LOG_FILE} || do_exit

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
    echo "Installing Mercurial."
    do_setup_py mercurial-2.4
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
    do_setup_py numpy-1.6.1 ${NUMPY_ARGS}
else
    if [ ! -e scipy-0.11.0/done ]
    then
	if [ ! -e BLAS/done ]
	then
	    tar xfz blas.tar.gz
	    echo "Building BLAS"
	    cd BLAS
	    gfortran -O2 -fPIC -fno-second-underscore -c *.f
	    ar r libfblas.a *.o 1>> ${LOG_FILE}
	    ranlib libfblas.a 1>> ${LOG_FILE}
	    rm -rf *.o
	    touch done
	    cd ..
	fi
	if [ ! -e lapack-3.4.2/done ]
	then
	    tar xfz lapack-3.4.2.tar.gz
	    echo "Building LAPACK"
	    cd lapack-3.4.2/
	    cp INSTALL/make.inc.gfortran make.inc
	    make lapacklib CFLAGS=-fPIC LDFLAGS=-fPIC 1>> ${LOG_FILE} || do_exit
	    touch done
	    cd ..
	fi
    fi
    export BLAS=$PWD/BLAS/libfblas.a
    export LAPACK=$PWD/lapack-3.4.2/liblapack.a    
    do_setup_py numpy-1.6.1 ${NUMPY_ARGS}
    do_setup_py scipy-0.11.0 ${NUMPY_ARGS}
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
mkdir -p ${DEST_DIR}/src/matplotlib-1.2.0
echo "[directories]" >> ${DEST_DIR}/src/matplotlib-1.2.0/setup.cfg
echo "basedirlist = ${DEST_DIR}" >> ${DEST_DIR}/src/matplotlib-1.2.0/setup.cfg
do_setup_py matplotlib-1.2.0
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
    if [ ! -e zeromq-2.2.0/done ]
    then
        [ ! -e zeromq-2.2.0 ] && tar xfz zeromq-2.2.0.tar.gz
        echo "Installing ZeroMQ"
        cd zeromq-2.2.0
        ( ./configure --prefix=${DEST_DIR}/ 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make install 2>&1 ) 1>> ${LOG_FILE} || do_exit
        ( make clean 2>&1) 1>> ${LOG_FILE} || do_exit
        touch done
        cd ..
    fi
    do_setup_py pyzmq-2.1.11 --zmq=${DEST_DIR}
    do_setup_py tornado-2.2
fi

do_setup_py ipython-0.13.1
do_setup_py h5py-2.1.0
do_setup_py Cython-0.17.1
do_setup_py Forthon-0.8.10
do_setup_py nose-1.2.1
do_setup_py python-hglib-0.2
do_setup_py sympy-0.7.2
[ $INST_PYX -eq 1 ] && do_setup_py PyX-0.11.1

# Now we build Rockstar and set its environment variable.
if [ $INST_ROCKSTAR -eq 1 ]
then
    if [ ! -e Rockstar/done ]
    then
        [ ! -e Rockstar ] && tar xfz rockstar-0.99.6.tar.gz
        echo "Building Rockstar"
        cd Rockstar
        ( make lib 2>&1 ) 1>> ${LOG_FILE} || do_exit
        cp librockstar.so ${DEST_DIR}/lib
        ROCKSTAR_DIR=${DEST_DIR}/src/Rockstar
        echo $ROCKSTAR_DIR > ${YT_DIR}/rockstar.cfg
        touch done
        cd ..
    fi
fi

echo "Doing yt update, wiping local changes and updating to branch ${BRANCH}"
MY_PWD=`pwd`
cd $YT_DIR
( ${HG_EXEC} pull 2>1 && ${HG_EXEC} up -C 2>1 ${BRANCH} 2>&1 ) 1>> ${LOG_FILE}

echo "Building Fortran kD-tree module."
cd yt/utilities/kdtree
( make 2>&1 ) 1>> ${LOG_FILE}
cd ../../..

echo "Installing yt"
echo $HDF5_DIR > hdf5.cfg
[ $INST_PNG -eq 1 ] && echo $PNG_DIR > png.cfg
[ $INST_FTYPE -eq 1 ] && echo $FTYPE_DIR > freetype.cfg
( ${DEST_DIR}/bin/python2.7 setup.py develop 2>&1 ) 1>> ${LOG_FILE} || do_exit
touch done
cd $MY_PWD

if !(${DEST_DIR}/bin/python2.7 -c "import readline" >> ${LOG_FILE})
then
    echo "Installing pure-python readline"
    ${DEST_DIR}/bin/pip install readline 1>> ${LOG_FILE}
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
    echo "    http://yt-project.org/doc/orientation/"
    echo
    echo "or just activate your environment and run 'yt serve' to bring up the"
    echo "yt GUI."
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
