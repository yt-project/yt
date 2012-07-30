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

DEST_SUFFIX="yt-`uname -m`"
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
INST_FTYPE=1    # Install FreeType2 locally?
INST_ENZO=0     # Clone a copy of Enzo?
INST_SQLITE3=1  # Install a local version of SQLite3?
INST_PYX=0      # Install PyX?  Sometimes PyX can be problematic without a
                # working TeX installation.
INST_0MQ=1      # Install 0mq (for IPython) and affiliated bindings?

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
        echo
        echo "You can accomplish this by executing:"
        echo
        echo "$ sudo apt-get install libssl-dev build-essential libncurses5 libncurses5-dev zip"
        echo
    fi
    if [ ! -z "${CFLAGS}" ]
    then
        echo "******************************************"
        echo "******************************************"
        echo "**                                      **"
        echo "**    Your CFLAGS is not empty.         **"
        echo "**    This can beak h5py compilation.   **"
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

ORIG_PWD=`pwd`

if [ -z "${DEST_DIR}" ]
then
    echo "Edit this script, set the DEST_DIR parameter and re-run."
    exit 1
fi

mkdir -p ${DEST_DIR}/src
cd ${DEST_DIR}/src

# Now we dump all our SHA512 files out.

echo '2c1933ab31246b4f4eba049d3288156e0a72f1730604e3ed7357849967cdd329e4647cf236c9442ecfb06d0aff03e6fc892a7ba2a5c1cf5c011b7ab9c619acec  Cython-0.16.tar.gz' > Cython-0.16.tar.gz.sha512
echo 'b8a12bf05b3aafa71135e47da81440fd0f16a4bd91954bc5615ad3d3b7f9df7d5a7d5620dc61088dc6b04952c5c66ebda947a4cfa33ed1be614c8ca8c0f11dff  PhiloGL-1.4.2.zip' > PhiloGL-1.4.2.zip.sha512
echo '44eea803870a66ff0bab08d13a8b3388b5578ebc1c807d1d9dca0a93e6371e91b15d02917a00b3b20dc67abb5a21dabaf9b6e9257a561f85eeff2147ac73b478  PyX-0.11.1.tar.gz' > PyX-0.11.1.tar.gz.sha512
echo '1a754d560bfa433f0960ab3b5a62edb5f291be98ec48cf4e5941fa5b84139e200b87a52efbbd6fa4a76d6feeff12439eed3e7a84db4421940d1bbb576f7a684e  Python-2.7.2.tgz' > Python-2.7.2.tgz.sha512
echo 'c017d3d59dd324ac91af0edc178c76b60a5f90fbb775cf843e39062f95bd846238f2c53705f8890ed3f34bc0e6e75671a73d13875eb0287d6201cb45f0a2d338  bzip2-1.0.5.tar.gz' > bzip2-1.0.5.tar.gz.sha512
echo 'a296dfcaef7e853e58eed4e24b37c4fa29cfc6ac688def048480f4bb384b9e37ca447faf96eec7b378fd764ba291713f03ac464581d62275e28eb2ec99110ab6  reason-js-20120623.zip' > reason-js-20120623.zip.sha512
echo 'b519218f93946400326e9b656669269ecb3e5232b944e18fbc3eadc4fe2b56244d68aae56d6f69042b4c87c58c881ee2aaa279561ea0f0f48d5842155f4de9de  freetype-2.4.4.tar.gz' > freetype-2.4.4.tar.gz.sha512
echo '1531789e0a77d4829796d18552a4de7aecae7e8b63763a7951a8091921995800740fe03e72a7dbd496a5590828131c5f046ddead695e5cba79343b8c205148d1  h5py-2.0.1.tar.gz' > h5py-2.0.1.tar.gz.sha512
echo '9644896e4a84665ad22f87eb885cbd4a0c60a5c30085d5dd5dba5f3c148dbee626f0cb01e59a7995a84245448a3f1e9ba98687d3f10250e2ee763074ed8ddc0e  hdf5-1.8.7.tar.gz' > hdf5-1.8.7.tar.gz.sha512
echo 'dbefad00fa34f4f21dca0f1e92e95bd55f1f4478fa0095dcf015b4d06f0c823ff11755cd777e507efaf1c9098b74af18f613ec9000e5c3a5cc1c7554fb5aefb8  libpng-1.5.12.tar.gz' > libpng-1.5.12.tar.gz.sha512
echo 'f5ab95c29ef6958096970265a6079f0eb8c43a500924346c4a6c6eb89d9110eeeb6c34a53715e71240e82ded2b76a7b8d5a9b05a07baa000b2926718264ad8ff  matplotlib-1.1.0.tar.gz' > matplotlib-1.1.0.tar.gz.sha512
echo 'ec7416729f99f5eef6700507e740552e771d6dd8863f757311538d7d67a0eecd3426381bd42a7ddbf0771bdde8bba5cb943f60031ae3567d6a3dcac738facda8  mercurial-2.2.2.tar.gz' > mercurial-2.2.2.tar.gz.sha512
echo 'de3dd37f753614055dcfed910e9886e03688b8078492df3da94b1ec37be796030be93291cba09e8212fffd3e0a63b086902c3c25a996cf1439e15c5b16e014d9  numpy-1.6.1.tar.gz' > numpy-1.6.1.tar.gz.sha512
echo '5ad681f99e75849a5ca6f439c7a19bb51abc73d121b50f4f8e4c0da42891950f30407f761a53f0fe51b370b1dbd4c4f5a480557cb2444c8c7c7d5412b328a474  sqlite-autoconf-3070500.tar.gz' > sqlite-autoconf-3070500.tar.gz.sha512
echo 'edae735960279d92acf58e1f4095c6392a7c2059b8f1d2c46648fc608a0fb06b392db2d073f4973f5762c034ea66596e769b95b3d26ad963a086b9b2d09825f2  zlib-1.2.3.tar.bz2' > zlib-1.2.3.tar.bz2.sha512
echo '42021737c93cea513116e6051cff9b803e3f25d6019c74370b42f4c91d1af73e94ac2b7ace85b7565fa3d45b79231079bd48a242531beeafa33c36d7139ce838  ipython-0.13.tar.gz' > ipython-0.13.tar.gz.sha512
echo 'fb3cf421b2dc48c31956b3e3ee4ab6ebc743deec3bf626c2238a1996c8c51be87260bd6aa662793a1f0c34dcda9b3146763777bb162dfad6fec4ca7acc403b2e  zeromq-2.2.0.tar.gz' > zeromq-2.2.0.tar.gz.sha512
echo 'd761b492352841cdc125d9f0c99ee6d6c435812472ea234728b7f0fb4ad1048e1eec9b399df2081fbc926566f333f7780fedd0ce23255a6633fe5c60ed15a6af  pyzmq-2.1.11.tar.gz' > pyzmq-2.1.11.tar.gz.sha512
echo '57fa5e57dfb98154a42d2d477f29401c2260ae7ad3a8128a4098b42ee3b35c54367b1a3254bc76b9b3b14b4aab7c3e1135858f68abc5636daedf2f01f9b8a3cf  tornado-2.2.tar.gz' > tornado-2.2.tar.gz.sha512

# Individual processes
[ -z "$HDF5_DIR" ] && get_ytproject hdf5-1.8.7.tar.gz
[ $INST_ZLIB -eq 1 ] && get_ytproject zlib-1.2.3.tar.bz2 
[ $INST_BZLIB -eq 1 ] && get_ytproject bzip2-1.0.5.tar.gz
[ $INST_PNG -eq 1 ] && get_ytproject libpng-1.5.12.tar.gz
[ $INST_FTYPE -eq 1 ] && get_ytproject freetype-2.4.4.tar.gz
[ $INST_SQLITE3 -eq 1 ] && get_ytproject sqlite-autoconf-3070500.tar.gz
[ $INST_PYX -eq 1 ] && get_ytproject PyX-0.11.1.tar.gz
[ $INST_0MQ -eq 1 ] && get_ytproject zeromq-2.2.0.tar.gz
[ $INST_0MQ -eq 1 ] && get_ytproject pyzmq-2.1.11.tar.gz
[ $INST_0MQ -eq 1 ] && get_ytproject tornado-2.2.tar.gz
get_ytproject Python-2.7.2.tgz
get_ytproject numpy-1.6.1.tar.gz
get_ytproject matplotlib-1.1.0.tar.gz
get_ytproject mercurial-2.2.2.tar.gz
get_ytproject ipython-0.13.tar.gz
get_ytproject h5py-2.0.1.tar.gz
get_ytproject Cython-0.16.tar.gz
get_ytproject reason-js-20120623.zip

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
    if [ ! -e hdf5-1.8.7/done ]
    then
        [ ! -e hdf5-1.8.7 ] && tar xfz hdf5-1.8.7.tar.gz
        echo "Installing HDF5"
        cd hdf5-1.8.7
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

if [ ! -e Python-2.7.2/done ]
then
    echo "Installing Python.  This may take a while, but don't worry.  YT loves you."
    [ ! -e Python-2.7.2 ] && tar xfz Python-2.7.2.tgz
    cd Python-2.7.2
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
    do_setup_py mercurial-2.2.2
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
        ( ${HG_EXEC} --debug clone http://hg.yt-project.org/yt-supplemental/ 2>&1 ) 1>> ${LOG_FILE}
        # Recently the hg server has had some issues with timeouts.  In lieu of
        # a new webserver, we are now moving to a three-stage process.
        # First we clone the repo, but only up to r0.
        ( ${HG_EXEC} --debug clone http://hg.yt-project.org/yt/ ./yt-hg 2>&1 ) 1>> ${LOG_FILE}
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

do_setup_py numpy-1.6.1 ${NUMPY_ARGS}

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
mkdir -p ${DEST_DIR}/src/matplotlib-1.1.0
echo "[directories]" >> ${DEST_DIR}/src/matplotlib-1.1.0/setup.cfg
echo "basedirlist = ${DEST_DIR}" >> ${DEST_DIR}/src/matplotlib-1.1.0/setup.cfg
do_setup_py matplotlib-1.1.0
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

do_setup_py ipython-0.13
do_setup_py h5py-2.0.1
do_setup_py Cython-0.16
[ $INST_PYX -eq 1 ] && do_setup_py PyX-0.11.1

echo "Doing yt update, wiping local changes and updating to branch ${BRANCH}"
MY_PWD=`pwd`
cd $YT_DIR
( ${HG_EXEC} pull && ${HG_EXEC} up -C ${BRANCH} 2>&1 ) 1>> ${LOG_FILE}

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
    ${HG_EXEC} clone https://enzo.googlecode.com/hg/ ./enzo-hg-stable
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
    echo "    (yt)$ "
    echo
    echo "This modifies the environment variables YT_DEST, PATH, PYTHONPATH, and"
    echo "LD_LIBRARY_PATH to match your new yt install. But don't worry - as soon"
    echo "as you are done you can run 'deactivate' to return to your previous"
    echo "shell environment.  If you use csh, just append .csh to the above."
    echo
    echo "For interactive data analysis and visualization, we recommend running"
    echo "the IPython interface, which will become more fully featured with time:"
    echo
    echo "    $DEST_DIR/bin/iyt"
    echo
    echo "For command line analysis run:"
    echo
    echo "    $DEST_DIR/bin/yt"
    echo
    echo "To bootstrap a development environment for yt, run:"
    echo 
    echo "    $DEST_DIR/bin/yt bootstrap_dev"
    echo
    echo "Note of interest: this installation will use the directory:"
    echo "    $YT_DIR"
    echo "as the source for all the YT code.  This means you probably shouldn't"
    echo "delete it, but on the plus side, any changes you make there are"
    echo "automatically propagated."
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
    echo "    http://lists.spacepope.org/listinfo.cgi/yt-users-spacepope.org"
    echo
    echo "========================================================================"
    echo
    echo "Oh, look at me, still talking when there's science to do!"
    echo "Good luck, and email the user list if you run into any problems."
}

print_afterword
print_afterword >> ${LOG_FILE}
