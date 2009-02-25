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

DEST_DIR="`pwd`/yt-`uname -p`"   # Installation location
PY_DIR="/Library/Frameworks/Python.framework/Versions/Current/"

# Here's where you put the HDF5 path if you like; otherwise it'll download it
# and install it on its own
#HDF5_DIR=

# If you've got YT some other place, set this to point to it.
YT_DIR=""

#------------------------------------------------------------------------------#
#                                                                              #
# Okay, the script starts here.  Feel free to play with it, but hopefully      #
# it'll work as is.                                                            #
#                                                                              #
#------------------------------------------------------------------------------#

function do_exit
{
    echo "Failure.  Check ${LOG_FILE}."
    exit 1
}

function do_setup_py
{
    [ -e $1/done ] && return
    echo "Installing $1"
    echo "Please type superuser password"
    [ ! -e $1 ] && tar xfz $1.tar.gz
    cd $1
    shift
    ( ${PY_DIR}/bin/python2.5 setup.py install $* 2>&1 ) 1>> ${LOG_FILE} || do_exit
    touch done
    cd ..
}

function get_enzotools
{
    echo "Downloading $1 from yt.enzotools.org"
    [ -e $1 ] && return
    wget -nv "http://yt.enzotools.org/dependencies/osx/$1" || do_exit
    wget -nv "http://yt.enzotools.org/dependencies/osx/$1.md5" || do_exit
    ( which md5sum &> /dev/null ) || return # return if we don't have md5sum
    ( md5sum -c $1.md5 2>&1 ) 1>> ${LOG_FILE} || do_exit
}

function self_install
{
    echo "Installing ${1}.  You will need to handle this procedure."
    echo
    echo "Press enter to start, then return when finished."
    echo
    read LOKI
    open ${1}
    echo
    echo "Press enter when the installation is complete."
    read LOKI
}

ORIG_PWD=`pwd`

LOG_FILE="${DEST_DIR}/yt_install.log"

if [ -z "${DEST_DIR}" ]
then
    echo "Edit this script, set the DEST_DIR parameter and re-run."
    exit 1
fi

echo "Installing into ${DEST_DIR}"

mkdir -p ${DEST_DIR}/src
cd ${DEST_DIR}/src

# Individual processes
if [ -z "$HDF5_DIR" ]
then
    echo "Downloading HDF5"
    get_enzotools hdf5-1.6.8.tar.gz
fi

get_enzotools python-2.5.4-macosx.dmg
get_enzotools wxPython2.8-osx-unicode-2.8.9.2-universal-py2.5.dmg
get_enzotools numpy-1.2.1-py2.5-macosx10.5.dmg
get_enzotools matplotlib-0.98.5.2-py2.5-mpkg.zip
get_enzotools ipython-0.9.1.tar.gz
get_enzotools tables-2.1.tar.gz

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
    HDF5_DIR=${DEST_DIR}
fi

[ ! -e ${DEST_DIR}/src/py_done ] && self_install \
    python-2.5.4-macosx.dmg
touch ${DEST_DIR}/src/py_done

[ ! -e ${DEST_DIR}/src/wx_done ] && self_install \
    wxPython2.8-osx-unicode-2.8.9.2-universal-py2.5.dmg
touch ${DEST_DIR}/src/wx_done

echo "Installing setuptools"
( ${PY_DIR}/bin/python2.5 ${YT_DIR}/ez_setup.py 2>&1 ) 1>> ${LOG_FILE} || do_exit

[ ! -e ${DEST_DIR}/src/np_done ] && self_install \
    numpy-1.2.1-py2.5-macosx10.5.dmg
touch ${DEST_DIR}/src/np_done

if [ ! -e ${DEST_DIR}/src/mp_done ]
then
    [ ! -e matplotlib-0.98.5.2-py2.5-mpkg ] && unzip matplotlib-0.98.5.2-py2.5-mpkg.zip
    self_install matplotlib-0.98.5.2-py2.5
    touch ${DEST_DIR}/src/mp_done
fi

do_setup_py ipython-0.9.1
do_setup_py tables-2.1 --hdf5=${HDF5_DIR}

echo "Doing yt update"
MY_PWD=`pwd`
cd $YT_DIR
( svn up 2>&1 ) 1>> ${LOG_FILE}

echo "Installing yt"
echo $HDF5_DIR > hdf5.cfg
( ${PY_DIR}/bin/python2.5 setup.py develop 2>&1 ) 1>> ${LOG_FILE} || do_exit
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
echo "PATH            => $PY_DIR/bin/"
echo
echo "For interactive data analysis and visualization, we recommend running"
echo "the IPython interface, which will become more fully featured with time:"
echo
echo "$PY_DIR/bin/iyt"
echo
echo "For command line analysis run:"
echo
echo "$PY_DIR/bin/yt"
echo
echo "Note of interest: this installation will use the directory"
echo "$YT_DIR"
echo "as the source for all the YT code.  This means you probably shouldn't"
echo "delete it, but on the plus side, any changes you make there are"
echo "automatically propagated."
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
