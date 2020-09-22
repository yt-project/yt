set -x   # Show which command is being run

case $TRAVIS_OS_NAME in
linux)
    $PIP install --upgrade virtualenv
    python -m virtualenv venv
    source venv/bin/activate
    export PATH=/usr/lib/ccache:$PATH
    ;;
windows)
    # From https://github.com/trichter/conda4travis/blob/master/conda4travis.sh
    export MINICONDA=/c/miniconda
    MINICONDA_WIN=$(cygpath --windows $MINICONDA)
    choco install openssl.light
    choco install miniconda3 --params="'/AddToPath:0 /D:$MINICONDA_WIN'"
    # The following lines remove bin directories from chocolatey.
    # Otherwise f2py might find the pre-installed fortran compiler
    # instead of the conda fortran compiler.
    # See dicussion in https://github.com/trichter/toeplitz/pull/2
    echo original PATH $PATH
    PATH=$(echo "$PATH" | sed -e 's|:/c/ProgramData/chocolatey/bin||')
    PATH=$(echo "$PATH" | sed -e 's|:/c/ProgramData/chocolatey/lib/mingw/tools/install/mingw64/bin||')
    echo manipulated PATH $PATH
    # the following line is necessary since conda 4.7
    # see travis fail https://travis-ci.org/trichter/conda4travis/jobs/592665691
    # see fix https://github.com/conda/conda/issues/8836#issuecomment-506388019
    source $MINICONDA/Scripts/activate

    source $MINICONDA/etc/profile.d/conda.sh
    hash -r
    conda config --set always_yes yes
    conda info -a
    ;;
osx)
    sudo mkdir -p /usr/local/man
    sudo chown -R "${USER}:admin" /usr/local/man
    HOMEBREW_NO_AUTO_UPDATE=1 brew install hdf5 proj geos open-mpi netcdf ccache
    HOMEBREW_NO_AUTO_UPDATE=1 brew uninstall gdal postgis numpy  # WHY?
    export PATH=/usr/local/opt/ccache/libexec:$PATH
    ;;
esac
mkdir -p $HOME/.config/yt
echo "[yt]" > $HOME/.config/yt/ytrc
echo "suppressStreamLogging = True" >> $HOME/.config/yt/ytrc
cat $HOME/.config/yt/ytrc
cp tests/matplotlibrc .

set +x  # Do not show which command is being run (Travis default)