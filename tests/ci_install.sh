set -x   # Show which command is being run

case ${RUNNER_OS} in
linux|Linux)
    sudo apt-get -qqy update
    sudo apt-get -qqy install \
      libhdf5-serial-dev \
      libnetcdf-dev \
      libproj-dev \
      proj-data \
      proj-bin \
      libgeos-dev \
      libopenmpi-dev
    ;;
osx|macOS)
    sudo mkdir -p /usr/local/man
    sudo chown -R "${USER}:admin" /usr/local/man
    brew update
    # proj can be unpinned when upstream incompatibility issue is resolved. See
    # https://github.com/SciTools/cartopy/issues/1140
    export LDFLAGS="$LDFLAGS -L/usr/local/opt/proj@7/lib"
    export CPPFLAGS="$CPPFLAGS -I/usr/local/opt/proj@7/include"
    export CFLAGS="$CFLAGS -I/usr/local/opt/proj@7/include"
    export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:/usr/local/opt/proj@7/lib/pkgconfig"
    HOMEBREW_NO_AUTO_UPDATE=1 brew install hdf5 proj@7 geos open-mpi netcdf ccache
    ;;
esac

# Disable excessive output
mkdir -p $HOME/.config/yt
echo "[yt]" > $HOME/.config/yt/yt.toml
echo "suppress_stream_logging = true" >> $HOME/.config/yt/yt.toml
cat $HOME/.config/yt/yt.toml
# Sets default backend to Agg
cp tests/matplotlibrc .

# Step 1: pre-install required packages
if [[ "${RUNNER_OS}" == "Windows" ]] && [[ ${dependencies} != "minimal" ]]; then
    # Install some dependencies using conda (if not doing a minimal run)
    CYTHON=$(grep cython tests/test_prerequirements.txt)
    NUMPY=$(grep numpy tests/test_prerequirements.txt)

    CARTOPY=$(grep cartopy tests/test_requirements.txt)
    H5PY=$(grep h5py tests/test_requirements.txt)
    MATPLOTLIB=$(grep matplotlib tests/test_requirements.txt)
    SCIPY=$(grep scipy tests/test_requirements.txt)
    conda config --set always_yes yes
    conda info -a
    conda install --quiet --yes -c conda-forge \
      $CYTHON $NUMPY $CARTOPY $H5PY $MATPLOTLIB $SCIPY
else
    python -m pip install --upgrade pip
    python -m pip install --upgrade wheel
    python -m pip install --upgrade setuptools
fi

# Step 2: install required packages (depending on whether the build is minimal)
if [[ ${dependencies} == "minimal" ]]; then
    # Ensure numpy and cython are installed so dependencies that need to be built
    # don't error out
    # The first numpy to support py3.6 is 1.12, but numpy 1.13 matches
    # unyt so we'll match it here.
    python -m pip install numpy==1.13.3 cython==0.26.1
    python -m pip install -r tests/test_minimal_requirements.txt
else
    # Getting cartopy installed requires getting cython and numpy installed
    # first; this is potentially going to be fixed with the inclusion of
    # pyproject.toml in cartopy.
    # These versions are pinned, so we will need to update/remove them when
    # the hack is no longer necessary.
    python -m pip install -r tests/test_prerequirements.txt
    CFLAGS="$CFLAGS -DACCEPT_USE_OF_DEPRECATED_PROJ_API_H" python -m pip install -r tests/test_requirements.txt
fi

# Step 3: install yt
python -m pip install -e .

set +x
