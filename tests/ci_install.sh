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
    # windows_conda_requirements.txt is a survivance of test_requirements.txt
    # keep in sync: setup.cfg
    while read requirement; do conda install --yes $requirement; done < tests/windows_conda_requirements.txt
else
    python -m pip install --upgrade pip
    python -m pip install --upgrade wheel
    python -m pip install --upgrade setuptools
fi

# Step 2: install deps and yt
if [[ ${dependencies} == "minimal" ]]; then
    python -m pip install -e .[test,minimal]
else
    # Cython and numpy are build-time requirements to the following optional deps in yt
    # - cartopy
    # - netcdf4
    # - pyqt5
    # The build system is however not specified properly in these projects at the moment
    # which means we have to install the build-time requirements first.
    # It is possible that these problems will be fixed in the future if upstream projects
    # include a pyproject.toml file or use any pip-comptatible solution to remedy this.
    python -m pip install numpy>=1.19.4 cython~=0.29.21

    # this is required for cartopy. It should normally be specified in our setup.cfg as
    # cartopy[plotting]
    # However it doesn't work on Ubuntu 18.04 (used in CI at the time of writing)
    python -m pip install shapely --no-binary=shapely
    CFLAGS="$CFLAGS -DACCEPT_USE_OF_DEPRECATED_PROJ_API_H" python -m pip install -e .[test,full]
fi

set +x
