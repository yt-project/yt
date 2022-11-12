set -x   # Show which command is being run

# geos is needed for cartopy, see
# https://scitools.org.uk/cartopy/docs/latest/installing.html?highlight=install#building-from-source
case ${RUNNER_OS} in
linux|Linux)
    sudo apt-get -qqy update
    sudo apt-get -qqy install \
      libhdf5-serial-dev \
      libnetcdf-dev \
      libgeos-dev \
      libopenmpi-dev \
      libfuse2
    ;;
osx|macOS)
    sudo mkdir -p /usr/local/man
    sudo chown -R "${USER}:admin" /usr/local/man
    HOMEBREW_NO_AUTO_UPDATE=1 brew install hdf5 geos open-mpi netcdf ccache macfuse
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
# installing in editable mode so this script may be used locally by developers
# but the primary intention is to embed this script in CI jobs
if [[ ${dependencies} == "minimal" ]]; then
    # test with minimal versions of runtime dependencies
    python -m pip install -e .[test,minimal]
elif [[ ${dependencies} == "full" ]]; then
    # test with all optional runtime dependencies

    # this is required for cartopy. see
    # https://scitools.org.uk/cartopy/docs/latest/installing.html?highlight=install#building-from-source
    python -m pip install shapely --no-binary=shapely
    CFLAGS="$CFLAGS -DACCEPT_USE_OF_DEPRECATED_PROJ_API_H" python -m pip install -e .[test,full]
else
   # test with no special requirements
   python -m pip install -e .[test]
fi

set +x
