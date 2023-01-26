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
echo "[yt]\nlog_level = 50" > yt.toml
cat yt.toml

# Sets default backend to Agg
cp tests/matplotlibrc .

# Step 1: pre-install required packages
if [[ "${RUNNER_OS}" == "Windows" ]] && [[ ${dependencies} != "minimal" ]]; then
    # windows_conda_requirements.txt is a survivance of test_requirements.txt
    # keep in sync: setup.cfg
    while read requirement; do conda install --yes $requirement; done < tests/windows_conda_requirements.txt
else
    # upgrading pip to guarantee installing extra dependencies with [full] is supported
    # this is only necessary for some versions of Python 3.8 and 3.9
    # see https://github.com/yt-project/yt/issues/4270
    python -m pip install --upgrade pip
fi

# Step 2: install deps and yt
# installing in editable mode so this script may be used locally by developers
# but the primary intention is to embed this script in CI jobs
if [[ ${dependencies} == "minimal" ]]; then
    # test with minimal versions of runtime dependencies
    python -m pip install -e ".[test,minimal]"
elif [[ ${dependencies} == "full" ]]; then
    # test with all optional runtime dependencies

    # this is required for cartopy. see
    # https://scitools.org.uk/cartopy/docs/latest/installing.html?highlight=install#building-from-source
    python -m pip install shapely --no-binary=shapely
    python -m pip install -e ".[test,full]"
else
   # test with no special requirements
   python -m pip install -e ".[test]"
fi

set +x
