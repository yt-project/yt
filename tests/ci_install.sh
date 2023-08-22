set -x   # Show which command is being run

# Sets default backend to Agg
cp tests/matplotlibrc .

# Step 1: pre-install required packages
if [[ ${dependencies} == "full" || ${dependencies} == "cartopy" ]]; then
    # upgrading pip to guarantee installing extra dependencies with [full] is supported
    # this is only necessary for some versions of Python 3.8 and 3.9
    # see https://github.com/yt-project/yt/issues/4270
    python -m pip install 'pip>=21.2'

    case ${RUNNER_OS} in
    linux|Linux)
        sudo apt-get -qqy update
        sudo apt-get -qqy install \
        libhdf5-serial-dev \
        libnetcdf-dev \
        libopenmpi-dev \
        libfuse2
        ;;
    osx|macOS)
        sudo mkdir -p /usr/local/man
        sudo chown -R "${USER}:admin" /usr/local/man
        HOMEBREW_NO_AUTO_UPDATE=1 brew install hdf5 open-mpi netcdf ccache macfuse
        ;;
    esac
fi

# Step 2: install deps and yt
# installing in editable mode so this script may be used locally by developers
# but the primary intention is to embed this script in CI jobs
if [[ ${dependencies} == "minimal" ]]; then
    # test with minimal versions of runtime dependencies
    python -m pip install -e ".[test,minimal]"
elif [[ ${dependencies} == "cartopy" ]]; then
    python -m pip install 'cartopy>=0.22'
    # scipy is an optional dependency to cartopy
    python -m pip install scipy
    python -m pip install -e ".[test]"
elif [[ ${dependencies} == "full" ]]; then
    # test with all optional runtime dependencies
    python -m pip install -e ".[test,full]"
else
   # test with no special requirements
   python -m pip install -e ".[test]"
fi

# Disable excessive output
yt config set --local yt log_level 50
cat yt.toml

set +x
