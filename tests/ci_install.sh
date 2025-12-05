set -x   # Show which command is being run

# Sets default backend to Agg
cp tests/matplotlibrc .

# Step 1: pre-install required packages
if [[ ${sync_args} == *"full"* || ${sync_args} == *"mapping"* ]]; then
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
        # uninstalling pkg-config to workaround a bug in macOS image
        # https://github.com/Homebrew/homebrew-core/pull/198691#issuecomment-2495500991
        # this can be cleaned-up once the following patch is released:
        # https://github.com/actions/runner-images/pull/11015
        HOMEBREW_NO_AUTO_UPDATE=1 brew uninstall pkg-config@0.29.2 || true
        HOMEBREW_NO_AUTO_UPDATE=1 brew install hdf5 open-mpi netcdf ccache macfuse
        ;;
    esac
fi

# Step 2: install deps and yt
# installing in editable mode so this script may be used locally by developers
# but the primary intention is to embed this script in CI jobs
uv sync --extra=test ${sync_args}

# Disable excessive output
uv run --no-sync yt config set --local yt log_level 50
cat yt.toml

set +x
