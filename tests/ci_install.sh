set -x   # Show which command is being run

# Step 1: install system level dependencies for ratarmount
if [[ ${sync_args} == *"full"* ]]; then
    case ${RUNNER_OS} in
    linux|Linux)
        sudo apt-get -qqy update
        sudo apt-get -qqy install libfuse2
        ;;
    osx|macOS)
        sudo mkdir -p /usr/local/man
        sudo chown -R "${USER}:admin" /usr/local/man
        HOMEBREW_NO_AUTO_UPDATE=1 brew install macfuse
        ;;
    esac
fi


# TEMPORARY STEP 2 for ewah-bool-utils  testing. ROLLBACK BEFORE MERGE.
# manually install build requirements and use no-build-isolation so that yt picks up the ewah-bool-utils
# build during its own build.

# need a venv first (python version should be picked up from the astral-sh/setup-uv python config)
uv venv
# build requirements for both yt and ewah-bool-utils
uv pip install meson-python meson ninja "cython>=3.0.3" "numpy>=2" "setuptools>=77" wheel
# build ewah-bool-utils dev branch (--no-config to skip the no-build-package list)
uv pip install --no-config --no-build-isolation --reinstall-package ewah-bool-utils "ewah-bool-utils @ git+https://github.com/chrishavlin/ewah_bool_utils.git@opaque_classes"

# now build yt (will use the dev ewah-bool-utils) and manually pick up the extras.
extras="test"
groups=""
[[ ${sync_args} == *"full"* ]] && extras="test,full"
[[ ${sync_args} == *"nosetest"* ]] && groups="--group nosetest"
[[ ${sync_args} == *"mapping"* ]] && groups="--group mapping"
uv pip install --no-build-isolation ${groups} -e ".[${extras}]"

# Step 2: install deps and yt
# installing in editable mode so this script may be used locally by developers
# but the primary intention is to embed this script in CI jobs
# uv sync --extra=test ${sync_args}

# Disable excessive output
uv run --no-sync yt config set --local yt log_level 50
cat yt.toml

set +x
