set -x   # Show which command is being run

# Step 1: pre-install required packages
if [[ ${sync_args} == *"full"* || ${sync_args} == *"mapping"* ]]; then
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

# Step 2: install deps and yt
# installing in editable mode so this script may be used locally by developers
# but the primary intention is to embed this script in CI jobs
uv sync --extra=test ${sync_args}

# Disable excessive output
uv run --no-sync yt config set --local yt log_level 50
cat yt.toml

set +x
