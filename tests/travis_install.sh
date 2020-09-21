set -x   # Show which command is being run

echo "Travis build stage: $TRAVIS_BUILD_STAGE_NAME"

# Step 1: pre-install required packages
if [[ "$TRAVIS_OS_NAME" == "windows" ]] && [[ $MINIMAL != 1 ]]; then
    # Install some dependencies using conda (if not doing a minimal run)
    CYTHON=$(grep cython tests/test_prerequirements.txt)
    NUMPY=$(grep numpy tests/test_prerequirements.txt)

    CARTOPY=$(grep cartopy tests/test_requirements.txt)
    H5PY=$(grep h5py tests/test_requirements.txt)
    MATPLOTLIB=$(grep matplotlib tests/test_requirements.txt)
    SCIPY=$(grep scipy tests/test_requirements.txt)
    conda install --yes -c conda-forge $CYTHON $NUMPY $CARTOPY $H5PY $MATPLOTLIB $SCIPY
else
    ccache -s
    # Upgrade pip and setuptools and wheel to get clean install
    $PIP install --upgrade pip
    $PIP install --upgrade wheel
    $PIP install --upgrade setuptools
fi

# Step 2: install required packages (depending on whether the build is minimal)
if [[ $MINIMAL == 1 ]]; then
    # Ensure numpy and cython are installed so dependencies that need to be built
    # don't error out
    # The first numpy to support py3.6 is 1.12, but numpy 1.13 matches
    # unyt so we'll match it here.
    $PIP install numpy==1.13.3 cython==0.26.1
    $PIP install -r tests/test_minimal_requirements.txt
else
    # Getting cartopy installed requires getting cython and numpy installed
    # first; this is potentially going to be fixed with the inclusion of
    # pyproject.toml in cartopy.
    # These versions are pinned, so we will need to update/remove them when
    # the hack is no longer necessary.
    $PIP install -r tests/test_prerequirements.txt
    CFLAGS="$CFLAGS -DACCEPT_USE_OF_DEPRECATED_PROJ_API_H" $PIP install -r tests/test_requirements.txt
fi

# Step 3: install yt
$PIP install -e .

set +x  # Do not show which command is being run (Travis default)