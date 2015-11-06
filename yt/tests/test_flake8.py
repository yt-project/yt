import subprocess
import yt
import os

from yt.testing import requires_module


@requires_module('flake8')
def test_flake8():
    yt_dir = os.path.dirname(os.path.abspath(yt.__file__))
    initial_dir = os.getcwd()
    os.chdir(yt_dir)
    output_file = os.path.sep.join([os.path.dirname(initial_dir), 'flake8.out'])
    if os.path.exists(output_file):
        os.remove(output_file)
    output_string = "--output-file=%s" % output_file
    subprocess.call(['flake8', output_string, os.curdir])
    os.chdir(initial_dir)
    with open(output_file) as f:
        flake8_output = f.readlines()
    if flake8_output != []:
        raise AssertionError(
            "flake8 found style errors:\n\n%s" % "\n".join(flake8_output))
