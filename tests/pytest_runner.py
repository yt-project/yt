"""This is a helper script for running answer tests on CI services.

It's currently used on:
  * Jenkins
  * GHA
for executing answer tests and optionally generating new answers.
"""


import glob
import os

import pytest

if __name__ == "__main__":
    os.environ["OMP_NUM_THREADS"] = "1"
    pytest_args = [
        "-s",
        "-v",
        "-rsfE",  # it means -r "sfE" (show skipped, failed, errors), no -r -s -f -E
        "--with-answer-testing",
        "-m answer_test",
        f"-n {int(os.environ.get('NUM_WORKERS', 1))}",
        "--dist=loadscope",
    ]
    pytest.main(pytest_args + ["--local-dir=answer-store", "--junitxml=answers.xml"])

    if files := glob.glob("generate_test*.txt"):
        tests = set()
        for fname in files:
            with open(fname) as fp:
                tests |= set(fp.read().splitlines())
        output_dir = "artifacts"
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        pytest.main(
            pytest_args + [f"--local-dir={output_dir}", "--answer-store"] + list(tests)
        )
