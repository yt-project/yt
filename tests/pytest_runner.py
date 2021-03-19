import os

import pytest

from yt.config import ytcfg

if __name__ == "__main__":
    os.environ["OMP_NUM_THREADS"] = "1"
    pytest_args = [
        f"--local-dir={os.path.join(ytcfg.get('yt', 'test_data_dir'), 'answers')}",
        "-c=pytest_answer.ini",
        "--junitxml=unittests.xml",
        "--answer-big-data",
        f"-n {int(os.environ.get('NUM_WORKERS', 6))}",
        "--dist=loadscope",
    ]
    pytest.main(pytest_args)
