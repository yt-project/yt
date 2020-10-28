import os

import pytest

from yt.config import ytcfg

if __name__ == "__main__":
    pytest_args = [
        f"--local-dir={os.path.join(ytcfg.get('yt', 'test_data_dir'), 'answers')}",
        "-c=pytest_answer.ini",
        "--junitxml=unittests.xml",
        "--answer-big-data",
        f"-n {int(os.environ.get('NUM_WORKERS', 6))}",
    ]
    pytest.main(pytest_args)
