from yt.loaders import load_sample2 as load_sample


def test_load_small_sample_files():
    # todo:
    # - mock the test data dir to force download
    # - inspect logger history to check that second call doesn't require download
    load_sample("ToroShockTube")
