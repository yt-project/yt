import os

import yt


def test_save_to_path(tmp_path):
    ds = yt.testing.fake_amr_ds()
    p = yt.SlicePlot(ds, "z", "Density")
    p.save(f"{tmp_path}/")
    assert len(list((tmp_path).glob("*.png"))) == 1


def test_save_to_missing_path(tmp_path):
    # the missing layer should be created
    ds = yt.testing.fake_amr_ds()
    p = yt.SlicePlot(ds, "z", "Density")

    # using forward slashes should work even on windows !
    save_path = os.path.join(tmp_path / "out") + "/"
    p.save(save_path)
    assert os.path.exists(save_path)
    assert len(list((tmp_path / "out").glob("*.png"))) == 1


def test_save_to_missing_path_with_file_prefix(tmp_path):
    # see issue
    # https://github.com/yt-project/yt/issues/3210
    ds = yt.testing.fake_amr_ds()
    p = yt.SlicePlot(ds, "z", "Density")
    p.save(tmp_path.joinpath("out", "saymyname"))
    assert (tmp_path / "out").exists()
    output_files = list((tmp_path / "out").glob("*.png"))
    assert len(output_files) == 1
    assert output_files[0].stem.startswith("saymyname")  # you're goddamn right
