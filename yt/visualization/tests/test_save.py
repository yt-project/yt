import os
import re

import pytest
from PIL import Image

import yt
from yt.testing import fake_amr_ds

# avoid testing every supported format as some backends may be buggy on testing platforms
FORMATS_TO_TEST = [".eps", ".pdf", ".svg"]


@pytest.fixture(scope="session")
def simple_sliceplot():
    ds = fake_amr_ds()
    p = yt.SlicePlot(ds, "z", "Density")
    yield p


def test_save_to_path(simple_sliceplot, tmp_path):
    p = simple_sliceplot
    p.save(f"{tmp_path}/")
    assert len(list((tmp_path).glob("*.png"))) == 1


def test_metadata(simple_sliceplot, tmp_path):
    simple_sliceplot.save(tmp_path / "ala.png")
    with Image.open(tmp_path / "ala.png") as img:
        assert "Software" in img.info
        assert "yt-" in img.info["Software"]
    simple_sliceplot.save(tmp_path / "ala.pdf")
    with open(tmp_path / "ala.pdf", "rb") as f:
        assert b"|yt-" in f.read()


def test_save_to_missing_path(simple_sliceplot, tmp_path):
    # the missing layer should be created
    p = simple_sliceplot

    # using forward slashes should work even on windows !
    save_path = os.path.join(tmp_path / "out") + "/"
    p.save(save_path)
    assert os.path.exists(save_path)
    assert len(list((tmp_path / "out").glob("*.png"))) == 1


def test_save_to_missing_path_with_file_prefix(simple_sliceplot, tmp_path):
    # see issue
    # https://github.com/yt-project/yt/issues/3210
    p = simple_sliceplot
    p.save(tmp_path.joinpath("out", "saymyname"))
    assert (tmp_path / "out").exists()
    output_files = list((tmp_path / "out").glob("*.png"))
    assert len(output_files) == 1
    assert output_files[0].stem.startswith("saymyname")  # you're goddamn right


@pytest.mark.parametrize("ext", FORMATS_TO_TEST)
def test_suffix_from_filename(ext, simple_sliceplot, tmp_path):
    p = simple_sliceplot

    target = (tmp_path / "myfile").with_suffix(ext)
    # this shouldn't raise a warning, see issue
    # https://github.com/yt-project/yt/issues/3667
    p.save(target)
    assert target.is_file()


@pytest.mark.parametrize("ext", FORMATS_TO_TEST)
def test_suffix_clashing(ext, simple_sliceplot, tmp_path):
    if ext == ".png":
        pytest.skip()

    p = simple_sliceplot

    target = (tmp_path / "myfile").with_suffix(ext)
    expected_warning = re.compile(
        rf"Received two valid image formats {ext.removeprefix('.')!r} "
        r"\(from filename\) and 'png' \(from suffix\)\. The former is ignored\."
    )

    with pytest.warns(UserWarning, match=expected_warning):
        p.save(target, suffix="png")
    output_files = list(tmp_path.glob("*.png"))
    assert len(output_files) == 1
    assert output_files[0].stem.startswith("myfile")
    assert not list((tmp_path / "out").glob(f"*.{ext}"))


def test_invalid_format_from_filename(simple_sliceplot, tmp_path):
    p = simple_sliceplot
    target = (tmp_path / "myfile").with_suffix(".nope")

    p.save(target)
    output_files = list(tmp_path.glob("*"))
    assert len(output_files) == 1
    # the output filename may contain a generated part
    # it's not exactly clear if it's desirable or intended in this case
    # so we just check conditions that should hold in any case
    assert output_files[0].name.startswith("myfile.nope")
    assert output_files[0].name.endswith(".png")


def test_invalid_format_from_suffix(simple_sliceplot, tmp_path):
    p = simple_sliceplot
    target = tmp_path / "myfile"
    with pytest.raises(ValueError, match=r"Unsupported file format 'nope'"):
        p.save(target, suffix="nope")
