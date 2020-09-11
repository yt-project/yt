from yt.testing import assert_equal
from yt.utilities.sample_data import PoochHandle

names = {
    "t1": {
        "load_name": "IsolatedGalaxy.tar.gz",
        "answers": {
            "fileext": "IsolatedGalaxy.tar.gz",
            "basename": "IsolatedGalaxy",
            "extension": "tar",
        },
    },
    "t2": {
        "load_name": "IsolatedGalaxy",
        "answers": {
            "fileext": "IsolatedGalaxy.tar.gz",
            "basename": "IsolatedGalaxy",
            "extension": "tar",
        },
    },
    "t3": {
        "load_name": "apec_emissivity_v3.h5",
        "answers": {
            "fileext": "apec_emissivity_v3.h5",
            "basename": "apec_emissivity_v3",
            "extension": "h5",
        },
    },
    "t4": {
        "load_name": "apec_emissivity_v3.hdf5",
        "answers": {
            "fileext": "apec_emissivity_v3.hdf5",
            "basename": "apec_emissivity_v3",
            "extension": "h5",
        },
    },
    "t5": {
        "load_name": "solution-00027.0000.vtu",
        "answers": {
            "fileext": "solution-00027.0000.vtu",
            "basename": "solution-00027.0000",
            "extension": ".vtu",
        },
    },
}


def test_name_validator():
    for test in names:
        fileext, bname, ext = PoochHandle._validate_sample_fname(
            names[test]["load_name"]
        )
        expected_answers = names[test]["answers"]
        assert_equal(fileext, expected_answers["fileext"])
        assert_equal(bname, expected_answers["basename"])
        assert_equal(ext, expected_answers["extension"])
