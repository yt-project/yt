import pytest

from yt.config import ytcfg_defaults
from yt.utilities.command_line import run_main


# todo: maybe set this fixture to auto use with a module scope
@pytest.fixture
def tmp_config_file(tmp_path, monkeypatch):
    config_file = tmp_path / "yt.toml"
    config_file.touch()
    monkeypatch.setattr("yt.config.YTConfig.get_local_config_file", lambda: config_file)
    monkeypatch.setattr("yt.config.YTConfig.get_global_config_file", lambda: "")
    return config_file


CONFIG_SET_HELP_MESSAGE = """
usage: yt config set [-h] [--local | --global] option value

set a config value

positional arguments:
  option      The option to set.
  value       The value to set the option to.

optional arguments:
  -h, --help  show this help message and exit
  --local     Store the configuration in the global configuration file.
  --global    Store the configuration in the global configuration file.
"""


@pytest.mark.parametrize("key", ["test_data_dir", "unknown_param"])
def test_set_missing_value(tmp_config_file, capsys, key):
    retval = run_main(["config", "set", key])
    out, err = capsys.readouterr()
    assert type(retval) is int
    assert retval != 0
    assert out == ""
    # skipping the first line because it differs between macos and linux
    assert err.splitlines()[1:] == CONFIG_SET_HELP_MESSAGE.splitlines()[2:]
    assert tmp_config_file.read_text() == ""


@pytest.mark.skip(reason="no time for this now")
def test_set_wrong_type_value(tmp_config_file, capsys):
    retval = run_main(["config", "set", "logging.level", "10"])
    out, err = capsys.readouterr()

    assert type(retval) is int
    assert retval != 0
    assert out == ""
    err = err.splitlines()
    assert "INFO" in err[0]
    assert str(tmp_config_file) in err[0]
    assert (
        err[1]
        == "Error: tried to assign a value with type `<class 'int'>`' to `yt.logging.level`. Excepted type `<class 'str'>`."
    )
    assert err[2] == "This entry was last set in defaults."
    assert len(err) == 3
    assert tmp_config_file.read_text() == ""


# @pytest.mark.skip(reason="still in the shop")
def test_set_unknown_param_with_value(tmp_config_file, capsys):
    retval = run_main(["config", "set", "log_level", "0"])
    out, err = capsys.readouterr()
    assert type(retval) is int
    assert retval != 0
    assert out == ""

    err = err.splitlines()
    assert "INFO" in err[0]
    assert str(tmp_config_file) in err[0]
    assert err[1] == "Error: unknown parameter `log_level`."
    assert len(err) == 2
    assert tmp_config_file.read_text() == ""


def test_set_nested_param(tmp_config_file, capsys):
    assert tmp_config_file.read_text() == ""
    retval = run_main(["config", "set", "logging.level", "CRITICAL"])
    out, err = capsys.readouterr()

    assert retval == 0  # this fails, I need to fix it
    assert out == ""
    err = err.splitlines()
    assert "INFO" in err[0]
    assert str(tmp_config_file) in err[0]
    assert len(err) == 1
    assert tmp_config_file.read_text() == '[yt.logging]\nlevel = "CRITICAL"\n'


# Getter tests
@pytest.mark.parametrize("key", ("test_data_dir", "logging.level"))
def test_get_default_params(tmp_config_file, capsys, key):
    assert tmp_config_file.read_text() == ""
    retval = run_main(["config", "get", key])
    out, err = capsys.readouterr()

    assert retval == 0  # this fails, I need to fix it
    err = err.splitlines()
    assert "INFO" in err[0]
    assert str(tmp_config_file) in err[0]
    assert (
        err[1]
        == f"`{key}` not found in configuration file. Falling back to default value."
    )
    assert out == f"{ytcfg_defaults['yt', key]}\n"
    assert tmp_config_file.read_text() == ""


def test_get_unknown_parameter(tmp_config_file, capsys):
    assert tmp_config_file.read_text() == ""
    retval = run_main(["config", "get", "unknown_parameter"])
    out, err = capsys.readouterr()

    assert type(retval) is int
    assert retval != 0
    err = err.splitlines()
    assert "INFO" in err[0]
    assert str(tmp_config_file) in err[0]
    assert err[1] == "Error: unknown parameter `unknown_parameter`."
    assert out == ""
    assert tmp_config_file.read_text() == ""


# TODO: yt config rm tests...
