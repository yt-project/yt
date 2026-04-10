import os
import sys
import tempfile
import types
import warnings
from io import StringIO
from pathlib import Path
from unittest import mock

import pytest

import yt.utilities.configure as cfg_module
from yt.utilities.configuration_tree import ConfigNode


def test_config_dir_and_cast_helpers_cover_bool_numeric_and_string_cases():
    with tempfile.TemporaryDirectory() as tmpdir:
        with mock.patch.dict(os.environ, {"XDG_CONFIG_HOME": tmpdir}, clear=False):
            assert cfg_module.config_dir() == os.path.join(tmpdir, "yt")

        with mock.patch.dict(os.environ, {"HOME": tmpdir}, clear=False):
            assert cfg_module._expand_all("~/yt-data") == os.path.join(tmpdir, "yt-data")

    assert cfg_module._cast_bool_helper("true") is True
    assert cfg_module._cast_bool_helper("False") is False
    with pytest.raises(ValueError, match="Cannot safely cast to bool"):
        cfg_module._cast_bool_helper("maybe")

    assert cfg_module._cast_value_helper("True") is True
    assert cfg_module._cast_value_helper("10") == 10
    assert cfg_module._cast_value_helper("3.5") == 3.5
    assert cfg_module._cast_value_helper("plain-text") == "plain-text"
    assert cfg_module._cast_value_helper("plain-text", types=(int,)) is None


def test_ytconfig_tree_access_specificity_and_item_protocol():
    cfg = cfg_module.YTConfig(defaults={"yt": {}})
    assert not cfg.has_section("yt")
    assert cfg.remove_section("yt") is False

    cfg.add_section("yt")
    assert cfg.has_section("yt")
    assert isinstance(cfg.get("yt"), ConfigNode)

    cfg.update({"analysis": {"width": 2.5}})
    assert cfg.get("analysis", "width") == 2.5

    cfg.set("yt", "answer", 42, metadata={"source": "manual"})
    cfg.set("yt", "internals", "answer", 10)
    cfg["yt", "internals", "parallel"] = True

    assert cfg.get("yt", "answer") == 42
    assert (
        cfg.get("yt", "answer", callback=lambda leaf: leaf.extra_data["source"])
        == "manual"
    )
    assert cfg["yt", "internals", "parallel"] is True
    assert "yt" in cfg
    assert cfg._repr_json_() == cfg.config_root._repr_json_()

    assert cfg.get_most_specific("yt", "internals", "missing", "answer") == 10
    assert (
        cfg.get_most_specific("yt", "internals", "missing", "parallel", fallback=False)
        is True
    )
    assert (
        cfg.get_most_specific("yt", "internals", "missing", "not_here", fallback="none")
        == "none"
    )
    with pytest.raises(KeyError):
        cfg.get_most_specific("yt", "internals", "missing", "not_here")

    assert cfg.remove_section("analysis") is True
    assert not cfg.has_section("analysis")


def test_ytconfig_read_covers_missing_valid_invalid_and_legacy_tomli_branch():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        missing = tmp_path / "missing.toml"
        valid = tmp_path / "valid.toml"
        invalid = tmp_path / "invalid.toml"

        valid.write_text(
            "[yt]\nanswer = 42\n\n[yt.internals]\nparallel = true\n",
            encoding="ascii",
        )
        invalid.write_text("[yt\nanswer = 1\n", encoding="ascii")

        cfg = cfg_module.YTConfig()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            read_files = cfg.read([str(missing), str(valid), str(invalid)])

        assert read_files == [str(valid)]
        assert cfg.get("yt", "answer") == 42
        assert (
            cfg.get("yt", "answer", callback=lambda leaf: leaf.extra_data["source"])
            == f"file: {valid}"
        )
        assert any("Could not load configuration file" in str(w.message) for w in caught)

        import tomllib

        fake_tomli = types.SimpleNamespace(
            load=tomllib.load, TOMLDecodeError=tomllib.TOMLDecodeError
        )
        with (
            mock.patch.object(cfg_module.sys, "version_info", (3, 10, 0)),
            mock.patch.dict(sys.modules, {"tomli": fake_tomli}),
        ):
            cfg_old = cfg_module.YTConfig()
            assert cfg_old.read(str(valid)) == [str(valid)]
            assert cfg_old.get("yt", "internals", "parallel") is True


def test_ytconfig_write_covers_file_objects_paths_and_invalid_targets():
    cfg = cfg_module.YTConfig()
    cfg.set("yt", "answer", 42)
    cfg.set("yt", "internals", "parallel", True)

    out = StringIO()
    cfg.write(out)
    assert out.getvalue() == "[yt]\nanswer = 42\n\n[yt.internals]\nparallel = true\n"

    with tempfile.TemporaryDirectory() as tmpdir:
        target = Path(tmpdir) / "nested" / "yt.toml"
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            cfg.write(target)

        assert target.read_text(encoding="utf-8") == out.getvalue()
        assert any("does not exist, creating it" in str(w.message) for w in caught)

    with pytest.raises(TypeError, match="Expected a path to a file, or a writable object"):
        cfg.write(object())


def test_config_file_helpers_and_module_level_wrappers():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        config_root = tmp_path / "config-home"
        project_root = tmp_path / "project"
        nested = project_root / "subdir" / "deeper"
        nested.mkdir(parents=True)

        local_file = project_root / "yt.toml"
        local_file.write_text("[yt]\n", encoding="ascii")

        original_cwd = Path.cwd()
        os.chdir(nested)
        try:
            with mock.patch.dict(
                os.environ, {"XDG_CONFIG_HOME": str(config_root)}, clear=False
            ):
                assert cfg_module.YTConfig.get_global_config_file() == os.path.join(
                    str(config_root), "yt", "yt.toml"
                )
                assert cfg_module.YTConfig.get_local_config_file() == str(
                    local_file.resolve()
                )

                local_file.unlink()
                assert cfg_module.YTConfig.get_local_config_file() == os.path.join(
                    os.path.abspath(os.curdir), "yt.toml"
                )
        finally:
            os.chdir(original_cwd)

        fresh_config = cfg_module.YTConfig()
        wrapped_file = tmp_path / "wrapped.toml"
        with mock.patch.object(cfg_module, "CONFIG", fresh_config):
            cfg_module.set_config("yt", "internals.parallel", "True", wrapped_file)
            assert cfg_module.get_config("yt", "internals.parallel") is True
            cfg_module.write_config(wrapped_file)
            assert "parallel = true" in wrapped_file.read_text(encoding="utf-8")

            cfg_module.rm_config("yt", "internals.parallel", wrapped_file)
            assert wrapped_file.read_text(encoding="utf-8") == ""
            with pytest.raises(KeyError):
                cfg_module.get_config("yt", "internals.parallel")
