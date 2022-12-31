import os

from yt.utilities.configure import YTConfig, config_dir, configuration_callbacks

ytcfg_defaults = {}

ytcfg_defaults["yt"] = dict(
    serialize=False,
    only_deserialize=False,
    time_functions=False,
    colored_logs=False,
    suppress_stream_logging=False,
    stdout_stream_logging=False,
    log_level=20,
    inline=False,
    num_threads=-1,
    store_parameter_files=False,
    parameter_file_store="parameter_files.csv",
    maximum_stored_datasets=500,
    skip_dataset_cache=True,
    load_field_plugins=False,
    plugin_filename="my_plugins.py",
    parallel_traceback=False,
    pasteboard_repo="",
    reconstruct_index=True,
    test_storage_dir="/does/not/exist",
    test_data_dir="/does/not/exist",
    enzo_db="",
    notebook_password="",
    answer_testing_tolerance=3,
    answer_testing_bitwise=False,
    gold_standard_filename="gold311",
    local_standard_filename="local001",
    answer_tests_url="http://answers.yt-project.org/{1}_{2}",
    sketchfab_api_key="None",
    imagebin_api_key="e1977d9195fe39e",
    imagebin_upload_url="https://api.imgur.com/3/image",
    imagebin_delete_url="https://api.imgur.com/3/image/{delete_hash}",
    curldrop_upload_url="http://use.yt/upload",
    thread_field_detection=False,
    ignore_invalid_unit_operation_errors=False,
    chunk_size=1000,
    xray_data_dir="/does/not/exist",
    supp_data_dir="/does/not/exist",
    default_colormap="cmyt.arbre",
    ray_tracing_engine="yt",
    internals=dict(
        within_testing=False,
        within_pytest=False,
        parallel=False,
        strict_requires=False,
        global_parallel_rank=0,
        global_parallel_size=1,
        topcomm_parallel_rank=0,
        topcomm_parallel_size=1,
        command_line=False,
    ),
)


# For backward compatibility, do not use these vars internally in yt
CONFIG_DIR = config_dir()


_global_config_file = YTConfig.get_global_config_file()
_local_config_file = YTConfig.get_local_config_file()

# Load the config
ytcfg = YTConfig()
ytcfg.update(ytcfg_defaults, metadata={"source": "defaults"})

# Try loading the local config first, otherwise fall back to global config
if os.path.exists(_local_config_file):
    ytcfg.read(_local_config_file)
elif os.path.exists(_global_config_file):
    ytcfg.read(_global_config_file)


def _setup_postinit_configuration():
    """This is meant to be run last in yt.__init__"""
    for callback in configuration_callbacks:
        callback(ytcfg)
