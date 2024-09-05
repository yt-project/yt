from yt.loaders import load

DEFAULT_ATTRS = (
    "dimensionality",
    "refine_by",
    "domain_dimensions",
    "current_time",
    "domain_left_edge",
    "domain_right_edge",
    "unique_identifier",
    "current_redshift",
    "cosmological_simulation",
    "omega_matter",
    "omega_lambda",
    "hubble_constant",
    "dataset_type",
)


def get_metadata(path, full_output=False, attrs=DEFAULT_ATTRS):
    ds = load(path)
    metadata = {"filename": path}
    for a in attrs:
        v = getattr(ds, a, None)
        if v is None:
            continue
        if hasattr(v, "tolist"):
            v = v.tolist()
        metadata[a] = v
    if full_output:
        params = {}
        for p, v in ds.parameters.items():
            if hasattr(v, "tolist"):
                v = v.tolist()
            params[p] = v
        metadata["params"] = params
    ds.close()
    return metadata
