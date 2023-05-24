from yt._maintenance.deprecation import issue_deprecation_warning

from .unions import ParticleUnion  # noqa: F401

issue_deprecation_warning(
    "Importing ParticleUnion from yt.data_objects.particle_unions is deprecated. "
    "Please import this class from yt.data_objects.unions instead",
    since="4.2",
)
