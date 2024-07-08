from importlib.metadata import version
from importlib.util import find_spec

from packaging.version import Version

__all__ = [
    "IS_IPYTHON",
    "IPYWIDGETS_ENABLED",
]

IS_IPYTHON: bool
HAS_IPYWIDGETS_GE_8: bool
IPYWIDGETS_ENABLED: bool


try:
    # this name is only defined if running within ipython/jupyter
    __IPYTHON__  # type: ignore [name-defined] # noqa: B018
except NameError:
    IS_IPYTHON = False
else:
    IS_IPYTHON = True


HAS_IPYWIDGETS_GE_8 = (
    Version(version("ipywidgets")) >= Version("8.0.0")
    if find_spec("ipywidgets") is not None
    else False
)

IPYWIDGETS_ENABLED = IS_IPYTHON and HAS_IPYWIDGETS_GE_8
