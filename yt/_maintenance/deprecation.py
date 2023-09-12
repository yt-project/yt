import warnings
from functools import wraps
from types import FunctionType
from typing import Optional


def issue_deprecation_warning(
    msg: str,
    *,
    stacklevel: int,
    since: str,
    removal: Optional[str] = None,
):
    """
    Parameters
    ----------
    msg : str
        A text message explaining that the code surrounding the call to this function is
        deprecated, and what should be changed on the user side to avoid it.

    stacklevel: int
        Number of stack frames to be skipped when pointing at caller code, starting from
        *this* function's frame. In general 3 is a minimum.

    since and removal: str version numbers, indicating the anticipated removal date

    Notes
    -----

    removal can be left empty if it is not clear how many minor releases are expected to
    happen before the next major.

    removal and since arguments are keyword-only to forbid accidentally swapping them.

    Examples
    --------
    >>> issue_deprecation_warning(
    ...     "This code is deprecated.", stacklevel=3, since="4.0"
    ... )
    """

    msg += f"\nDeprecated since yt {since}"
    if removal is not None:
        msg += f"\nThis feature is planned for removal in yt {removal}"
    warnings.warn(msg, DeprecationWarning, stacklevel=stacklevel)


def future_positional_only(positions2names: dict[int, str], /, **depr_kwargs):
    """Warn users when using a future positional-only argument as keyword.
    Note that positional-only arguments are available from Python 3.8
    See https://www.python.org/dev/peps/pep-0570/
    """

    def outer(func: FunctionType):
        @wraps(func)
        def inner(*args, **kwargs):
            for no, name in sorted(positions2names.items()):
                if name not in kwargs:
                    continue
                value = kwargs[name]
                issue_deprecation_warning(
                    f"Using the {name!r} argument as keyword (on position {no}) "
                    "is deprecated. "
                    "Pass the argument as positional to suppress this warning, "
                    f"i.e., use {func.__name__}({value!r}, ...)",
                    stacklevel=3,
                    **depr_kwargs,
                )
            return func(*args, **kwargs)

        return inner

    return outer
