import warnings
from typing import Optional


class VisibleDeprecationWarning(UserWarning):
    """Visible deprecation warning, adapted from NumPy

    The nose runner does not show users DeprecationWarning.
    This ensures that a deprecation warning is visible to users
    if that is desired.
    """

    # this class becomes useless after the tests are migrated from nose to pytest

    pass


def issue_deprecation_warning(
    msg: str, *, since: str, removal: Optional[str] = None, stacklevel: int = 3
):
    """
    Parameters
    ----------
    msg : str
        A text message explaining that the code surrounding the call to this function is
        deprecated, and what should be changed on the user side to avoid it.

    since and removal: str version numbers, indicating the anticipated removal date

    Notes
    -----

    removal can be left empty if it is not clear how many minor releases are expected to
    happen before the next major.

    removal and since arguments are keyword-only to forbid accidentally swapping them.

    Examples
    --------
    >>> issue_deprecation_warning(
    ...     "This code is deprecated.", since="4.0.0", removal="4.2.0"
    ... )
    """

    msg += f"\nDeprecated since yt {since}\nThis feature is planned for removal "
    if removal is None:
        msg += "two minor releases later (anticipated)"
    else:
        msg += f"in yt {removal}"
    warnings.warn(msg, VisibleDeprecationWarning, stacklevel=stacklevel)
