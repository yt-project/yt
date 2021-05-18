import warnings

# This is a singleton object meant to be used as a default value for a deprecated
# argument in a function signature.
# Its name signals the deprecataion much better than a None and it is helpful to
# have a clear distinction with None, which is a perfectly fine default value, in particular
# for mutables.
#
# Intended usage:
#
# def foo(a, b_new=None, b_old=DEPRECATED_DEFAULT):
#     if b_old is not DEPRECATED_DEFAULT:
#         issue_deprecation_warning(...)
#         # possibly do something else as redirecting the value
#         # to another arguement, e.g.,
#         b_new = b_old
DEPRECATED_DEFAULT = object()


class VisibleDeprecationWarning(UserWarning):
    """Visible deprecation warning, adapted from NumPy

    The nose runner does not show users DeprecationWarning.
    This ensures that a deprecation warning is visible to users
    if that is desired.
    """

    # this class becomes useless after the tests are migrated from nose to pytest

    pass


def issue_deprecation_warning(msg, *, removal, since=None, stacklevel=3):
    """
    Parameters
    ----------
    msg : str
        A text message explaining that the code surrounding the call to this function is
        deprecated, and what should be changed on the user side to avoid it.

    since and removal: str version numbers, indicating the anticipated removal date

    Crucial note:
    beware that `removal` is required (it doesn't have a default value). This is
    vital since deprecated code is typically untested and not specifying a required
    keyword argument will turn the warning into a TypeError.
    What it gets us however is that the release manager will know for a fact whether it
    is safe to remove a feature at any given point, and users have a better idea when
    their code will become incompatible.

    `since` is optional only because it was introduced in the 4.0.0 release, and it
    should become mandatory in the future.
    Both `since` and `removal` are keyword-only arguments so that their order can be
    swapped in the future without introducing bugs.

    Examples
    --------
    >>> issue_deprecation_warning(
    ...     "This code is deprecated.", since="4.0.0", removal="4.2.0"
    ... )
    """
    msg += "\n"
    if since is not None:
        msg += f"Deprecated since v{since}. "

    msg += f"This feature will be removed in v{removal}"
    warnings.warn(msg, VisibleDeprecationWarning, stacklevel=stacklevel)
