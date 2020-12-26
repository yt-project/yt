import warnings
from functools import wraps


def deprecated_class(cls):
    @wraps(cls)
    def _func(*args, **kwargs):
        # Note we use SyntaxWarning because by default, DeprecationWarning is
        # not shown.
        warnings.warn(
            f"This usage is deprecated.  Please use {cls.__name__} instead.",
            SyntaxWarning,
            stacklevel=2,
        )
        return cls(*args, **kwargs)

    return _func


class VisibleDeprecationWarning(UserWarning):
    """Visible deprecation warning, adapted from NumPy

    By default python does not show users deprecation warnings.
    This ensures that a deprecation warning is visible to users
    if that is desired.
    """

    pass


def issue_deprecation_warning(msg, *, removal, since=None, stacklevel=3):
    """
    `since` and `removal`: str version numbers, indicating the anticipated removal date

    `since` is optional only because it was introduced in the 4.0.0 release, and it
    should become mandatory in the future.
    Both `since` and `removal` are keyword-only arguments so that their order can be
    swapped in the future without introducing bugs.
    """
    msg += "\n"
    if since is not None:
        msg += f"Deprecated since v{since}."

    msg += f"this feature will be removed in v{removal}"
    warnings.warn(msg, VisibleDeprecationWarning, stacklevel=stacklevel)


def deprecate(replacement):
    def real_deprecate(func):
        """
        This decorator issues a deprecation warning.

        This can be used like so:

        .. code-block:: python

        @deprecate("new_function")
        def some_really_old_function(...):

        """

        @wraps(func)
        def run_func(*args, **kwargs):
            message = "%s has been deprecated and may be removed without notice!"
            if replacement is not None:
                message += f" Use {replacement} instead."

            issue_deprecation_warning(message % func.__name__, stacklevel=2)
            func(*args, **kwargs)

        return run_func

    return real_deprecate
