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


def issue_deprecation_warning(msg, stacklevel=3):
    warnings.warn(msg, VisibleDeprecationWarning, stacklevel=stacklevel)


def issue_demeshening_deprecation_warning(msg):
    msg = " ".join(
        [
            msg,
            "Since yt-4.0, it's no longer necessary to add a field specifically for "
            "smoothing, because the global octree is removed. The old behavior of "
            "interpolating onto a grid structure can be recovered through data objects "
            "like ds.arbitrary_grid, ds.covering_grid, and most closely ds.octree. The "
            "visualization machinery now treats SPH fields properly by smoothing onto "
            "pixel locations. See this page to learn more: "
            "https://yt-project.org/doc/yt4differences.html",
        ]
    )
    issue_deprecation_warning(msg)


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
