# future division is important to divide integers and get as
# a result precise floating numbers (instead of truncated int)
from __future__ import absolute_import

__author__ = "github.com/casperdcl"
__all__ = ["tqdm_pandas"]


def tqdm_pandas(t):  # pragma: no cover
    """
    Registers the given `tqdm` instance with
    `pandas.core.groupby.DataFrameGroupBy.progress_apply`.
    It will even close() the `tqdm` instance upon completion.

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> from tqdm import tqdm, tqdm_pandas
    >>>
    >>> df = pd.DataFrame(np.random.randint(0, 100, (100000, 6)))
    >>> tqdm_pandas(tqdm())  # can use tqdm_gui, optional kwargs, etc
    >>> # Now you can use `progress_apply` instead of `apply`
    >>> df.groupby(0).progress_apply(lambda x: x**2)

    References
    ----------
    https://stackoverflow.com/questions/18603270/
    progress-indicator-during-pandas-operations-python
    """
    from pandas.core.groupby import DataFrameGroupBy

    def inner(groups, func, *args, **kwargs):
        """
        Parameters
        ----------
        groups  : DataFrameGroupBy
            Grouped data.
        func  : function
            To be applied on the grouped data.

        *args and *kwargs are transmitted to DataFrameGroupBy.apply()
        """
        t.total = len(groups) + 1  # pandas calls update once too many

        def wrapper(*args, **kwargs):
            t.update()
            return func(*args, **kwargs)

        result = groups.apply(wrapper, *args, **kwargs)

        t.close()

        return result

    # Enable custom tqdm progress in pandas!
    DataFrameGroupBy.progress_apply = inner
