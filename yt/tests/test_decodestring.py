def test__rdbeta_py3_does_not_crash():
    from yt.funcs import _rdbeta

    _rdbeta("k")
