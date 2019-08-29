def requires_ds(x):
    def ffalse(func):
        return lambda: None
    def ftrue(func):
        return func
    if x:
        return ftrue
    else:
        return ffasle


def my_func(fn):
    my_other_func(fn)
