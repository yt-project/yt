import inspect
from textwrap import TextWrapper

import yt

ds = yt.load("RD0005-mine/RedshiftOutput0005")

output = open("source/analyzing/_dq_docstrings.inc", "w")

template = """

.. function:: %(funcname)s%(sig)s:

   (This is a proxy for :func:`~%(funcproxy)s`.)
%(docstring)s

"""

tw = TextWrapper(initial_indent="   ", subsequent_indent="   ", width=60)


def write_docstring(f, name, func):
    docstring = inspect.getdoc(func)
    funcname = name
    sig = inspect.formatargspec(*inspect.getargspec(func))
    sig = sig.replace("data, ", "")
    sig = sig.replace("(data)", "()")
    funcproxy = f"yt.data_objects.derived_quantities.{func.__name__}"
    docstring = "\n".join("   %s" % line for line in docstring.split("\n"))
    f.write(
        template
        % dict(funcname=funcname, sig=sig, funcproxy=funcproxy, docstring=docstring)
    )
    # docstring = "\n".join(tw.wrap(docstring))))


dd = ds.all_data()
for n, func in sorted(dd.quantities.functions.items()):
    print(n, func)
    write_docstring(output, n, func[1])
