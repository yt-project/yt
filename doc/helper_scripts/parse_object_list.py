import inspect
from textwrap import TextWrapper

import yt

ds = yt.load("RD0005-mine/RedshiftOutput0005")

output = open("source/analyzing/_obj_docstrings.inc", "w")

template = """

.. class:: %(clsname)s%(sig)s:

   For more information, see :ref:`%(docstring)s`
   (This is a proxy for :class:`~%(clsproxy)sBase`.)
"""

tw = TextWrapper(initial_indent="   ", subsequent_indent="   ", width=60)


def write_docstring(f, name, cls):
    for clsi in inspect.getmro(cls):
        docstring = inspect.getdoc(clsi.__init__)
        if docstring is not None:
            break
    clsname = name
    sig = inspect.formatargspec(*inspect.getargspec(cls.__init__))
    sig = sig.replace("**kwargs", "**field_parameters")
    clsproxy = f"yt.data_objects.data_containers.{cls.__name__}"
    f.write(
        template
        % dict(
            clsname=clsname, sig=sig, clsproxy=clsproxy, docstring="physical-object-api"
        )
    )


for n, c in sorted(ds.__dict__.items()):
    if hasattr(c, "_con_args"):
        print(n)
        write_docstring(output, n, c)
