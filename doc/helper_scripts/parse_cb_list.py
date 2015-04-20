from yt.mods import *
import inspect
from textwrap import TextWrapper

ds = load("RD0005-mine/RedshiftOutput0005")

output = open("source/visualizing/_cb_docstrings.inc", "w")

template = """

.. function:: %(clsname)s%(sig)s:

   (This is a proxy for :class:`~%(clsproxy)s`.)

%(docstring)s

"""

tw = TextWrapper(initial_indent='   ', subsequent_indent='   ', width=60)

def write_docstring(f, name, cls):
    if not hasattr(cls, '_type_name') or cls._type_name is None: return
    for clsi in inspect.getmro(cls):
        docstring = inspect.getdoc(clsi.__init__)
        if docstring is not None: break
    clsname = cls._type_name
    sig = inspect.formatargspec(*inspect.getargspec(cls.__init__))
    sig = sig.replace("**kwargs", "**field_parameters")
    clsproxy = "yt.visualization.plot_modifications.%s" % (cls.__name__)
    #docstring = "\n".join(["   %s" % line for line in docstring.split("\n")])
    #print docstring
    f.write(template % dict(clsname = clsname, sig = sig, clsproxy=clsproxy,
                            docstring = "\n".join(tw.wrap(docstring))))
                            #docstring = docstring))

for n,c in sorted(callback_registry.items()):
    write_docstring(output, n, c)
    print(".. autoclass:: yt.visualization.plot_modifications.%s" % n)
    print("   :members:")
    print()
