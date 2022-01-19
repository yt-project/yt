import collections

templates = dict(
    autoclass=r"""
%(name)s
%(header)s

.. autoclass:: %(name)s
   :members:
   :inherited-members:
   :undoc-members:

""",
    autofunction=r"""
%(name)s
%(header)s

.. autofunction:: %(name)s

""",
    index_file=r"""
%(title)s
%(header)s

.. autosummary::
   :toctree: generated/%(dn)s

""",
)

file_names = dict(
    ft=("Field Types", "source/api/field_types/%s.rst"),
    pt=("Plot Types", "source/api/plot_types/%s.rst"),
    cl=("Callback List", "source/api/callback_list/%s.rst"),
    ee=("Extension Types", "source/api/extension_types/%s.rst"),
    dd=("Derived Datatypes", "source/api/derived_datatypes/%s.rst"),
    mt=("Miscellaneous Types", "source/api/misc_types/%s.rst"),
    fl=("Function List", "source/api/function_list/%s.rst"),
    ds=("Data Sources", "source/api/data_sources/%s.rst"),
    dq=("Derived Quantities", "source/api/derived_quantities/%s.rst"),
)

to_include = collections.defaultdict(list)

for line in open("auto_generated.txt"):
    ftype, name, file_name = (s.strip() for s in line.split("::"))
    cn = name.split(".")[-1]
    if cn[0] == "_":
        cn = cn[1:]  # For leading _
    fn = file_names[file_name][1] % cn
    # if not os.path.exists(os.path.dirname(fn)):
    #    os.mkdir(os.path.dirname(fn))
    header = "-" * len(name)
    dd = dict(header=header, name=name)
    # open(fn, "w").write(templates[ftype] % dd)
    to_include[file_name].append(name)

for key, val in file_names.items():
    title, file = val
    fn = file.rsplit("/", 1)[0] + ".rst"
    print(fn)
    f = open(fn, "w")
    dn = fn.split("/")[-1][:-4]
    dd = dict(header="=" * len(title), title=title, dn=dn)
    f.write(templates["index_file"] % dd)
    for obj in sorted(to_include[key]):
        f.write(f"   {obj}\n")
