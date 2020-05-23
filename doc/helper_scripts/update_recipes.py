header = """Cookbook
========

yt scripts can be a bit intimidating, and at times a bit obtuse.  But there's a
lot you can do, and this section of the manual will assist with figuring out
how to do some fairly common tasks -- which can lead to combining these, with
other Python code, into more complicated and advanced tasks.

.. note::
   All of these scripts are located in the mercurial repository at
   http://bitbucket.org/yt_analysis/cookbook/

"""
footer = """ """

import cStringIO, sys, glob, os
from mercurial import hg, ui, commands
uii = ui.ui()

if "--index" in sys.argv:
    recipes = open("source/cookbook/recipes.rst", "w")
else:
    recipes = cStringIO.StringIO()
recipes.write(header)

url = "here: http://bitbucket.org/yt_analysis/cookbook/raw/tip/%s ."

def cond_output(f, v):
    if not v:
        f.write(".. rubric:: Sample Output\n\n")
    return True

repo = hg.repository(uii, "../cookbook/")
commands.pull(uii, repo, "http://bitbucket.org/yt_analysis/cookbook/")
ctx = repo["tip"]
for file in ctx:
    if not file.startswith("recipes/"): continue
    print(f"Parsing {file}")
    lines = ctx[file].data().split("\n")
    fn = file[8:-3]
    title = fn.replace("_", " ").capitalize()
    title += "\n" + "-" * len(title) + "\n"*2
    title = f".. _cookbook-{fn}:\n\n{title}"
    if lines[0] != '"""':
        print("    Bad docstring: breaking.")
        print(file)
    di = lines[1:].index('"""')
    docstring = lines[1:di+1]
    recipe = lines[di+2:]
    recipes.write(f".. include:: {fn}.inc\n")
    output = open(f"source/cookbook/{fn}.inc", "w")
    output.write(title)
    output.write("\n".join(docstring))
    output.write("\n")
    output.write("\nThe latest version of this recipe can be downloaded ")
    output.write(url % file)
    output.write("\n\n.. code-block:: python\n\n")
    for line in recipe:
        output.write("   " + line + "\n")
    if os.path.isdir(f"../cookbook/images/{fn}"):
        output.write("\n")
        ndir = f"source/cookbook/_{fn}/"
        if not os.path.isdir(ndir): os.mkdir(ndir)
        written = False
        for ifn in sorted(glob.glob(f"../cookbook/images/{fn}/*.png")):
            written = cond_output(output, written)
            ofn = f"{ndir}/{fn}_{os.path.basename(ifn)}"
            open(ofn, "wb").write(open(ifn, "rb").read())
            output.write(f".. image:: _{fn}/{fn}_{os.path.basename(ifn)}\n" +
                         "   :width: 240\n" +
                         f"   :target: ../_images/{fn}_{os.path.basename(ifn)}\n"
                        )
        for ifn in sorted(glob.glob(f"../cookbook/images/{fn}/*.txt")):
            written = cond_output(output, written)
            open(f"{ndir}/{fn}_{os.path.basename(ifn)}", "w").write(open(ifn, "r").read())
            output.write("\n``%s``\n\n" % os.path.basename(ifn))
            output.write(f".. literalinclude:: _{fn}/{fn}_{os.path.basename(ifn)}\n")
    output.write("\n\n")
    output.close()

recipes.write(footer)
recipes.close()
