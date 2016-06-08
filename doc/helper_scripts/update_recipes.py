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
    print("Parsing %s" % (file))
    lines = ctx[file].data().split("\n")
    fn = file[8:-3]
    title = fn.replace("_", " ").capitalize()
    title += "\n" + "-" * len(title) + "\n"*2
    title = ".. _cookbook-%s:\n\n%s" % (fn, title)
    if lines[0] != '"""':
        print("    Bad docstring: breaking.")
        print(file)
    di = lines[1:].index('"""')
    docstring = lines[1:di+1]
    recipe = lines[di+2:]
    recipes.write(".. include:: %s.inc\n" % fn)
    output = open("source/cookbook/%s.inc" % fn, "w")
    output.write(title)
    output.write("\n".join(docstring))
    output.write("\n")
    output.write("\nThe latest version of this recipe can be downloaded ")
    output.write(url % file)
    output.write("\n\n.. code-block:: python\n\n")
    for line in recipe:
        output.write("   " + line + "\n")
    if os.path.isdir("../cookbook/images/%s" % (fn)):
        output.write("\n")
        ndir = "source/cookbook/_%s/" % (fn)
        if not os.path.isdir(ndir): os.mkdir(ndir)
        written = False
        for ifn in sorted(glob.glob("../cookbook/images/%s/*.png" % (fn))):
            written = cond_output(output, written)
            ofn = "%s/%s_%s" % (ndir, fn, os.path.basename(ifn))
            open(ofn, "wb").write(open(ifn, "rb").read())
            output.write(".. image:: _%s/%s_%s\n" % (fn, fn, os.path.basename(ifn)) +
                         "   :width: 240\n" +
                         "   :target: ../_images/%s_%s\n" % (fn, os.path.basename(ifn))
                        )
        for ifn in sorted(glob.glob("../cookbook/images/%s/*.txt" % (fn))):
            written = cond_output(output, written)
            open("%s/%s_%s" % (ndir, fn, os.path.basename(ifn)), "w").write(open(ifn, "r").read())
            output.write("\n``%s``\n\n" % os.path.basename(ifn))
            output.write(".. literalinclude:: _%s/%s_%s\n" % (fn, fn, os.path.basename(ifn)))
    output.write("\n\n")
    output.close()

recipes.write(footer)
recipes.close()
