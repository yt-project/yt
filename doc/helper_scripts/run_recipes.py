#!/usr/bin/env python3

import glob
import os
import shutil
import subprocess
import sys
import tempfile
import traceback
from multiprocessing import Pool

import matplotlib

from yt.config import ytcfg

matplotlib.use("Agg")

FPATTERNS = ["*.png", "*.txt", "*.h5", "*.dat", "*.mp4"]
DPATTERNS = ["LC*", "LR", "DD0046"]
BADF = [
    "cloudy_emissivity.h5",
    "apec_emissivity.h5",
    "xray_emissivity.h5",
    "AMRGridData_Slice_x_density.png",
]
CWD = os.getcwd()
ytcfg["yt", "serialize"] = False
BLACKLIST = ["opengl_ipython", "opengl_vr"]


def prep_dirs():
    for directory in glob.glob(f"{ytcfg.get('yt', 'test_data_dir')}/*"):
        os.symlink(directory, os.path.basename(directory))


def run_recipe(payload):
    (recipe,) = payload
    module_name, ext = os.path.splitext(os.path.basename(recipe))
    dest = os.path.join(os.path.dirname(recipe), "_static", module_name)
    if module_name in BLACKLIST:
        return 0
    if not os.path.exists(f"{CWD}/_temp/{module_name}.done"):
        sys.stderr.write(f"Started {module_name}\n")
        tmpdir = tempfile.mkdtemp()
        os.chdir(tmpdir)
        prep_dirs()
        try:
            subprocess.check_call(["python", recipe])
        except Exception:
            trace = "".join(traceback.format_exception(*sys.exc_info()))
            trace += f" in module: {module_name}\n"
            trace += f" recipe: {recipe}\n"
            raise Exception(trace)
        open(f"{CWD}/_temp/{module_name}.done", "wb").close()
        for pattern in FPATTERNS:
            for fname in glob.glob(pattern):
                if fname not in BADF:
                    shutil.move(fname, f"{dest}__{fname}")
        for pattern in DPATTERNS:
            for dname in glob.glob(pattern):
                shutil.move(dname, dest)
        os.chdir(CWD)
        shutil.rmtree(tmpdir, True)
        sys.stderr.write(f"Finished with {module_name}\n")
    return 0


for path in [
    "_temp",
    "source/cookbook/_static",
    "source/visualizing/colormaps/_static",
]:
    fpath = os.path.join(CWD, path)
    if os.path.exists(fpath):
        shutil.rmtree(fpath)
    os.makedirs(fpath)

os.chdir("_temp")
recipes = []
for rpath in ["source/cookbook", "source/visualizing/colormaps"]:
    fpath = os.path.join(CWD, rpath)
    sys.path.append(fpath)
    recipes += glob.glob(f"{fpath}/*.py")
WPOOL = Pool(processes=6)
RES = WPOOL.map_async(run_recipe, ((recipe,) for recipe in recipes))
RES.get()
os.chdir(CWD)
