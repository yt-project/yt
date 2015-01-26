#!/usr/bin/python

import os
import glob
import sys
import shutil
import tempfile
import traceback
import subprocess
import matplotlib
matplotlib.use('Agg')
from multiprocessing import Pool
from yt.config import ytcfg

FPATTERNS = ['*.png', '*.txt', '*.h5', '*.dat']
DPATTERNS = ['LC*', 'LR', 'DD0046', 'halo_analysis']
BADF = ['cloudy_emissivity.h5', 'apec_emissivity.h5',
        'xray_emissivity.h5', 'AMRGridData_Slice_x_density.png']
CWD = os.getcwd()
ytcfg["yt", "serialize"] = "False"
PARALLEL_TEST = {"rockstar_nest": "3"}
BLACKLIST = []


def prep_dirs():
    for directory in glob.glob('%s/*' % ytcfg.get("yt", "test_data_dir")):
        os.symlink(directory, os.path.basename(directory))


def run_receipe((receipe,)):
    module_name, ext = os.path.splitext(os.path.basename(receipe))
    dest = os.path.join(os.path.dirname(receipe), '_static', module_name)
    if module_name in BLACKLIST:
        return 0
    if not os.path.exists("%s/_temp/%s.done" % (CWD, module_name)):
        sys.stderr.write('Started %s\n' % module_name)
        tmpdir = tempfile.mkdtemp()
        os.chdir(tmpdir)
        prep_dirs()
        if module_name in PARALLEL_TEST:
            cmd = ["mpiexec", "-n", PARALLEL_TEST[module_name],
                   "python", receipe]
        else:
            cmd = ["python", receipe]
        try:
            subprocess.check_call(cmd)
        except:
            trace = "".join(traceback.format_exception(*sys.exc_info()))
            trace += " in module: %s\n" % module_name
            trace += " receipe: %s\n" % receipe
            raise Exception(trace)
        open("%s/_temp/%s.done" % (CWD, module_name), 'wb').close()
        for pattern in FPATTERNS:
            for fname in glob.glob(pattern):
                if fname not in BADF:
                    shutil.move(fname, "%s__%s" % (dest, fname))
        for pattern in DPATTERNS:
            for dname in glob.glob(pattern):
                shutil.move(dname, dest)
        os.chdir(CWD)
        shutil.rmtree(tmpdir, True)
        sys.stderr.write('Finished with %s\n' % module_name)
    return 0

for path in ['_temp', 'source/cookbook/_static',
             'source/visualizing/colormaps/_static']:
    fpath = os.path.join(CWD, path)
    if os.path.exists(fpath):
        shutil.rmtree(fpath)
    os.makedirs(fpath)

os.chdir('_temp')
receipes = []
for rpath in ['source/cookbook', 'source/visualizing/colormaps']:
    fpath = os.path.join(CWD, rpath)
    sys.path.append(fpath)
    receipes += glob.glob('%s/*.py' % fpath)
WPOOL = Pool(processes=6)
RES = WPOOL.map_async(run_receipe, ((receipe,) for receipe in receipes))
RES.get()
os.chdir(CWD)
