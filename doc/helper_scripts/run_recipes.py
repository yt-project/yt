#!/usr/bin/python

import os
import glob
import sys
import shutil
import tempfile
from multiprocessing import Pool
from yt.config import ytcfg

FPATTERNS = ['*.png', '*.txt', '*.h5', '*.dat']
DPATTERNS = ['LC*', 'LR', 'DD0046', 'halo_analysis']
CWD = os.getcwd()
DIR = 'source/cookbook/_static'
DESTDIR = os.path.join(CWD, DIR)
#NEEDS_SERIAL = ["light_cone_with_halo_mask"]
ytcfg["yt", "serialize"] = "False"


def prep_dirs():
    for directory in glob.glob('%s/*' % ytcfg.get("yt", "test_data_dir")):
        os.symlink(directory, os.path.basename(directory))


def run_receipe((receipe,)):
    module_name, ext = os.path.splitext(os.path.basename(receipe))
    if not os.path.exists("%s/_temp/%s.done" % (CWD, module_name)):
        sys.stderr.write('Started %s\n' % module_name)
        cwd = os.getcwd()
        tmpdir = tempfile.mkdtemp()
        os.chdir(tmpdir)
        prep_dirs()
        module = __import__(module_name)
        open("%s/_temp/%s.done" % (CWD, module_name), 'wb').close()
        for pattern in FPATTERNS:
            for fname in glob.glob(pattern):
                dst = os.path.join(DESTDIR, module_name)
                shutil.move(fname, "%s__%s" % (dst, fname))
        for pattern in DPATTERNS:
            for dname in glob.glob(pattern):
                dst = os.path.join(DESTDIR, module_name)
                shutil.move(dname, dst)
        os.chdir(cwd)
        shutil.rmtree(tmpdir, True)
        sys.stderr.write('Finished with %s\n' % module_name)
    return 0

for path in ['_temp', DESTDIR]:
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)

os.chdir('_temp')
sys.path.append(os.path.join(CWD, 'source/cookbook'))
WPOOL = Pool(processes=8)
RES = WPOOL.map_async(run_receipe, (
    (receipe,) for receipe in glob.glob('%s/*.py' % sys.path[-1])))
RES.get()
os.chdir(CWD)
