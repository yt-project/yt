"""
A simple, Enzo-aware interface to drop.io based on another drop.io API.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2009 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from yt.config import ytcfg
from yt.mods import *
import zipfile, tempfile, glob, os, stat
from dropio.client import DropIoClient

API_KEY = ytcfg.get("yt", "dropiokey")

class SharedDataProduct(object):
    _upload_filename = None
    _temp_file = None
    def upload(self, drop_file_name, drop_name = None):
        self._create_archive()
        # This is where our drop.io API calls go
        client = DropIoClient(API_KEY)
        if drop_name is None:
            drop = client.create_drop
        else:
            drop = client.get_drop(drop_name)
            if drop is None: drop = client.create_drop(drop_name)
        drop_name = drop.name
        for asset in client.get_asset_list(drop_name):
            if asset.name != drop_file_name: continue
            print "FOUND THE DROP FILE NAME IN THE ASSET LIST"
            client.delete_asset(drop_name, drop_file_name)
        fsize = os.stat(self._temp_file.name)[stat.ST_SIZE]
        print "Uploading.  Total size is %0.2fMB." % (fsize/(1024*1024))
        print "This may take a considerable amount of time."
        asset = client.create_file(
                    drop_name, self._temp_file.name, 
                    remote_file_name = drop_file_name)
        return asset

    def _create_archive(self):
        self._temp_file = tempfile.NamedTemporaryFile(mode="wb")
        zf = zipfile.ZipFile(self._temp_file, mode='w')
        for arcname, fn in self._filenames():
            # Note that this defaults to zipfile.ZIP_STORED as a compression
            # method.  This might be the safest method, as some platforms have
            # trouble with zlib.
            print "Adding", arcname, fn
            zf.write(fn, arcname)
        zf.close()
        self._temp_file.flush()

class SharedStaticOutput(SharedDataProduct):
    def __init__(self, pf):
        self.pf = pf

    def _filenames(self):
        gpatt = glob.glob("%s/%s*" % ( self.pf.fullpath, self.pf.basename))
        arcprefix = os.path.basename(self.pf.fullpath)
        print arcprefix
        for fn in sorted(gpatt):
            aname = os.path.join(arcprefix, os.path.basename(fn))
            yield aname, fn

class SharedDataStore(SharedDataProduct):
    pass
