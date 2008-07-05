"""
All file-handling takes place in here.

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2007-2008 Matthew Turk.  All Rights Reserved.

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


from yt.fido import *

def copy_output(basename, new_location, extra_files=None):
    if extra_files == None: extra_files = []
    copy_glob("%s*" % (basename), new_location)
    for file in extra_files:
        copy_glob(file, location)
    return os.path.join(new_location, os.path.basename(basename))

def move_output(basename, new_location, extra_files=None):
    if extra_files == None: extra_files = []
    move_glob("%s*" % (basename), new_location)
    for file in extra_files:
        move_glob(file, location)

def delete_output(basename, extra_files=None):
    if extra_files == None: extra_files = []
    delete_glob("%s*" % os.path.normpath(filename), location)
    for file in extra_files:
        delete_glob(file)

def bury_output(filename, new_prefix = None, new_location=None, extra_files=None):
    if extra_files == None: extra_files = []
    if not new_location:
        dirname = NewDirectoryPattern % filename
        if not new_prefix: new_prefix = os.getcwd()
        new_location = os.path.abspath(os.path.join(new_prefix, dirname))
    move_glob("%s*" % filename, new_location)
    for file in extra_files:
        move_glob(file, new_location)
    return os.path.join(new_location, filename)

def digup_output(filename, new_location=None, extra_files=None):
    if extra_files == None: extra_files = []
    if not new_location:
        new_location = get_parent_dir(filename)
    else:
        # Make sure we have no relative references
        new_location = os.path.normpath(new_location)
    print "%s*" % (filename), new_location
    move_glob("%s*" % filename, new_location)

def copy_glob(glob_pattern, new_location):
    if not os.path.isdir(new_location):
        os.makedirs(new_location)
    for file in glob.glob(glob_pattern):
        mylog.debug("Copying %s to %s", file, new_location)
        shutil.copy(file, new_location)
    
def move_glob(glob_pattern, new_location):
    if not os.path.isdir(new_location):
        os.makedirs(new_location)
    # This is slower than dumping to shell.
    # But, I trust it.  And it should work on boh
    # BSD and GNU util systems.
    for file in glob.glob(glob_pattern):
        if os.path.abspath(file) == os.path.abspath(new_location): continue
        nl = os.path.join(new_location, os.path.basename(file))
        mylog.debug("Moving %s to %s", file, nl)
        shutil.move(file, nl)

def delete_glob(glob_pattern):
    for file in glob.glob(glob_pattern):
        print "Removing %s" % (file)
        shutil.unlink(file)
