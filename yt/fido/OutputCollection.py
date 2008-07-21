"""
Handling of sets of outputs.

We don't do any instantiation, or touching of lagos, etc.

Very simple nowadays.

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

from stat import ST_CTIME

class OutputCollection:
    def __init__(self, title):
        self.title = title
        self.output_names = na.array(())
        self.output_time_ids = na.array((), dtype='int64')
        self.output_times = na.array((), dtype='float64')

    def read_in(self, filename):
        self._last_read_filename = filename
        lines = open(filename).readlines()
        output_lines = [line for line in lines if line.startswith("Output:")]
        outputs = [line.split(":")[1] for line in output_lines]
        output_time_ids = [int(line.split(":")[2]) for line in output_lines]
        output_times = [float(line.split(":")[3]) for line in output_lines]

        self.output_names = na.array(outputs)
        self.output_time_ids = na.array(output_time_ids, dtype='int64')
        self.output_times = na.array(output_times, dtype='float64')

    def write_out(self):
        path=ytcfg.get("fido","rundir")
        if not os.path.isdir(path): os.makedirs(path)
        fn = os.path.join(path, "runF_%s" % (self.title))
        f = open(fn, "w")
        for i, output in enumerate(self.output_names):
            f.write("Output:%s:%s:%s\n" % ( \
                        os.path.abspath(output), \
                        self.output_time_ids[i], \
                        self.output_times[i]))
        f.close()

    def sort(self):
        order = na.argsort(self.output_times)
        self.output_names = self.output_names[order]
        self.output_times = self.output_times[order]
        self.output_time_ids = self.output_time_ids[order]

    def add_output(self, filename):
        # We're passing in *just* filenames here.  So, we simply snag the
        # two appropriate lines in the parameter file.
        if filename in self.output_names: return
        time = get_parameter_line(filename, "InitialTime").split()[-1]
        # Implement exception catching
        try:
            timeID = get_parameter_line(filename, "CurrentTimeIdentifier").split()[-1]
        except:
            timeID = int(os.stat(filename)[ST_CTIME])
        self.output_names = \
            na.array(self.output_names.tolist() + [filename])
        self.output_time_ids = \
            na.array(self.output_time_ids.tolist() + [int(timeID)], dtype='int64')
        self.output_times = \
            na.array(self.output_times.tolist() + [float(time)], dtype='float64')
        self.sort()

    def get_before(self, time):
        return na.where(self.output_times <= time)[0]

    def get_after(self, time):
        return na.where(self.output_times > time)[0]

    def __delitem__(self, key):
        if isinstance(key, types.StringType):
            id = self.index(key)
        elif isinstance(key, types.IntType):
            id = key
        self.output_names = na.array(self.output_names[:id].tolist() \
                                  + self.output_names[id+1:].tolist())
        self.output_times = na.array(self.output_times[:id].tolist() \
                                  + self.output_times[id+1:].tolist())
        self.output_time_ids = na.array(self.output_time_ids[:id].tolist() \
                                  + self.output_time_ids[id+1:].tolist())

    def __getitem__(self, key):
        """
        Based on the type of input, we either return based on index or
        basename.
        """
        if isinstance(key, types.StringType):
            index = self.index(key)
        elif isinstance(key, types.IntType):
            index = key
        # This fails, but I don't know how to fix it.
        a = (self.output_names[index],  \
             self.output_times[index],  \
             self.output_time_ids[index] )
        return a

    def index(self, key):
        t = os.path.basename(key)
        # Find out the index
        index = None
        for i in range(self.output_names.shape[0]):
            if os.path.basename(self.output_names[i]) \
                == os.path.basename(t):
                index = i
                break
        if index == None:
            raise KeyError
        return index

    def __repr__(self):
        return self.title

    def __len__(self):
        return len(self.output_names)

    def __iter__(self):
        for i in range(len(self)):
            yield self.output_names[i]

    def keys(self):
        return self.output_names.tolist()

    def __convert_args(self, args, kwargs, my_dict):
        new_args = list(args)[:]
        new_kwargs = {}
        for j, arg in enumerate(args):
            if isinstance(arg, types.StringType):
                new_args[j] = arg % myDict
        for key in kwargs.keys():
            if isinstance(kwargs[key], types.StringType):
                new_kwargs[key] = kwargs[key] % myDict
        return new_args, new_kwargs

    def run_functions(self, function, args = None, kwargs = None):
        if args is None: args = []
        if kwargs is None: kwargs = {}
        import yt.lagos as lagos # We import *here* so that we only import if we need it
        for i,o in enumerate(self):
            # Now we format string the various args and kwargs
            myDict = {'fn':str(o), 'index':i, 'time':self.output_times[i],
                      'timeID':self.output_time_ids[i]}
            newa, newk = self.__convert_args(args, kwargs)
            print args, kwargs
            function(*newa, **newk)

def GrabCollections(path=None):
    if not path: path=ytcfg.get("fido","rundir")
    if not os.path.isdir(path): os.makedirs(path)
    ocs = []
    for file in glob.glob(os.path.join(path,"runF_*")):
        title=os.path.basename(file)
        if title.startswith("runF_"): title = title[5:]
        ocs.append(OutputCollection(title))
        ocs[-1].read_in(file)
    return ocs
