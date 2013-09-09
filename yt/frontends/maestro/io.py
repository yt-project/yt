"""
Maestro data-file handling functions - identical to Orion frontend,
but with all instances of "orion" changed to "maestro".


Authors:
 * Matthew Turk 
 * J. S. Oishi 
 * Chris Malone 


"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
import numpy as np
from yt.utilities.io_handler import \
           BaseIOHandler

from definitions import \
    yt2maestroFieldsDict

class IOHandlerNative(BaseIOHandler):

    _data_style = "maestro"

    def modify(self, field):
        return field.swapaxes(0,2)

    def _read_data(self,grid,field):
        """
        reads packed multiFABs output by BoxLib in "NATIVE" format.

        """
        filen = os.path.expanduser(grid.filename[field])
        off = grid._offset[field]
        inFile = open(filen,'rb')
        inFile.seek(off)
        header = inFile.readline()
        header.strip()

        if grid._paranoid:
            mylog.warn("Maestro Native reader: Paranoid read mode.")
            headerRe = re.compile(maestro_FAB_header_pattern)
            bytesPerReal,endian,start,stop,centerType,nComponents = headerRe.search(header).groups()

            # we will build up a dtype string, starting with endian
            # check endianness (this code is ugly. fix?)
            bytesPerReal = int(bytesPerReal)
            if bytesPerReal == int(endian[0]):
                dtype = '<'
            elif bytesPerReal == int(endian[-1]):
                dtype = '>'
            else:
                raise ValueError("FAB header is neither big nor little endian. Perhaps the file is corrupt?")

            dtype += ('f%i'% bytesPerReal) #always a floating point

            # determine size of FAB
            start = np.array(map(int,start.split(',')))
            stop = np.array(map(int,stop.split(',')))

            gridSize = stop - start + 1

            error_count = 0
            if (start != grid.start).any():
                print "Paranoia Error: Cell_H and %s do not agree on grid start." %grid.filename
                error_count += 1
            if (stop != grid.stop).any():
                print "Paranoia Error: Cell_H and %s do not agree on grid stop." %grid.filename
                error_count += 1
            if (gridSize != grid.ActiveDimensions).any():
                print "Paranoia Error: Cell_H and %s do not agree on grid dimensions." %grid.filename
                error_count += 1
            if bytesPerReal != grid.hierarchy._bytesPerReal:
                print "Paranoia Error: Cell_H and %s do not agree on bytes per real number." %grid.filename
                error_count += 1
            if (bytesPerReal == grid.hierarchy._bytesPerReal and dtype != grid.hierarchy._dtype):
                print "Paranoia Error: Cell_H and %s do not agree on endianness." %grid.filename
                error_count += 1

            if error_count > 0:
                raise RunTimeError("Paranoia unveiled %i differences between Cell_H and %s." % (error_count, grid.filename))

        else:
            start = grid.start_index
            stop = grid.stop_index
            dtype = grid.hierarchy._dtype
            bytesPerReal = grid.hierarchy._bytesPerReal

        nElements = grid.ActiveDimensions.prod()

        # one field has nElements*bytesPerReal bytes and is located
        # nElements*bytesPerReal*field_index from the offset location
        if yt2maestroFieldsDict.has_key(field):
            fieldname = yt2maestroFieldsDict[field]
        else:
            fieldname = field
        field_index = grid.field_indexes[fieldname]
        inFile.seek(int(nElements*bytesPerReal*field_index),1)
        field = np.fromfile(inFile,count=nElements,dtype=dtype)
        field = field.reshape(grid.ActiveDimensions[::-1]).swapaxes(0,2)

        # we can/should also check against the max and min in the header file

        inFile.close()
        return field

