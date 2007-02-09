#
# gallery:
#   A module for dealing with raven images
#   TOTALLY STAND ALONE
#
# Written by: Matthew Turk (mturk@stanford.edu) Jan 2007
# Modified:
#

import os, os.path, pickle, types, urllib

class RavenImageDB:
    def __init__(self, filename, RO=False, local=True):
        # We load in the whole file if it doesn't exist
        self.filename = filename
        if local:
            if os.path.exists(filename):
                self.gallery = pickle.load(open(filename))
            elif RO==True:
                print "File does not exist, and in RO-mode"
                return None
            else:
                self.gallery = [] # List of dicts
        else:
            u = urllib.urlopen(filename)
            try:
                self.gallery = pickle.load(u)
            except:
                print "Problem loading %s - %s" % (u, filename)
                return None
    def getCharacteristics(self):
        chars = []
        # There's probably a better way to do this
        for im in self.gallery:
            for key in im.keys():
                if key not in chars:
                    chars.append(key)
        return chars

    def getValues(self, char):
        vals = []
        for im in self.gallery:
            if im.has_key(char):
                if im[char] not in vals:
                    vals.append(im[char])
        return vals

    def selectValue(self, char):
        ims = []
        for im in self.gallery:
            append = 1
            for k in char.keys():
                val = char[k]
                if (not im.has_key(k)) or (str(im[k]) != val):
                    append = 0
            if append:
                ims.append(im)
        return ims

    def add(self, im):
        if not isinstance(im, types.DictType):
            raise TypeError
        self.gallery.append(im.copy())
        self.writeOut()

    
    def writeOut(self):
        pickle.dump(self.gallery, open(self.filename,"w"))

