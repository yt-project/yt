"""
Some functions to include the output from enzo_anyl, and generate it if need
be.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@todo: Create an AnalyzeCluster file from a given set of parameters, and a
hierarchy file.
@todo: Run enzo_anyl and parse the resultant data.
"""

from yt.lagos import *
import types, exceptions

ACOtranslation = {
        'r' : 'bin central radius (Mpc)',
        'R' : 'bin outer radius (Mpc)',
        'rho' : 'd_gas (Ms/Mpc^3)',
        'v_rms' : 'v_rms_gas_3d (km/s)',
        'T' : 'temp_gas_mass_weighted (K)',
        'v_r' : 'vr_gas (km/s)',
        'v_r_rms' : 'vr_rms_gas (km/s)',
        'T_cool' : 'T_cool (s)',
        'T_dyn' : 'T_dyn (s)',
        'm_enc' : 'm_gas (Ms)',
        'mu' : 'mu (mh) mean mass per particle'
                 }

CFfieldInfo = {}
CFfieldInfo["T"] = (None, r"$T (\rm{K})$", False)
CFfieldInfo["m_enc"] = (None, r"$M_{\rm{enc}} (M_\odot)$", True)
CFfieldInfo["v_r"] = (None, r"$v_r (km/s)$", True)
CFfieldInfo["H2I fraction"] = (None, r"$\rho({\rm{H}}_2) / \rho$", True)

class AnalyzeClusterOutput:
    # This is a class for storing the results of enzo_anyl runs
    def __init__(self, filename, id=None):
        """
        We receive a filename, and then parse it to get the key of columns, and
        then parse the columns to get the array of values
        """
        self.filename = filename
        if isinstance(filename, types.StringType):  # For the __add__ method
            self.lines = open(self.filename).readlines()
            self.parseKey()
            self.parseData()
        elif isinstance(filename, types.DictType):
            self.filename = id
            ks = filename.keys()
            ks.sort()
            numBins = filename[ks[0]].shape[0]
            self.data = na.zeros((len(ks), numBins), nT.Float64)
            self.columns = {}
            self.rcolumns = {}
            i = 0
            for column in ks:
                #print column, self.data.shape, self.data[i,:].shape, filename[column].shape
                self.data[i,:] = filename[column]
                self.columns[i] = column
                self.rcolumns[column] = i
                i += 1

    def __getitem__(self, item):
        if isinstance(item, types.StringType):
            if self.rcolumns.has_key(item):
                i = self.rcolumns[item]
            elif ACOtranslation.has_key(item):
                return self[ACOtranslation[item]]
            else:
                return self.generateField(item)
            #print item, i
        else:
            i = int(item)
        return self.data[i]

    def generateField(self, item):
        if CFfieldInfo.has_key(item):
            return CFfieldInfo[item][0](self)
        raise KeyError

    def parseKey(self):
        """
        Scans, grabs the column names
        """
        self.columns = {}
        self.rcolumns = {}
        index = self.lines.index("# COLUMN   DESCRIPTION\n")
        for line in self.lines[index+1:]:
            if line.rstrip() == "#":
                break
            col, name = line[1:].strip().split(" ",1)
            col = int(col) - 1
            self.columns[col] = name.strip().rstrip()
            self.rcolumns[name.strip().rstrip()] = col

    def parseData(self):
        """
        Parse the data into self.data

        @todo: Convert to use a regular expression and counting the number of
        matches
        """
        numBins = 0
        for line in self.lines:
            if line[0] != "#" and len(line.strip()) > 0:
                numBins += 1
        self.data = na.zeros((len(self.columns), numBins), nT.Float64)
        i = 0
        nc = len(self.columns)
        for line in self.lines:
            if line[0] != "#" and len(line.strip()) > 0:
                self.data[:,i] = map(float, line.split()[:nc])
                i += 1

    def outputHippo(self, filename):
        """
        This outputs the columns in a format that HippoDraw can read easily
        """
        fs = "\t".join(["%0.15e"] * len(self.columns)) + "\n"
        f = open(filename,"w")
        # Dicts are unsorted, so we do this.  Could probably be better.
        fields = [] 
        keys = self.columns.keys()
        keys.sort()
        for key in keys:
            fields.append(self.columns[key])
        header = "\t".join(fields) + "\n"
        f.write(header)
        for i in xrange(self.data.shape[1]):
            f.write(fs % tuple(self.data[:,i]))
        f.close()

    def __add__(self, other):
        """
        This will add all non-duplicated columns.

        @note: We kind of have to do this the slow way, to ensure that no
        double-columning goes on.  I couldn't come up with a better way than
        this.  Maybe I will, some day.  But for now, it works!  Hooray!
        """
        comb = AnalyzeClusterOutput(None)
        comb.columns = self.columns.copy()
        comb.rcolumns = self.rcolumns.copy()
        cols = other.columns.keys()
        offset = len(comb.columns) # 0-indexed
        i = offset
        for col in cols:
            if other.columns[col] not in self.columns.values():
                comb.columns[i] = other.columns[col]
                comb.rcolumns[other.columns[col]] = i
                i += 1
        comb.data = na.zeros((len(comb.columns), self.data.shape[1]), nT.Float64)
        comb.data[:len(self.columns),:] = self.data
        i = offset
        for col in cols:
            if other.columns[col] not in self.columns.values():
                comb.data[i,:] = other.data[col,:] # data is zero-indexed
                i += 1
        #print comb.data[offset:,:].shape, other.data[1:,:].shape
        #comb.data[offset:,:] = other.data[1:,:]
        return comb

def cfCS(self):
    # Assume gamma = 5./3.
    k = 1.380e-16 # ergs / K
    cs = na.sqrt((5./3. * k * self["T"] * self["NumberDensity"] \
            / self["RhoCGS"])) \
            / 1.0e5
    return cs
CFfieldInfo["cs"] = (cfCS, r"$c_s (\rm{km / s})$", False)

def cfMu(self):
    # This is in case mu craps out in enzo_anyl
    # mu = mass / num particles
    # mu = density / number density
    return (self["RhoCGS"] / self["NumberDensity"]) / 1.67e-24
CFfieldInfo["Mu"] = (cfMu, r"$\mu$", False)

def cfMach(self):
    return abs(self["v_r"]/self["cs"])
CFfieldInfo["Mach"] = (cfMach, r"$\rm{Radial Mach Number}$", False)

def cfTurbMach(self):
    return abs(self["vr_rms_gas (km/s)"]/self["cs"])
CFfieldInfo["TurbMach"] = (cfTurbMach, r"$\rm{Turbulent Mach Number}$", False)

def cfRAU(self):
    return self["r"] * unitList["au"]
CFfieldInfo["RAU"] = (cfRAU, r"$R (\rm{AU})$", True)

def cfRpc(self):
    return self["r"] * unitList["pc"]
CFfieldInfo["Rpc"] = (cfRpc, r"$R (\rm{pc})$", True)

def cfRhoCGS(self):
    return self["rho"] * 6.769e-41
CFfieldInfo["RhoCGS"] = (cfRhoCGS, r"$\rho (\rm{g} \rm{cm}^{-3})$", True)

def cfCoolDynComp(self):
    return self["T_cool"]/self["T_dyn"]
CFfieldInfo["CoolDynComp"] = (cfCoolDynComp, r"$T_{\rm{cool}} / T_{\rm{dyn}}$", True)

def cfNumberDensity(self):
    t= (self["RhoCGS"] / mh) * (\
           (self["HI fraction"] + self["HII fraction"] + self["H- fraction"]) / 1.0 \
         + (self["HeI fraction"] + self["HeII fraction"] + self["HeIII fraction"]) / 4.0 \
         + (self["H2I fraction"] + self["H2II fraction"]) / 2.0 \
         + (self["e- fraction"])) #\
         #+ (self["DI fraction"] + self["DII fraction"])/2.0 + self["HDI fraction"]/3.0 )
    #return (self["RhoCGS"]/mh) * self["H2I fraction"] / 2.0
    return t
CFfieldInfo["NumberDensity"] = (cfNumberDensity, r"$n (\rm{cm}^{-3})$", True)
