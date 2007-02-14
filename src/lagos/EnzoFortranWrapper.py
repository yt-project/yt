#
# src_wrapper:
#   Using pyFort, we can make fortran calls
#   Useful for testing chemistry, getting cooling time, etc
#
# Written by: Matthew Turk (mturk@stanford.edu) Nov 2006
# Modified:
#

from yt.lagos import *
import Numeric # Hate doing this, but we have to for inout ability
from ravenDefs import *
from numarray import *

import enzo_routines

# So I will first write a wrapper for the solve_rate_cool function
# This will take an actual grid and its actual data, and then get some results
# back.
#
#      subroutine solve_rate_cool(d, e, ge, u, v, w, de,
#     &                HI, HII, HeI, HeII, HeIII,
#     &                in, jn, kn, nratec, iexpand, imethod,
#     &                idual, ispecies, imetal, idim,
#     &                is, js, ks, ie, je, ke, ih2co, ipiht,
#     &                dt, aye, temstart, temend, 
#     &                utem, uxyz, uaye, urho, utim, uvel,
#     &                eta1, eta2, gamma, fh, dtoh,
#     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
#     &                k11a, k12a, k13a, k13dda, k14a, k15a,
#     &                k16a, k17a, k18a, k19a, k21a, k22a, k23a,
#     &                k24, k25, k26, k27, k28, k29, k30, k31,
#     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
#     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
#     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
#     &                reHeII2a, reHeIIIa, brema, compa,
#     &                comp_xraya, comp_temp, piHI, piHeI, piHeII,
#     &                HM, H2I, H2II, DI, DII, HDI, metal,
#     &                hyd01ka, h2k01a, vibha, rotha, rotla, 
#     &                gpldla, gphdla, hdltea, hdlowa, hdcoola, ciecoa, 
#     &                inutot, iradtype, nfreq, imetalregen,
#     &                iradshield, avgsighp, avgsighep, avgsighe2p,
#     &                iciecool, ih2optical, errcode, omaskflag, subgridmask )
#

def runSolveRateCool(g, dt, omaskflag=0):
    a = g.hierarchy
    # First we will make the grid read all the data in
    #print "Reading all data and feeding to solve_rate_cool"
    g.readAllData()
    # Now we copy, so that we can transpose it to row-major order
    dataCopy = {}
    for ds in g.data.keys():
        try:
            t = Numeric.array(g.data[ds].copy(), Numeric.Float32)
            dataCopy[ds] = t
        except TypeError:
            print "Skipping %s" % (ds)
            del g.data[ds]
        #t.transpose()
        #dataCopy[ds] = t
    # Let's get all the rates
    #print "Getting chemistry rates from rate file"
    for rate in rates_out_key:
        #print rate,
        exec("%s = a.rates['%s']" % (rate, rate))
    for rate in a.rates.params.keys():
        exec("%s = a.rates.params['%s']" % (rate, rate))
    #print "\n\nGetting cooling rates from cooling file"
    for rate in cool_out_key:
        #print rate,
        exec("%s = a.cool['%s']" % (rate, rate))
    for rate in a.cool.params.keys():
        exec("%s = a.cool.params['%s']" % (rate, rate))
    utim = 2.52e17 / sqrt(a.parameters["CosmologyOmegaMatterNow"]) \
                   / a.parameters["CosmologyHubbleConstantNow"] \
                   / (1+a.parameters["CosmologyInitialRedshift"])**1.5
    urho = 1.88e-29 * a.parameters["CosmologyOmegaMatterNow"] \
                    * a.parameters["CosmologyHubbleConstantNow"]**2 \
                    * (1.0 + a["CosmologyCurrentRedshift"])**3
    uxyz = 3.086e24 * \
           a.parameters["CosmologyComovingBoxSize"] / \
           a.parameters["CosmologyHubbleConstantNow"] / \
           (1.0 + a.parameters["CosmologyCurrentRedshift"])
    uaye = 1.0/(1.0 + a.parameters["CosmologyInitialRedshift"])
    uvel = 1.225e7*a.parameters["CosmologyComovingBoxSize"] \
                  *sqrt(a.parameters["CosmologyOmegaMatterNow"]) \
                  *sqrt(1+ a.parameters["CosmologyInitialRedshift"])
    #uvel = a["x-velocity"]
    utem = 1.88e6 * (a.parameters["CosmologyComovingBoxSize"]**2) \
                  * a.parameters["CosmologyOmegaMatterNow"] \
                  * (1.0 + a.parameters["CosmologyInitialRedshift"])
    aye  = (1.0 + a.parameters["CosmologyInitialRedshift"]) / \
           (a.parameters["CosmologyCurrentRedshift"] - 1.0)
    # Now we have all the units!  We're almost done...
    blank_field = Numeric.zeros(g.data["Total_Energy"].shape, Numeric.Float32)
    hdc = array([hdc_1, hdc_2, hdc_3, hdc_4, hdc_5], Float32)
    hdc.transpose()
    k13dd = array([k13_1, k13_2, k13_3, k13_4, k13_5, k13_6, k13_7], Float32)
    k13dd.transpose()
    inutot = array([0, 0, 1, 0], Float32)
    inutot.transpose()
    comp_xray = 0
    comp_temp = 0
    errcode = 0
    g.generateChildMask()
    #subgridmask = Numeric.ones(g.data["Density"].shape, Numeric.Int32) * array(g.myChildMask, Int32)
    #subgridmask = Numeric.array(g.myChildMask, Numeric.Int64)
    subgridmask = g.myChildMask * Numeric.ones(g.ActiveDimensions, Numeric.Int32)
    #subgridmask = Numeric.reshape(Numeric.transpose(subgridmask), subgridmask.shape)
    print subgridmask.shape
    #return subgridmask
    enzo_routines.solve_rate_cool( \
        dataCopy["Density"], dataCopy["Total_Energy"], dataCopy["Gas_Energy"],
        dataCopy["x-velocity"], dataCopy["y-velocity"], dataCopy["z-velocity"], 
        dataCopy["Electron_Density"], dataCopy["HI_Density"], dataCopy["HII_Density"],
        dataCopy["HeI_Density"], dataCopy["HeII_Density"], dataCopy["HeIII_Density"],
        g.ActiveDimensions[0], g.ActiveDimensions[1], g.ActiveDimensions[2],
        len(tgas), a.parameters["ComovingCoordinates"], a.parameters["HydroMethod"],
        a.parameters["DualEnergyFormalism"], a.parameters["MultiSpecies"],
        0, 3, 0, 0, 0,
        g.ActiveDimensions[0]-1, g.ActiveDimensions[1]-1, g.ActiveDimensions[2]-1,
        1, 1, dt, aye, tgas[0], tgas[-1],
        utem, uxyz, uaye, urho, utim, uvel,
        a.parameters["DualEnergyFormalismEta1"], a.parameters["DualEnergyFormalismEta2"], 
        a.parameters["Gamma"], 0.76, 2.0*3.4e-5,
        k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
        k11, k12, k13, k13dd, k14, k15,
        k16, k17, k18, k19, k21, k22, k23,
        k24, k25, k26, k27, k28, k29, k30, k31,
        k50, k51, k52, k53, k54, k55, k56,
        ceHI, ceHeI, ceHeII, ciHI, ciHeI, 
        ciHeIS, ciHeII, reHII, reHeII1, 
        reHeII2, reHeIII, brem, comp,
        comp_xray, comp_temp, piHI, piHeI, piHeII,
        dataCopy["HM_Density"], dataCopy["H2I_Density"], dataCopy["H2II_Density"],
        blank_field, blank_field, blank_field, blank_field,
        hyd01k, h2k01, vibh, roth, rotl, 
        gpldl, gphdl, hdlte, hdlow, hdc, cieco, 
        inutot[0], int(inutot[1]), int(inutot[2]), int(inutot[3]),
        0, 0, 0, 0,
        1, 1, errcode, omaskflag, subgridmask)
    # iciecool, ih2optical
    # Okay, now we're done.  Note the couple blank fields, especially the HD
    # ones!
    #print dataCopy["H2I_Density"] - g.data["H2I_Density"]#, dataCopy["H2I_Density"]
    #print dataCopy["HM_Density"] - g.data["HM_Density"]#, dataCopy["H2I_Density"]
    return dataCopy

#
# We call this to get the cooling time.  Hooray!
#
#      subroutine cool_multi_time(
#     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
#     &                  cooltime,
#     &                in, jn, kn, nratec, iexpand, imethod,
#     &                idual, ispecies, imetal, idim,
#     &                is, js, ks, ie, je, ke, ih2co, ipiht,
#     &                dt, aye, temstart, temend,
#     &                utem, uxyz, uaye, urho, utim,
#     &                eta1, eta2, gamma,
#     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa,
#     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a,
#     &                reHeII2a, reHeIIIa, brema, compa,
#     &                comp_xraya, comp_temp, piHI, piHeI, piHeII,
#     &                HM, H2I, H2II, DI, DII, HDI, metal,
#     &                hyd01ka, h2k01a, vibha, rotha, rotla,
#     &                gpldla, gphdla, hdltea, hdlowa, hdcoola, ciecoa,
#     &                inutot, iradtype, nfreq, imetalregen,
#     &                iradshield, avgsighp, avgsighep,
#     &                avgsighe2p, iciecool, ih2optical )
#

class Reshaper:
    def __init__(self, g):
        self.g = g
        self.data = {}
        # So let's figure out how big to make this thing
        # We have a maximum of NUMPERROW on a side
        
        
    def __getitem__(self, key):
        self.data[key] = array(reshape(self.g[key], self.g.ActiveDimensions), Numeric.Float32)
        return self.data[key]

NUMPERROW = 256

def runCoolMultiTime(g):
    a = g.hierarchy
    # First we will make the grid read all the data in
    #print "Reading all data from grid %s (%s, %s, %s) and feeding to cool1d_multi_bdf" % \
        #(g.id, g.ActiveDimensions[0], g.ActiveDimensions[1], g.ActiveDimensions[2])
    g.readAllData()
    dataCopy = {}
    # Now we copy, so that we can transpose it to row-major order
    for rate in cool_out_key:
        exec("%s = Numeric.array(a.cool['%s'])" % (rate, rate))
    for rate in a.cool.params.keys():
        exec("%s = a.cool.params['%s']" % (rate, rate))
    dataCopy = {}
    # Let's figure out what to reshape it to
    if max(g.ActiveDimensions[0], g.ActiveDimensions[1], g.ActiveDimensions[2]) > NUMPERROW:
        # Now we generate the new sizes
        numElements = g.ActiveDimensions[0] \
                    * g.ActiveDimensions[1] \
                    * g.ActiveDimensions[2]
        # We want to partition this into NUMPERROW size chunks
        AD = ones(3)
        #print "AD:", AD, numElements
        for i in range(3):
            #j = 2-i
            j = i
            AD[j] = min(NUMPERROW,numElements)
            #print AD
            numElements = ceil(numElements / float(AD[j]))
            print "AD:", AD, numElements
        #AD.reverse()
        #AD = array(AD)
    else:
        AD = g.ActiveDimensions
    for ds in g.data.keys():
        try:
            dataCopy[ds] = Numeric.resize(g.data[ds], tuple(AD.tolist()))
            #dataCopy[ds] = Numeric.reshape(Numeric.transpose(dataCopy[ds]), dataCopy[ds].shape)
            print "Copying %s" % (ds)
        except KeyError:
            print "Skipping %s" % (ds)
            del g.data[ds]
    utim = a.conversionFactors["Time"]
    urho = a.conversionFactors["Density"]
    uxyz = 3.086e24 * \
           a.parameters["CosmologyComovingBoxSize"] / \
           a.parameters["CosmologyHubbleConstantNow"] / \
           (1.0 + a.parameters["CosmologyCurrentRedshift"])
    uaye = 1.0/(1.0 + a.parameters["CosmologyInitialRedshift"])
    uvel = a.conversionFactors["x-velocity"]
    utem = a.conversionFactors["Temp"]
    aye  = (1.0 + a.parameters["CosmologyInitialRedshift"]) / \
           (a.parameters["CosmologyCurrentRedshift"] - 1.0)
    # Now we have all the units!  We're almost done...
    blank_field = zeros(AD, Numeric.Float32)
    hdc = Numeric.transpose(Numeric.array([hdc_1, hdc_2, hdc_3, hdc_4, hdc_5], Numeric.Float32))
    inutot = array([0, 0, 1, 0], Float32)
    inutot.transpose()
    comp_xray = 0
    comp_temp = 0
    dt = 0 # Doesn't matter
    #cooltime = Numeric.zeros(g.ActiveDimensions,Numeric.Float32)
    k = g.data["Density"].shape
    cooltime = Numeric.zeros(dataCopy["Density"].shape, Numeric.Float32)
    #cooltime = Numeric.resize(cooltime, AD) 
    #cooltime = Numeric.reshape(Numeric.transpose(cooltime), cooltime.shape)
    # cool1d_multi is one-d, and works on a slice at a time
    #h = Reshaper(g)
    print "Going to call it!"
    dc = Numeric.reshape(Numeric.transpose(dataCopy["Density"]), dataCopy["Density"].shape)
    #print dataCopy["Density"].shape, dc.shape
    enzo_routines.cool_multi_time( \
        dataCopy["Density"], dataCopy["Total_Energy"], dataCopy["Gas_Energy"],
        dataCopy["x-velocity"], dataCopy["y-velocity"], dataCopy["z-velocity"], 
        dataCopy["Electron_Density"], dataCopy["HI_Density"], dataCopy["HII_Density"],
        dataCopy["HeI_Density"], dataCopy["HeII_Density"], dataCopy["HeIII_Density"],
        cooltime,
        AD[0], AD[1], AD[2],
        len(tgas), a.parameters["ComovingCoordinates"], a.parameters["HydroMethod"],
        a.parameters["DualEnergyFormalism"], a.parameters["MultiSpecies"],
        0, 3, 0, 0, 0, 
        AD[0]-1, AD[1]-1, AD[2]-1,
        1, 1, dt, aye, tgas[0], tgas[-1],
        utem, uxyz, uaye, urho, utim,
        a.parameters["DualEnergyFormalismEta1"], a.parameters["DualEnergyFormalismEta2"], 
        a.parameters["Gamma"], 
        ceHI, ceHeI, ceHeII, ciHI, ciHeI, 
        ciHeIS, ciHeII, reHII, reHeII1, 
        reHeII2, reHeIII, brem, comp,
        comp_xray, comp_temp, piHI, piHeI, piHeII,
        dataCopy["HM_Density"], dataCopy["H2I_Density"], dataCopy["H2II_Density"],
        blank_field, blank_field, blank_field, blank_field, blank_field,
        hyd01k, h2k01, vibh, roth, rotl, 
        gpldl, gphdl, hdlte, hdlow, hdc, cieco, 
        inutot[0], int(inutot[1]), int(inutot[2]), int(inutot[3]),
        0, 0, 0, 0,
        1, 1)
    # On some platforms you need the following three lines:
    #print hdc, cieco
    print "SHAPE",hdc.shape,hdc[3000,:],cieco[3000]
    t = Numeric.transpose(cooltime)
    cooltime = Numeric.reshape(cooltime, t.shape)
    cooltime = Numeric.transpose(cooltime)
    cooltime = abs(cooltime)
    ct = resize(cooltime, k)
    #i = (where(ct == 0))[0].size()
    #j = (where(cooltime == 0))[0].size()
    #print "STUFF", i, j, ct.shape, cooltime.shape, AD
    #print "START"
    #print min(ct), max(ct)
    #print min(abs(ct)), max(abs(ct))
    #print "STOP"
    return abs(ct)
