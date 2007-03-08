import yt.lagos as lagos
import yt.raven as raven
import os, os.path, hippo
import yt.deliverator as deliverator_upload
import yt.fido

import logging

f=logging.getLogger("yt.fido")
s=logging.StreamHandler()
f.addHandler(s)

app = hippo.HDApp( )

fields = ["NumberDensity", "Temperature", "H2I_Fraction", "RadialVelocity"]
imageSkel = "/u/ki/mturk/public_html/ravenimages/%s/"
httpPrefix = "http://www.slac.stanford.edu/~mturk/ravenimages/%s/"

def makePlots(h):
    bn=h.basename
    #h=lagos.EnzoHierarchy(os.path.join(path, bn))
    #h.rates = chemistry.EnzoTable(os.path.join(path,"rates.out"),rates_out_key)
    #h.cool = chemistry.EnzoTable(os.path.join(path,"cool_rates.out"),cool_out_key)
    user=os.getenv("USER")
    if user == "matthewturk":
        user = "mturk"
    RunID, Response = deliverator_upload.SubmitRun(h["MetaDataString"], user)
    print "Submitted '%s' as user '%s' and got %s and '%s'" % \
            (h["MetaDataString"],user, RunID, Response)
    md = h["MetaDataString"]
    imageDir = imageSkel % (md)
    if not os.path.isdir(imageDir):
        os.makedirs(imageDir)
    dx = h.getSmallestDx()
    plot=raven.EnzoHippo(h, submitToDeliverator=RunID, httpPrefix=httpPrefix % (md), app=app)
    plot.canvas.setPlotMatrix(1,1)
    i = 0
    for field in fields:
        if i == 0:
            plots = plot.addSlice(field)
        else:
            for p in plots:
                p.switchField(field)
        for w in [10000, 1000, 100, 10, 1]:#, 0.1, 0.01, 0.001]:
            #print w/h["au"], dx
            if (w/h["au"] < 10*dx):
                continue
            plot.setWidth(w,"au")
            fn = plot.saveImages(os.path.join(imageDir,"%s_%010iau" % (bn,w)),"png")
            #print fn
        for w in [10000, 1000, 100, 10, 1]:#, 0.1, 0.01, 0.001]:
            #print w/h["rsun"], dx
            if (w/h["rsun"] < 10*dx):
                continue
            plot.setWidth(w,"rsun")
            fn = plot.saveImages(os.path.join(imageDir,"%s_%010irsun" % (bn,w)),"png")
            #print fn
        i += 1
    if (1000/h["au"] >= 10*dx):
        p = plot.addThreePhase(["NumberDensity","Temperature","H2I_Fraction"], 1000,"au")
        plot.saveImages(os.path.join(imageDir,"%s_%010iau" % (bn,1000)),"png",-1)
    if (100/h["au"] >= 10*dx):
        p = plot.addThreePhase(["NumberDensity","Temperature","H2I_Fraction"], 100,"au")
        plot.saveImages(os.path.join(imageDir,"%s_%010iau" % (bn,100)),"png",-1)
    if (1000/h["au"] >= 10*dx):
        p = plot.addRadialProfilePlot(["NumberDensity","H2I_Fraction"], 1000,"au")
        for field in fields:
            p.switchField(field)
            plot.saveImages(os.path.join(imageDir,"%s_%010iau" % (bn,1000)),"png",-1)


yt.fido.watchDir(funcHandler=makePlots)
