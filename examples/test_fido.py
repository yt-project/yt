import yt.lagos as lagos
import yt.raven as raven
import yt.deliverator as deliverator_upload
import yt.fido
import os, os.path, hippo

import logging

f=logging.getLogger("yt.fido")
s=logging.StreamHandler()
f.addHandler(s)

fields = ["NumberDensity", "Temperature", "H2I_Fraction", "RadialVelocity"]
imageSkel = "/u/ki/mturk/public_html/ravenimages/%s/"
httpPrefix = "http://www.slac.stanford.edu/~mturk/ravenimages/%s/"

MINLEVEL=15

def makePlots(h):
    # Now that we fork, we need to create the HDApp in each function call
    app = hippo.HDApp( )
    if h.maxLevel < MINLEVEL:
        return
    bn=h.basename
    user=os.getenv("USER")
    if user == "matthewturk":
        user = "mturk"
    md = h["MetaDataString"]
    RunID, Response = deliverator_upload.SubmitRun(md, user)
    print "Submitted '%s' as user '%s' and got %s and '%s'" % \
            (md,user, RunID, Response)
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
            if (w/h["au"] < 10*dx):
                continue
            plot.setWidth(w,"au")
            fn = plot.saveImages(os.path.join(imageDir,"%s_%010iau" % (bn,w)),"png")
            print fn
        for w in [10000, 1000, 100, 10, 1]:#, 0.1, 0.01, 0.001]:
            if (w/h["rsun"] < 10*dx):
                continue
            plot.setWidth(w,"rsun")
            fn = plot.saveImages(os.path.join(imageDir,"%s_%010irsun" % (bn,w)),"png")
        i += 1
    for w in [10000, 1000, 100, 10, 1]:
        if (w/h["au"] < 10*dx):
            continue
        p = plot.addThreePhase(["NumberDensity","Temperature","H2I_Fraction"], w,"au")
        plot.saveImages(os.path.join(imageDir,"%s_%010iau" % (bn,w)),"png",-1)
    if (1000/h["au"] >= 10*dx):
        p = plot.addRadialProfilePlot(["NumberDensity","H2I_Fraction"], 1000,"au")
        for field in fields:
            p.switchField(field)
            plot.saveImages(os.path.join(imageDir,"%s_%010iau" % (bn,1000)),"png",-1)
    return


yt.fido.watchDir(funcHandler=makePlots)
