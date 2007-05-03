# This uses EnzoRun to make some plots

import yt.lagos as lagos
import yt.fido as fido
import yt.raven as raven
import yt.deliverator as deliverator

import os, os.path, hippo

# These will probably end up being in .yt/config

imageSkel = "/u/ki/mturk/public_html/ravenimages/%s/"
httpPrefix = "http://www.slac.stanford.edu/~mturk/ravenimages/%s/"

def makePlots(h, string):
    # Now that we fork, we need to create the HDApp in each function call
    app = hippo.HDApp( )
    bn=h.basename
    user=os.getenv("USER")
    if user == "matthewturk":
        user = "mturk"
    md = h["MetaDataString"]
    RunID, Response = deliverator.SubmitRun(md, user)
    print "Submitted '%s' as user '%s' and got %s and '%s'" % \
            (md,user, RunID, Response)
    imageDir = imageSkel % (md)
    if not os.path.isdir(imageDir):
        os.makedirs(imageDir)
    dx = h.getSmallestDx()
    plot=raven.EnzoHippo(h, submitToDeliverator=RunID, httpPrefix=httpPrefix % (md), app=app)
    plot.canvas.setPlotMatrix(1,1)
    plots = plot.addNewProj("Density", weight="Temperature")
    for w in [1000, 100, 10, 1]:
        if (w/h["pc"] < 10*dx):
            continue
        plot.setWidth(w,"pc")
        fn = plot.saveImages(os.path.join(imageDir,"%s_%010ipc" % (bn,w)),"png")
    for w in [10000, 1000, 100, 10, 1]:
        if (w/h["au"] < 10*dx):
            continue
        plot.setWidth(w,"au")
        fn = plot.saveImages(os.path.join(imageDir,"%s_%010iau" % (bn,w)),"png")
    for w in [10000, 1000, 100, 10, 1]:
        if (w/h["au"] < 10*dx):
            continue
        p = plot.addTwoPhase(["NumberDensity","Temperature"], w,"au")
        plot.saveImages(os.path.join(imageDir,"%s_%010iau" % (bn,w)),"png",-1)
    return

# Note that we assume an existing run is already around...  And there's
# currently a slight disconnect between the directory name and the metadata.
# I aim to fix this...

myRun = fido.fetchRun(os.path.basename(os.getcwd()))
myRun.runFunction(makePlots, "%04i")
