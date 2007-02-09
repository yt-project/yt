import yt.lagos as lagos
import yt.raven as raven
import os, os.path
import yt.deliverator as deliverator_upload

fields = ["NumberDensity", "Temperature", "H2I_Fraction"]
imageSkel = "/u/ki/mturk/public_html/ravenimages/%s/slices/"

def makePlots(path, bn):
    h=lagos.EnzoHierarchy(os.path.join(path, bn))
    #h.rates = chemistry.EnzoTable(os.path.join(path,"rates.out"),rates_out_key)
    #h.cool = chemistry.EnzoTable(os.path.join(path,"cool_rates.out"),cool_out_key)
    user=os.getenv("USER")
    if user == "matthewturk":
        user = "mturk"
    #RunID, Response = deliverator_upload.SubmitRun(h["MetaDataString"], user)
    #print "Submitted '%s' as user '%s' and got %s and '%s'" % \
            #(h["MetaDataString"],user, RunID, Response)
    md = h["MetaDataString"]
    imageDir = imageSkel % (md)
    if not os.path.isdir(imageDir):
        os.makedirs(imageDir)
    dx = h.getSmallestDx()
    for axis in range(3):
        plot=raven.EnzoHippo(h)#, submitToDeliverator=RunID)
        plot.canvas.setPlotMatrix(1,1)
        i = 0
        for field in fields:
            if i == 0:
                p = plot.addSlice(field, axis)
            else:
                p[0].switchField(field)
            for w in [10000, 1000, 100, 10, 1]:#, 0.1, 0.01, 0.001]:
                print w/h["au"], dx
                if (w/h["au"] < dx):
                    continue
                plot.setWidth(w,"au")
                fn = plot.saveImages(os.path.join(imageDir,"slices/%s_%010iau" % (bn,w)),"png")
                print fn
            for w in [10000, 1000, 100, 10, 1]:#, 0.1, 0.01, 0.001]:
                print w/h["rsun"], dx
                if (w/h["rsun"] < dx):
                    continue
                plot.setWidth(w,"rsun")
                fn = plot.saveImages(os.path.join(imageDir,"slices/%s_%010irsun" % (bn,w)),"png")
                print fn
            i += 1
        del plot
