import yt.raven as raven
import yt.fido as fido
import os, os.path

imageSkel = "/u/ki/mturk/public_html/ravenimages/%s/"
imgPrefix = "%(bn)s_%(width)05i_%(unit)s"

def fidoPlots(h):
    md = h["MetaDataString"]
    imageDir = imageSkel % (md)
    if not os.path.isdir(imageDir):
        os.makedirs(imageDir)
    prefix = os.path.join(imageDir, imgPrefix)
    plot=raven.EnzoHippo(h)
    raven.MakePlots(plot, "sample_raven_config.xml", prefix)
    return

fido.watchDir(funcHandler=fidoPlots)

