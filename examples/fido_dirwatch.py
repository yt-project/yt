import yt.fido as fido
import os, os.path

image_skel = "ravenimages/%s/"
image_prefix = "%(bn)s_%(width)05i_%(unit)s"

def fido_plot_wrap():
    import yt.raven as raven
    import yt.lagos as lagos
    def fido_plots(fn):
        pf = lagos.EnzoStaticOutput(fn)
        md = pf["MetaDataString"]
        image_dir = image_skel % (md)
        if not os.path.isdir(image_dir):
            os.makedirs(image_dir)
        prefix = os.path.join(image_dir, image_prefix)
        raven.MakePlots(pf, "automated_plotconfig.xml", prefix, -1)
        return
    return fido_plots

Giles = fido.Watcher(function_handler=fido_plot_wrap)
Giles.run()
