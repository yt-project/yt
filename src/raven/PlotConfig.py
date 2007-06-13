"""
Here we define a method for parsing a plot-descriptor configuration file, and
then using that to create and save out the plots.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

from xml.etree.ElementTree import ElementTree as et

def MakePlots(eh, fn, prefix):
    """
    This is the simplest way I can think of to do this.  We can get basically
    everything we need from the EnzoHippo instance.

    @param eh: The EnzoHippo instance we're going to use
    @param fn: the XML file we're going to interpret to get
                           plot-info
    @todo: Add "center" argument
    """
    pi = et(file=fn)
    r = pi.getroot()
    h = eh.hierarchy
    bn = h.basename
    dx = h.getSmallestDx()
    for plot in r.getchildren():
        # Now we have a plot type, so let's get going, peoples!
        if plot.attrib.has_key("mindx"):
            mindx = float(plot.attrib["mindx"])
        else:
            mindx = 10 # 10 cells by default
        fields = []
        widths = []
        for fe in plot.findall("field"):
            fields.append(fe.text)
        for wi in plot.findall("width"):
            #print wi.text, wi.attrib["unit"]
            widths.append((float(wi.text), \
                           wi.attrib["unit"]))
        #print "widths:", widths
        plotType = plot.tag.lower()
        prefixDict = {'bn':bn, 'type':plotType}
        # Note that we never set fields in here, as the saveImages takes care
        # of that.
        if plotType == "slice":
            f = 0
            center = None
            if plot.attrib.has_key("center"):
                center = map(float(plot.attrib["center"].split()))
            for field in fields:
                if f == 0:
                    pis = eh.addSlice(field, center=center)
                    pisl = range(-1*len(pis),0)
                else:
                    for pi in pis:
                        pi.switchField(field)
                f += 1
                for width, unit in widths:
                    if (width/h[unit] < mindx*dx):
                        #print "char:", width/h[unit]
                        #print "morechar:",  mindx*dx, mindx, dx
                        continue
                    eh.setWidth(width, unit, pisl)
                    prefixDict['width'] = width
                    prefixDict['unit'] = unit
                    eh.saveImages(prefix % prefixDict, "png", pisl)
        if plotType == "threephase":
            for width, unit in widths:
                if (width/h[unit] < mindx*dx):
                    continue
                prefixDict['width'] = width
                prefixDict['unit'] = unit
                eh.addThreePhase(fields, width, unit)
                eh.saveImages(prefix % prefixDict, "png", -1)
        elif plotType == "twophase":
            for width, unit in widths:
                if (width/h[unit] < mindx*dx):
                    continue
                prefixDict['width'] = width
                prefixDict['unit'] = unit
                eh.addTwoPhase(fields, width, unit)
                eh.saveImages(prefix % prefixDict, "png", -1)
        elif plotType == "proj":
            for fe in plot.findall("field"):
                field = fe.text
                weight = None
                if fe.attrib.has_key("weight"):
                    weight = fe.attrib["weight"]
                pis = eh.addNewProj(field, weight=weight)
                pisl = range(-1*len(pis),0)
                for width, unit in widths:
                    if (width/h[unit] < mindx*dx):
                        continue
                    prefixDict['width'] = width
                    prefixDict['unit'] = unit
                    eh.setWidth(width, unit, pisl)
                    eh.saveImages(prefix % prefixDict, "png", pisl)
