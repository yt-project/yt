"""
Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

This file is part of yt.

yt is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


"""
Here we define a method for parsing a plot-descriptor configuration file, and
then using that to create and save out the plots.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
"""

from xml.etree.ElementTree import ElementTree as et

from yt.raven import *

def MakePlots(h, fn, prefix, deliverator):
    """
    This is the simplest way I can think of to do this.  We can get basically
    everything we need from the EnzoHippo instance.

    @param h: hierarchy instance
    @param fn: the XML file we're going to interpret to get
                           plot-info
    @param prefix: image file prefix
    @todo: Add "center" argument
    """
    pi = et(file=fn)
    r = pi.getroot()
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
        # Note that we never set fields in here, as the save takes care
        # of that.
        if plotType == "slice":
            pc = PlotCollection(h, deliverator)
            f = 0
            center = None
            if plot.attrib.has_key("center"):
                center = map(float(plot.attrib["center"].split()))
            for field in fields:
                if f == 0:
                    pc.addSlice(field, 0, center=center)
                    pc.addSlice(field, 1, center=center)
                    pc.addSlice(field, 2, center=center)
                else:
                    for plot in pc.plots:
                        plot.switch_z(field)
                f += 1
                for width, unit in widths:
                    if (width/h[unit] < mindx*dx):
                        #print "char:", width/h[unit]
                        #print "morechar:",  mindx*dx, mindx, dx
                        continue
                    pc.set_width(width, unit)
                    prefixDict['width'] = width
                    prefixDict['unit'] = unit
                    pc.save(prefix % prefixDict, "png")
        elif plotType == "threephase":
            for width, unit in widths:
                pc = PlotCollection(h, deliverator)
                if (width/h[unit] < mindx*dx):
                    continue
                prefixDict['width'] = width
                prefixDict['unit'] = unit
                #print fields, width, unit
                pc.addThreePhaseSphere(width, unit, fields)
                pc.save(prefix % prefixDict, "png")
        elif plotType == "twophase":
            for width, unit in widths:
                pc = PlotCollection(h, deliverator)
                if (width/h[unit] < mindx*dx):
                    continue
                prefixDict['width'] = width
                prefixDict['unit'] = unit
                pc.addTwoPhaseSphere(width, unit, fields)
                pc.save(prefix % prefixDict, "png")
        elif plotType == "proj":
            pc = PlotCollection(h, deliverator)
            for fe in plot.findall("field"):
                field = fe.text
                weight = None
                if fe.attrib.has_key("weight"):
                    weight = fe.attrib["weight"]
                pc.addProjection(field, 0, weightField=weight)
                pc.addProjection(field, 1, weightField=weight)
                pc.addProjection(field, 2, weightField=weight)
                for width, unit in widths:
                    if (width/h[unit] < mindx*dx):
                        continue
                    prefixDict['width'] = width
                    prefixDict['unit'] = unit
                    pc.set_width(width, unit)
                    pc.save(prefix % prefixDict, "png")
