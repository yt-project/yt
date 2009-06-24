"""
Python kD Tree

Author: Michael Knight
Affiliation: ?
Note: Code is based on <http://sites.google.com/site/mikescoderama/Home/kd-tree-knn>,
with periodicity added and a few other cosmetic changes. There is no contact
infomation on that page, and this code's license is a bit uncertain.
Author: Stephen Skory <stephenskory@yahoo.com>
Affiliation: UCSD Physics/CASS
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Matthew Turk.  All Rights Reserved.

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

from yt.lagos import *

from bisect import insort

class Point:pass

class Node:

    def _printOut(self,depth = 0):
        """
        iteratively print out the kd tree nodes, depth really isn't a user-set
        parameter.
        """
        for i in range(0,depth):
            print "....",
        print "point count =", self.pointCount, "rect =", \
            self.hyperRect._nodeToString(), '\n' #"points =",self.points
        if (self.leftChild is not None):
            for i in range(0,depth):
                print "....",
            print "left: "
            self.leftChild._printOut(depth+1)
        if (self.rightChild is not None):
            for i in range(0,depth):
                print "....",
            print "right: "
            self.rightChild._printOut(depth+1)
    
    def _buildBoundingHyperRect(self,points):
        self.hyperRect = HyperRect()
        self.hyperRect._buildBoundingHyperRect(points)
 
def getFastDistance(a,b,period):
    """
    returns the square of the distance between points a and b,
    using periodic boundary conditions 'period'.
    """
    dim = len(b);
    total = 0;
    for i in range(0,dim):
        delta = min(abs(a[i] - b[i]), period[i] - abs(a[i] - b[i]));
        total = total + (delta *delta)
    return total
 
class Neighbors:
    
    def _addNeighbors(self,node,query,period):
        """
        for each point in this node, calculate the distance to the query point,
        if it's close enough add it to the list .points using insort
        """
        for i in range(0,node.pointCount):
            dist = getFastDistance(node.points[i].data,query,period)
            
            if (dist < self.minDistanceSquared):
                item = [dist,node.points[i]]
                insort(self.points,item)
                if (len(self.points) > self.k):
                    self.points = self.points[0:self.k]
 
        if (len(self.points) == self.k):
            self.minDistanceSquared = self.points[self.k-1][0]
        return;
 
class HyperRect:
    def _buildBoundingHyperRect(self,points):
        """
        find the extremities of the hypercube.
        """
        self.k = len(points[0].data)
        self.dims  = range(0,self.k) 
        high = points[0].data[:]
        low = points[0].data[:]
        for i in range(0,len(points)):
            for j in self.dims:
                point = points[i].data[j]
                if (high[j]  < point):
                    high[j] = point
                if (low[j] > point):
                    low[j] = point
        self.high = high
        self.low = low
        return 
 
    def _getWidestDimension(self):
        """
        get the widest dimension of the points in order to find the dimension
        to bisect.
        """
        widest =0
        widestDim =-1
        for i in self.dims:
            width = self.high[i] -  self.low[i]
            if (width > widest):
                widestDim =i
                widest = width
        self.widest = widest;
        self.widestDim = widestDim;    
        return self.widestDim
    
    def _getWidestDimensionWidth(self): # I don't know why this is here.
        return self.widest
    
    def _nodeToString(self):
        return "high =",self.high,"low =",self.low,
    
    def _getMinDistance(self,query,period):
        """
        find the minimum distance squared to a corner of the hypercube
        from the query point using periodicity.
        """
        total = 0.0
        for i in self.dims:
            delta = 0.0
            min_high = min(abs(query[i] - self.high[i]), \
                period[i] - abs(query[i] - self.high[i]))
            min_low = min(abs(query[i] - self.low[i]), \
                period[i] - abs(query[i] - self.low[i]))
            delta = min(min_low,min_high)
            total = total + (delta*delta)
        return total;

def buildKdHyperRectTree(points,rootMin=3):
    """
    Recursively build the kdTree, adding nodes as needed until all have no more
    than rootMin points. The final nodes are called leafs, which contain the
    point data.
    """
    if (points is None or len(points) ==0):
        return None
    n = Node() # make a new node
    n._buildBoundingHyperRect(points)  # build the hyper rect for these points
                        # this will find the top left and botom
                        # right of all given points.

    leaf = False
    # If the size of points is small enough, this node is a leaf
    if len(points) <= rootMin:
        leaf = True
    splitDim  = -1
 
    if (not leaf):
        # get the widest dimension to split n to maximize splitting affect
        splitDim = n.hyperRect._getWidestDimension()
        # do we have a bunch of children at the same point?
        if (n.hyperRect._getWidestDimensionWidth() == 0.0):
            left = True 
    #init the node
    n.pointCount = len(points)
    n.points = None
    n.leftChild = None
    n.rightChild = None
    n.points = None
 
    if (leaf or len(points)==0):
        n.points = points # we are a leaf so just store all points in the rect
    else:
        # sort by the best split dimension
        temp = []
        for index,p in enumerate(points):
            insort(temp,[p.data[splitDim],index])
        temp2 = points[:]
        for index,t in enumerate(temp):
            temp2[index] = points[t[1]]
        points = temp2[:]
        del temp, temp2
        #points.sort(key=lambda points: points.data[splitDim])
        median = len(points)/2     # get the median
        # and split left for smaller values in splitDim, right for larger
        n.leftChild = buildKdHyperRectTree(points[0:(median+1)], rootMin)
        if (median+1 < len(points)):
            n.rightChild = buildKdHyperRectTree(points[median+1:], rootMin)
    return n;
 
def getKNN(query,node, neighbors,distanceSquared,period):
    """
    Recursively walk the kd tree, limited by *distanceSquared* to the extrema of
    the hypercubes, only finding distances to the query point in leaf nodes.
    *neighbors* is a Neighbors object, and needs to be initialized bofore
    calling this. *period* is a list or array of the period for each dimension
    of the hypercube.
    """
    
    # test to see if the query point is inside this node
    # <= and >= on both ends is okay, it's better to be inclusive and it
    # prevents problems with particles on boundaries
    for i in node.hyperRect.dims:
        if query[i] <= node.hyperRect.high[i] and \
                query[i] >= node.hyperRect.low[i]:
            inside = True
        else:
            inside = False
            break
    
    # if this node is close enough (the distances are calculated in the previous
    # iteration), or if the query point is inside the node, continue on
    if (neighbors.minDistanceSquared > distanceSquared) or inside:
        # leafs don't have children, so this tests to see if this node
        # is a leaf, and if it is a leaf, calculate distances to the query point
        if (node.leftChild is None):
            # add to neighbors.points
            neighbors._addNeighbors(node,query,period)
        # if this node is not a leaf, find out the distance to its children,
        # and then continue the iteration down the kd tree.
        else:
            distLeft = node.leftChild.hyperRect._getMinDistance(query,period)
            distRight = node.rightChild.hyperRect._getMinDistance(query,period)
            if (distLeft < distRight):
                getKNN(query,node.leftChild,neighbors,distLeft,period)
                getKNN(query,node.rightChild,neighbors,distRight,period)
            else:
                getKNN(query,node.rightChild,neighbors,distRight,period)
                getKNN(query,node.leftChild,neighbors,distLeft,period)
