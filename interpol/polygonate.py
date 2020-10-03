# -*- coding: utf-8 -*-
"""
Module to turn a set of points into a set of nonoverlapping polygons.
Done by identifying the removable edges, and then removing one of them.
Rinse-repeat until no further walls can be removed.

Not optimized or anything. For example, all candidiates are recalculated
after removing a wall, and the Delaunay grid is calculated several times.

2020-10
rwijtvliet@gmail.com
"""

import numpy as np
from scipy.spatial import Delaunay, ConvexHull
from typing import Iterable


class Polygonate:
    """
    Turn a set of points into a set of polygons.
    
    Arguments:
        points: Iterable of (x, y)-points
        pickwall: 
            'long'to remove the longest walls first;
            'acute' to remove the most acute angles first (default);
            'round' to remove wall to create polygon with roundest corners.
        convex: True if resulting polygons must be convex (default).
    """
    
    def __init__(self, points:Iterable, pickwall:str='', convex:bool=True):
        self.__points = np.array(points)
        self.__convex = convex
        self.__shapes, self.__neighbors_of_shapes = \
            self.__polygonation(pickwall)
        self.__descendent = self.__find_descedent()
    
    @property
    def points(self):
        """The (x, y) coordinates of the points."""
        return self.__points
    
    @property
    def vertices(self):
        """The point-indices of the vertices of each shape."""
        return self.__shapes
    
    @property
    def shapes(self):        
        """The point-indices of the vertices of each shape."""
        return self.__shapes
    
    @property
    def neighbors_of_shapes(self):
        return self.__neighbors_of_shapes
    
    def __polygonation(self, pickwall):
        self.__delaunay = Delaunay(self.__points)
        shapes = self.__delaunay.simplices.tolist()
        neighbors_of_shapes = [[si for si in neighbors if si != -1] 
                               for neighbors in self.__delaunay.neighbors]
        
        if pickwall.startswith('long'):
            pickwallfunc = lambda cands: np.argmax([cand['len'] 
                                                    for cand in cands])
        elif pickwall.startswith('round'):
            pickwallfunc = lambda cands: np.argmin([cand['error'] 
                                                    for cand in cands])
        else:
            pickwallfunc = lambda cands: np.argmin([cand['angles_before'].min()
                                                    for cand in cands])
    
        def melt(si1, si2, shape3): #remove shapes with indices si1 and si2. Add shape with vertices shape3.
            nonlocal shapes, neighbors_of_shapes
            if si1 > si2: si1, si2 = si2, si1
            shapes.pop(si2)
            shapes.pop(si1)  
            si3 = len(shapes)
            shapes.append(shape3)
            nei3 = [*neighbors_of_shapes.pop(si2), *neighbors_of_shapes.pop(si1)]
            nei3 = [si for si in nei3 if si != si1 and si != si2]
            neighbors_of_shapes.append(nei3)
            def new_si(si):
                if si==si1 or si==si2: return si3
                if si<si1: return si
                if si<si2: return si-1
                return si-2
            neighbors_of_shapes = [[new_si(si) for si in neighbors] 
                                   for neighbors in neighbors_of_shapes]
        
        while True:
            cands = self._candidates(shapes, neighbors_of_shapes)
            if len(cands) == 0: break
            # Find which one to remove.
            picked = cands[pickwallfunc(cands)]
            melt(*picked['si'], picked['shape3'])
        
        return shapes, neighbors_of_shapes
    
    @staticmethod
    def point_inside_polygon(x, y, poly):
        n = len(poly)
        inside = False
        p1x,p1y = poly[0]
        for i in range(n+1):
            p2x,p2y = poly[i % n]
            if y > min(p1y,p2y):
                if y <= max(p1y,p2y):
                    if x <= max(p1x,p2x):
                        if p1y != p2y:
                            xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x,p1y = p2x,p2y
        return inside
    
    def __find_descedent(self): #For each delaunay simplex, find shape it went into.
        descendent = []
        for sim in self.__delaunay.simplices:
            for s, shape in enumerate(self.__shapes):
                if len(np.intersect1d(sim, shape)) == len(sim):
                    descendent.append(s)
                    break
            else:
                raise ValueError(f'Cannot find child shape of simplex {sim}.')
        return descendent
        
    def find_shape(self, point:Iterable):
        """Returns index of shape that contain the point."""  
        sim = self.__delaunay.find_simplex(point)
        if sim > -1:
            return self.__descendent[sim]
        return -1
    
    def _candidates(self, shapes, neighbors_of_shapes):
        """
        Find the edges that could be removed. Also store additional information, 
        such as wall length and existing angles.
        """
        def prepshape(shape, wall): #rotate/flip shape so, that wall[0] is at start and wall[1] is at end.
            while len(np.intersect1d(shape[0:2], wall)) != 2:
                shape = np.roll(shape, 1)
            shape = np.roll(shape, -1) #one vwall vertice at start, the other at the end.
            if shape[0] == wall[1]:
                shape = np.flip(shape) #vwall[0] is at beginning, vwall[1] is at end
            return shape
        def vec(*vi):
            return self.__points[vi[1]] - self.__points[vi[0]]
        def angle(vecA, vecB):
            cosangle = np.dot(vecA, vecB) / (np.linalg.norm(vecA) * np.linalg.norm(vecB))
            return np.arccos(np.clip(cosangle, -1, 1))
        def PolyArea(vi):
            x, y = self.__points[vi, 0], self.__points[vi, 1]
            return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
        
        candidates = []
        for si1, neighbors in enumerate(neighbors_of_shapes):
            shape1 = shapes[si1]
            for si2 in neighbors:
                if si1 > si2: continue #only add each wall once
                shape2 = shapes[si2] 
                # Find vertices of shared wall.
                wall = np.intersect1d(shape1, shape2)
                if len(wall) != 2: continue
                # Prepare by putting wall vertice 0 (1) at position 0 (-1) in each shape
                shape1, shape2 = prepshape(shape1, wall), prepshape(shape2, wall)
                # Get candidate-polygon
                shape3 = [*shape1[:-1], *shape2[::-1][:-1]]
                h3 = ConvexHull(self.__points[shape3]).vertices
                if self.__convex and len(h3) != len(shape3): continue
                # Add characteristics.
                wallvec = vec(*wall) # pointing 0->1
                # Vectors pointing along the edges, starting at where wall is.
                vecs = [[vec(*vi) for vi in [shape1[:2], shape2[:2]]],
                        [vec(*vi) for vi in [shape1[-2:], shape2[-2:]]]]
                angles = np.array([[angle(wallvec, v) for v in vec] for vec in vecs]) #angles at corner 0, angles at corner 1
                candidates.append({
                    'wall': [*wall],
                    'si': [si1, si2], 
                    'shape3': shape3,
                    'len': np.linalg.norm(wallvec),
                    'angles_before': angles,
                    'angles_after': angles.sum(axis=1),
                    'error': np.abs(angles.sum(axis=1) - (np.pi * (1 - 2/len(shape3)))).mean()
                    })
        return candidates

    def plotpoints(self, ax, **kwargs):
        ax.plot(*self.__points.T, 'ko', **kwargs)
    def plotdelaunay(self, ax, **kwargs):
        indptr, indices = self.__delaunay.vertex_neighbor_vertices
        for vi1 in np.arange(len(self.__points)):
            for vi2 in indices[indptr[vi1]:indptr[vi1+1]]:
                if vi1 < vi2:
                    ax.plot(*self.__points[[vi1, vi2],:].T, 'k', **{'alpha':1, **kwargs})
    def plotremovablewalls(self, ax, **kwargs):
        cands = self._candidates(self.__delaunay.simplices, self.__delaunay.neighbors)
        for w in [cand['wall'] for cand in cands]:
            ax.plot(*self.__points[w, :].T, **{'color':'k', **kwargs})
    def plotpolygons(self, ax, **kwargs):
        for shape in self.__shapes:
            for vi in zip(shape, np.roll(shape, 1)):
                ax.plot(*self.__points[vi,:].T, **{'color':'b', **kwargs})
