"""
Color interpolation functions.

2020-09 rwijtvliet@gmail.com
"""    
    
from matplotlib.colors import to_rgb, LinearSegmentedColormap
import numpy as np
from typing import Iterable, Callable
from scipy.spatial import Delaunay, ConvexHull
import interpol.polygonate as pg

SR3 = np.sqrt(3)


#%% Interpolation function in single polygon.

def polygon(anchorpoints:Iterable, anchorvalues:Iterable) -> Callable:
    """
    Generalised barycentric interpolation function for point in a single 
    polygon.
    
    Polygon has n vertices ("anchorpoints" of shape (n, 2)) in order
    and n values ("anchorvalues", e.g. iterable of n floats, or np-array of 
    shape (n, 3) with colors).
    
    The returned function takes a point (2-Iterable) and returns the inter-
    polated value (or extrapolated value if the point lies outside the shape).
    
    Source: https://cgvr.cs.uni-bremen.de/papers/barycentric/barycentric.pdf
    """
    anchorpoints = np.array(anchorpoints)
    anchorvalues = np.array(anchorvalues)
    if len(anchorpoints) != len(anchorvalues):
        raise ValueError("Parameters 'anchorpoints' and 'anchorvalues' must be of equal length.")
    if len(anchorpoints) < 3:
        raise ValueError("At least 3 anchorpoints must be specified.")
        
    def interp(point):
        dv = anchorpoints - point #distance vector to each anchorpoint
        r = np.linalg.norm(dv, axis=1) #distances to each anchorpoint
        for i, rr in enumerate(r): 
            if rr == 0: 
                return anchorvalues[i] #return anchor value if point is on anchor point
        dv_next = np.roll(dv, -1, axis=0) # neighbor to each dv
        A = np.array([np.linalg.det([d1, d2])/2 for d1, d2 in zip(dv, dv_next)]) #determinant of each displ. vector with neighbor
        D = np.array([np.dot(d1, d2) for d1, d2 in zip(dv, dv_next)])
        for i, (aa, dd) in enumerate(zip(A, D)):
            if aa == 0 and dd < 0:
                j = (i + 1) % len(anchorpoints)
                return (r[j] * anchorvalues[i] + r[i] * anchorvalues[j]) / (r[i] + r[j])
        W = 0
        value = np.zeros(anchorvalues[0].shape)
        for i, anchorvalue in enumerate(anchorvalues):
            i_prev = (i-1)%len(anchorvalues)
            i_next = (i+1)%len(anchorvalues)
            w = 0
            if A[i_prev] != 0:
                w += (r[i_prev] - D[i_prev]/r[i])/A[i_prev]
            if A[i] != 0:
                w += (r[i_next] - D[i]/r[i])/A[i]
            value += anchorvalue * w
            W += w
        if W == 0: #no weights - should be impossible
            return ValueError(f'Point {point} is outside the shape.')
        return value/W
    
    return interp


#%% Interpolation functions in plane with points.

def triangles(anchorpoints:Iterable, anchorvalues:Iterable, valoutside=None) -> Callable:
    """
    Standard barycentric interpolation function for point inside a single 
    triangle, or, more generally, for point inside the convex hull around a 
    set of anchor points that is tesselated with triangles.
    
    Triangle has 3 vertices ("anchorpoints" of shape (3, 2)) and 3 values
    ("anchorvalues", e.g. Iterable of 3 floats, or of 3 colors). When passed 
    >3 anchorpoints, the area is tesselated with triangles.
    
    The returned function takes a point (2-Iterable) and returns the inter-
    polated value. If the point lies outside the shape, "valoutside" is 
    returned (if specified) or an attempt at extrapolation is done (default).
    
    Source: https://stackoverflow.com/questions/57863618/how-to-vectorize-calculation-of-barycentric-coordinates-in-python
    """
    anchorpoints = np.array(anchorpoints)
    anchorvalues = np.array(anchorvalues)
    if len(anchorpoints) != len(anchorvalues):
        raise ValueError("Parameters 'anchorpoints' and 'anchorvalues' must be of equal length.")
    if len(anchorpoints) < 3:
        raise ValueError("At least 3 anchorpoints must be specified.")
    
    # Tesselate into simplexes (individual triangles).
    delaunay = Delaunay(anchorpoints) # each row has indices of the 3 anchorpoints that are the simplex corners.
    
    def interp(point):
        # Find simplex point is in.
        s = delaunay.find_simplex(point) #simplex-index that contains point. (-1 means point is in none)
        if s > -1: #normal point, inside the hull
            # Get barycentric coordinates of the triangle.
            b0 = delaunay.transform[s, :2].dot((point - delaunay.transform[s, 2]))
            weights = np.array([*b0, 1 - b0.sum()]) # add final coordinate / weight.
            indices = delaunay.simplices[s]
        else: #point outside the hull
            if valoutside:
                return valoutside
            # Find the 2 anchorpoints on the hull line nearest to the point
            hull = ConvexHull([*anchorpoints, point], qhull_options='QG' + str(len(anchorpoints))) 
            visible = hull.simplices[hull.good] #lines (anchorpoints) visible from the point
            for indices in visible: #anchor-indices of visible line
                p01 = point - anchorpoints[indices] #from line anchors to point
                lin = anchorpoints[indices[0]] - anchorpoints[indices[1]]
                dot12 = p01.dot(lin)
                if np.sign(dot12).sum() == 0: #inside line 'shadow' if one dot product <0, >0
                    lens = np.linalg.norm(p01, axis=1)
                    lens = np.abs(dot12)
                    weights = np.flip(lens) / lens.sum()
                    break
            else: #not in shadow of line - use value of nearest anchor.
                #Find nearest anchor (="anchor 0"). Must be included in 2 lines.
                indices = list(set(visible.flatten()))
                sd = ((anchorpoints[indices] - point)**2).sum(axis=1) #squared distance to each anchorpoint
                indices = [indices[np.argmin(sd)]] #keep only nearest one
                weights = [1]               

        # Get interpolated value.
        value = np.dot(anchorvalues[indices].T, weights)
        return value 
    return interp


def polygons(anchorpoints:Iterable, anchorvalues:Iterable) -> Callable:
    """
    Generalised barycentric interpolation function for point inside a single 
    polygon, or, more generally, for point inside the convex hull around a 
    set of anchor points that is tesselated with polygons.
    
    The function tessalates when necessary.
    
    The returned function takes a point (2-Iterable) and returns the inter-
    polated value (or extrapolated value if the point lies outside the shape).
    """
    anchorpoints = np.array(anchorpoints)
    anchorvalues = np.array(anchorvalues)
    if len(anchorpoints) != len(anchorvalues):
        raise ValueError("Parameters 'anchorpoints' and 'anchorvalues' must be of equal length.")
    if len(anchorpoints) < 3:
        raise ValueError("At least 3 anchorpoints must be specified.")
    
    # Tesselate into polygons.
    polygonate = pg.Polygonate(anchorpoints, convex=False)
    # Interpolation function for each polygon...
    interpf = [polygon(anchorpoints[shape], anchorvalues[shape])
               for shape in polygonate.shapes]
    # ...and inter(extra)polation function for the hull.
    hull = ConvexHull(anchorpoints).vertices
    interpf_hull = polygon(anchorpoints[hull], anchorvalues[hull])
    
    def interp(point):
        # Find simplex point is in.
        s = polygonate.find_shape(point) #simplex that contains point. (-1 means point is in none)
        if s > -1: #normal point, inside the hull
            return interpf[s](point)
        else: #point outside the hull
            return interpf_hull(point)
    return interp


#%% Interpolation on 2 axes.

def cmap2(colorsA:Iterable=[to_rgb('r')], colorsB:Iterable=[to_rgb('g')],
               colorCenter=None, maxdiff=1, ratio:bool=False) -> Callable:
    """
    Create function to find color from 2 float values.
    
    Arguments:
        colorsA is a list of colors, which are used if the value on axis A is 
        larger than than on axis B. colorCenter, if specified, is the color 
        used if values on both axes are equal.
        maxdiff: if the difference (if ratio == False) or the ratio between
        values is larger than this maximum, the colors are clipped.
        ratio: specifies if difference (False) or ratio (True, default) between 
        the values must be calculated. 

    Returns:
        Function that accepts values a and b, and returns the color that 
        best fits this datapoint. Smallest value is used as reference for the
        other.
    """
    lst = []
    for c, v in zip(colorsA, np.linspace(0, 0.5, len(colorsA), False)):
        lst.append((v, c))
    if colorCenter:
        lst.append((0.5, colorCenter))
    for c, v in zip(np.flip(colorsB, axis=0),
                    np.flip(np.linspace(1, 0.5, len(colorsB), False))):
        lst.append((v, c))
    cmap = LinearSegmentedColormap.from_list('2axes', lst)
    f = lambda a, b: cmap(0.5 + 0.5*(b - a))

    def color(*ab:float):
        ab = np.array(ab) / maxdiff #mapping full scale to [0, 1]
        return f(*ab)
        
        # # Turn into correct range: one component being 0, and others maximum of 1.
        # if ratio:
        #     minval = np.min(abc)
        #     if minval <= 0: minval = 0.01
        #     abc = mapper(abc / minval)
        # else:
        #     abc = mapper(abc - np.min(abc))
        # # 3D to 2D, inside the triangle
        # times = np.argmax(abc) #index of best value; determines how often to rotate to project onto A or axis a
        # point = np.linalg.matrix_power(R2, times) \
        #     @ T @ np.linalg.matrix_power(R3, times) @ abc
        # return f(point)

    return color


#%% Interpolation on 3 axes.

def interp_on_3axes(valuesA:Iterable, valuesB:Iterable, valuesC:Iterable,
                    valueCenter=None) -> Callable:
    """
    Create function for interpolation on three axes, with specified values on 
    the axes. The axes meet at the center, where a value can also be specified.
    NB: the word 'value' can be substituted with 'color' if used for color 
    interpolation.
    
    Arguments:
        valuesA (and valuesB and valuesC) is an Iterable of values. The first
        value is put at the extreme end of axis A. Other values (if any) are
        linearly distributed between the end and the center.
        valueCenter is a single value. If specified, it is put at the center.
        
    Returns:
        Function that accepts values (a, b, c) on each axis, and returns the
        interpolated value. (a, b, c) should be barycentric coordinates, but 
        are forced in correct domain if they are not.
    
    Implementation details:
        valueCenter specified:
            Interpolation in 3 triangular sections, namely (A, B, center), 
            (B, C, center), (C, A, center). 
        valueCenter not specified:
            Interpolation in central triangle with final elements in valuesA, 
            valuesB, valuesC. If valuesA (or valuesB or valuesC) contains >1
            element, there exist additional shapes between the corners ABC and 
            the central triangle, and interpolation is done for these as well.
        If a shape is triangular and has its value specified only at its corners,
        'normal' barycentric coordinates are used. Otherwise, 'generalised' bary-
        centric coordinates are used.
    """
    from scipy.spatial import Delaunay
    valuesABC = np.array([np.array(valuesA), np.array(valuesB), np.array(valuesC)])
    cornersABC = np.array([(0, 1), (-SR3/2, -0.5), (SR3/2, -0.5)])

    anchors_ax = [[{'point': corner * (1-i/len(values)), 'value': value} 
                   for i, value in enumerate(values)]
                  for corner, values in zip(cornersABC, valuesABC)]
    anchor_ctr = {'point': (0,0), 'value': valueCenter} if valueCenter else None
    
    def polygon(*anchors):
        points = [anchor['point'] for anchor in anchors]
        values = [anchor['value'] for anchor in anchors]
        return {'shape': Delaunay(points), 'f': polygon(points, values)}

    if anchor_ctr is not None:
        polys = []
        for i in range(3):
            a1, a2 = anchors_ax[(i+1)%3], anchors_ax[i]
            polys.append(polygon(*a1, anchor_ctr, *np.flip(a2)))
    else:
        polys = [polygon(*[anchors[-1] for anchors in anchors_ax])] #central triangle
        for i in range(3):
            a1, a2 = anchors_ax[(i+1)%3], anchors_ax[i]
            if len(a1) > 1 or len(a2) > 1:
                polys.append(polygon(*a1, *np.flip(a2))) #outside shapes (if any)
        
    def interp(*abc:float):
        abc = np.array(abc)
        abc = abc - abc.min() #does not change location of point
        if abc.sum() > 1:
            abc = abc / abc.sum() #make sure it stays inside triangle
        point = abc.dot(cornersABC)
        for p in polys:
            if p['shape'].find_simplex(point) > -1:
                return p['f'](point)
        return None #should never happen
    
    return interp


def cmap3(colorsA:Iterable=[to_rgb('r')], colorsB:Iterable=[to_rgb('g')], 
               colorsC:Iterable=[to_rgb('b')], colorCenter=None, 
               maxdiff=1, ratio:bool=False) -> Callable:
    """
    Create function to find color from 3 float values.
    
    Arguments:
        colorsA is a list of colors, which are used if the value on axis a is 
        largest (and those on b and c are equal). In other cases, there is an 
        interpolation between the colors of the axes with the largest values.
        colorCenter, if specified, is the color used if values on all axes are 
        equal.
        maxdiff: if the difference (if ratio == False) or the ratio between
        values is larger than this maximum, the colors are clipped.
        ratio: specifies if difference (False) or ratio (True, default) between 
        the values must be calculated. 

    Returns:
        Function that accepts values a, b, and c, and returns the color that 
        best fits this datapoint. Smallest value is used as reference for the
        other ones.
    """
    f = interp_on_3axes(colorsA, colorsB, colorsC, colorCenter)
    mindiff = 1 if ratio else 0
    # mapper = interp1d((mindiff, maxdiff), (0, 1), bounds_error=False, fill_value=(0, 1))

    # R2 = np.array([[-1, -SR3/2],
    #                 [SR3/2,  -1]]) #rotate 120 degrees in 2D (i.e. A->B->C)
    # R3 = np.array([[0, 1, 0],
    #                 [0, 0, 1],
    #                 [1, 0, 0]]) #rotate through axes in 3D a<-b<-c
    # T = np.array([[0, -SR3/4, SR3/4],
    #               [1,   -3/4,  -3/4]]) # (a, b, c) to (x, y)
    
    def color(*abc:float):
        abc = np.array(abc) / maxdiff #mapping full scale to [0, 1]
        return f(*abc)
        
        # # Turn into correct range: one component being 0, and others maximum of 1.
        # if ratio:
        #     minval = np.min(abc)
        #     if minval <= 0: minval = 0.01
        #     abc = mapper(abc / minval)
        # else:
        #     abc = mapper(abc - np.min(abc))
        # # 3D to 2D, inside the triangle
        # times = np.argmax(abc) #index of best value; determines how often to rotate to project onto A or axis a
        # point = np.linalg.matrix_power(R2, times) \
        #     @ T @ np.linalg.matrix_power(R3, times) @ abc
        # return f(point)

    return color



    