# -*- coding: utf-8 -*-
"""
Tests and illustrations of the functions.
"""

import interpol.point2color as p2c 
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import numpy as np
from numpy import random as rnd



# %% Sample use of interpolation function in/around single polygon.

anchors = [[0,0], [1,0], [1,1], [0,1]] #in order of appearance in polygon

# Interpolation with numeric values.
values = [1, 5, 0, 2]
f = p2c.polygon(anchors, values)
f((0.5, 0.6)) #1.8

# Interpolation with colors.
values = [to_rgba('w'), to_rgba('r'), to_rgba('g'), to_rgba('c')]
f = p2c.polygon(anchors, values)
f((0.5, 0.6)) #array([0.4  , 0.575, 0.425, 1.   ])



# %% Sample use of interpolation functions in/around set of points.

anchors = np.random.rand(10, 2) #random points in plane

# Interpolation with numeric values.
values = np.random.rand(10) #random value for each point
# using triangles
f = p2c.triangles(anchors, values)
f((0.5, 0.6)) #0.69...
# using polygons
f = p2c.polygons(anchors, values)
f((0.5, 0.6)) #0.71...

# Interpolation with colors: analogously


# %% Sample use of ColorMap2, short

cmap = p2c.ColorMap2(['#D87A11', '#F39935'], ['#6C2FB1'])
cmap.color(0.2, 0.5) #(0.6699, 0.3778, 0.4678, 1.0) (includes alpha)


# %% Sample use of ColorMap2, longer

# Set-up canvas
fig = plt.figure(figsize=(7.5, 9))
ax = fig.add_axes([0, 0.2, 1, 0.8])
cax = fig.add_axes([0.6, 0, 0.4, 0.3])
fig.suptitle('color of circle based on x-coordinate and y-coordinate')

# Set-up colormap
cmap = p2c.ColorMap2(['#7A11D8', '#9935F3'], ['#2FB16C'])

# Add sample data
for i in range(100):
    a, b = rnd.rand(2)
    col = cmap.color(a, b) #this is where the magic happens
    circle = plt.Circle((a, b), 0.05, color=col)
    ax.add_artist(circle)
    
# Add legends
cmap.colorbar(cax)
cax.set_ticks([(1,0), (0,0), (0,1)], '')
cax.set_ticklabels(['x', 'same', 'y'])


# %% Sample use of ColorMap3, short

cmap = p2c.ColorMap3(['#D87A11', '#F39935'], ['#6C2FB1'], ['blue'])
cmap.color(0.2, 0.5, -0.3) #(0.5864, 0.2974, 0.4528, 1.0)


# %% Sample use of ColorMap3, longer

# Set-up canvas
fig = plt.figure(figsize=(7.5, 9))
ax = fig.add_axes([0, 0.2, 1, 0.8])
cax1 = fig.add_axes([0.6, 0, 0.4, 0.3])
cax2 = fig.add_axes([0, 0, 0.2, 0.12])
fig.suptitle('color of circle based on x-coordinate, y-coordinate, and radius')

# Set-up colormap
cmap = p2c.ColorMap3(['#7A11D8', '#9935F3'], 
                     ['#2FB16C'], 
                     ['#751E22', '#95373C'],
                     'gray',
                     maxdiff=1.2)

# Add sample data
for i in range(100):
    a, b, c = rnd.rand(3)
    circle = plt.Circle((a, b), c/12, color=cmap.color(a, b, c))
    ax.add_artist(circle)
    
# Add legends
cmap.colortriangle(cax1)
cax1.set_ticks([(1.2,0,0), (0,1.2,0), (0,0,1.2)])
cax1.set_ticklabels(['x', 'y', 'r'], bbox={'alpha':0})
cmap.colorbars(cax2, orientation='vertical')
for i, ax in enumerate(cax2.subaxes):
    ax.set_title('xyr'[i])
    ax.set_ticks([])
    
