# -*- coding: utf-8 -*-
"""
Tests and illustrations of the functions.
"""

import interpol.point2color as p2c 
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.colors import to_rgba
import numpy as np

SR3 = np.sqrt(3)


# %% Illustration of interpolation function in/around single polygon.

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
fig.suptitle('Interpolation in/around polygon, with `polygon` function.')
pix = 150
count = 7

anchors = np.array([(np.cos(a), np.sin(a)) for a in 
                    np.pi*np.linspace(0.25, 2.25, count, False)]) + \
                    0.5*np.random.rand(count, 2)
extremes = np.array([anchors.min(axis=0), anchors.max(axis=0)])
extremes = np.array([[1.3, -0.3], [-0.3, 1.3]]) @ extremes # add 30% margin
extent = extremes.T.flatten() #xmin xmax ymin ymax

for j, ax in enumerate(axes):
    if j == 0:
        values = np.random.rand(count) + 0.1
        ax.set_title('Anchor points define floats')
    else:
        values = np.random.rand(count, 3)
        ax.set_title('Anchor points define colors')
    
    f = p2c.polygon(anchors, values)
    data = np.array([[f((x, y))
                      for x in np.linspace(extent[0], extent[1], pix)]
                     for y in np.linspace(extent[2], extent[3], pix)])
    ax.imshow(data, origin='lower', extent=extent, cmap='gist_heat',
              vmin=0, vmax=data.max())
    ax.set_xticks([])
    ax.set_yticks([])
    ax.plot(*anchors.T, 'ko', markersize=4)
    ax.plot(*anchors.T, 'wo', markersize=2)
    ax.plot(*np.array([*anchors, anchors[0]]).T, 'gray', linewidth=0.5)

fig.tight_layout()


# %% Illustration of interpolation functions in/around set of points.

fig, axes = plt.subplots(3, 3, figsize=(10, 10))
fig.suptitle('Interpolation around set of points with `triangles` and `polygons` functions.')
pix = 100
count = 10

anchors = np.random.rand(count, 2)
extremes = np.array([anchors.min(axis=0), anchors.max(axis=0)])
extremes = np.array([[1.3, -0.3], [-0.3, 1.3]]) @ extremes # add 30% margin
extent = extremes.T.flatten() #xmin xmax ymin ymax

for j in range(2):
    if j == 0:
        values = np.random.rand(count) + 0.1
    else:
        values = np.random.rand(count, 3)
        
    f0 = p2c.triangles(anchors, values)
    f1 = p2c.polygons(anchors, values)
    data = np.array([[[v0:=f0((x,y)), v1:=f1((x,y)), np.abs(v1-v0)]
                      for x in np.linspace(extent[0], extent[1], pix)]
                     for y in np.linspace(extent[2], extent[3], pix)])

    for i in range(3):
        ax = axes[i, j+1]
        ax.imshow(data[:, :, i], origin='lower', extent=extent, cmap='gist_heat',
                  vmin=0, vmax=data.max())  
        
for i in range(2):
    ax = axes[i, 0]
    if i == 0:
        shapes = f0.delaunay.simplices
    else:
        shapes = f1.polygonate.shapes
    for shape in shapes:
        for vi in zip(shape, np.roll(shape, 1)):
            ax.plot(*anchors[vi,:].T, 'gray', linewidth=0.5)

for j in range(3):
    for i in range(3):
        ax = axes[i, j]
        ax.set_xlim(extent[:2])
        ax.set_ylim(extent[2:])
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.plot(*anchors.T, 'ko', markersize=4)
        ax.plot(*anchors.T, 'wo', markersize=2)
        if i == 0: 
            title = ['Tesselation...', 'Interpolation with floats...', 'Interpolatation with colors...'][j] + '\n'
        else:
            title = ''
        title += ['...using `triangles`', '...using `polygons`', 'delta'][i]
        ax.set_title(title)
axes[2, 0].set_visible(False)
            
fig.tight_layout()


# %% Illustration for interpolation on 1 axis.

pix = 250

# 'Normal' on 1 axis.
cmap = plt.get_cmap('terrain')
img = [[cmap(s) for s in np.linspace(0, 1, pix)]]
fig, ax = plt.subplots(1, 1)
ax.imshow(img, extent=[0, 1, -0.05, 0.05])
ax.set_yticks([])
xvals = [t[0] for t in cmap._segmentdata['red']]
ax.plot(xvals, [0]*len(xvals), 'ko', markersize=4)
ax.plot(xvals, [0]*len(xvals), 'wo', markersize=2)
fig.tight_layout()


# %% Illustrations for interpolation on 2 axes.

# Income vs expenses.

colorsA = ['#40181B', 'red']
colorsB = ['#0D5D00', '#159E00', '#3AD622']
cmap = p2c.ColorMap2(colorsA, colorsB, '#aaa')
fig, ax = plt.subplots(1, 1)
cmap.colorbar(ax)
ax.set_ticks([(1,0), (0,0), (0,1)])
ax.set_ticklabels(['more expenses (a)', 'break even', 'more income (b)'])

#additionally, show anchors
xvals = np.array((*[cmap._uv_to_x((u, 0)) for u in np.linspace(1, 0, len(colorsA), False)], 
                  *[cmap._uv_to_x((0, v)) for v in np.linspace(1, 0, len(colorsB), False)], 
                  cmap._uv_to_x((0, 0))))
xy = np.array([xvals, np.zeros(len(xvals))]).T
ax.plot(xy[:, 0], xy[:, 1], 'ko', markersize=4)
ax.plot(xy[:, 0], xy[:, 1], 'wo', markersize=2)
fig.tight_layout()


# Bike vs car.

colorsA = ['#D87A11', '#F39935']
colorsB = ['#6C2FB1', '#8249C2']
cmap = p2c.ColorMap2(colorsA, colorsB)
fig, ax = plt.subplots(1, 1, figsize=(6,6))
cmap.colorbar(ax, orientation='vertical')
ax.set_ticks([(1, 0), (0.5, 0), (0, 0), (0, 0.5), (0, 1)])
ax.set_ticklabels(['10:00\nbike is fastest', '5:00', 'same-same', '5:00', '10:00\ncar is fastest'])

#additionally, show anchors
xvals = np.array((*[cmap._uv_to_x((u, 0)) for u in np.linspace(1, 0, len(colorsA), False)], 
                  *[cmap._uv_to_x((0, v)) for v in np.linspace(1, 0, len(colorsB), False)]))
xy = np.array([xvals, np.zeros(len(xvals))]).T
ax.plot(xy[:, 1], xy[:, 0], 'ko', markersize=4)
ax.plot(xy[:, 1], xy[:, 0], 'wo', markersize=2)
fig.tight_layout()


# %% Illustrations for interpolation on 3 axes.

colorsA = ['#D87A11', '#F39935']
colorsB = ['#6C2FB1', '#8249C2']
colorsC = ['#22751E', '#3C9537']
cmap = p2c.ColorMap3(colorsA, colorsB, colorsC)
pix = 300


# Color triangle legend.
fig, ax = plt.subplots()
cmap.colortriangle(ax, pix)
#additionally, show anchors
es = np.array([(0, 1), (-SR3/2, -0.5), (SR3/2, -0.5)]) #unit vectors in axes' directions
ax.text(*(es[0]+[0,0.01]), "bike is fastest", ha='center', va='bottom')
ax.text(*(es[1]+[0,-0.03]), "car is fastest", ha='center', va='top')
ax.text(*(es[2]+[0,-0.03]), "bus is fastest", ha='center', va='top')
anchors = np.array([es[k] * np.linspace(1, 0, len(c), False)[:, np.newaxis] 
                    for k, c in enumerate([colorsA, colorsB, colorsC])]).reshape(-1, 2)
ax.plot(*anchors.T, 'ko', markersize=4)
ax.plot(*anchors.T, 'wo', markersize=2)


# Color bars legend.
fig, ax = plt.subplots()
cmap.colorbars(ax, pix)


# Illustration of color bars in triangle shape.
shape = (int(pix*0.8), pix)
img = np.ones((*shape, 4)) #white image
es = np.array([(0, 1), (-SR3/2, -0.5), (SR3/2, -0.5)]) #unit vectors in axes' directions
W = 0.1
ep = es @ np.array([[0, -1], [1, 0]]) #unit vectors perpendicular to axes
for j, x in enumerate(np.linspace(-1, 1, shape[1])):
    for i, y in enumerate(np.linspace(-0.65, 1.05, shape[0])):
        v = np.array((x, y))
        dot = es.dot(v) #projections onto the axes
        k = np.argmax(dot)
        if not(0 < dot[k] < 1): continue
        if abs(ep[k].dot(v)) > W: continue
        img[i, j, :] = cmap.color(*np.roll([dot[k],0,0],k)) #partial image
fig, ax = plt.subplots(1, 1)
ax.imshow(img, extent=[-1, 1, -0.65, 1.05], origin='lower')
plt.axis(False)
line = np.array([[es[k] + W*ep[k], es[k] - W*ep[k], W*(-ep[k]+es[k]/SR3)] 
                 for k in range(3)]).reshape(-1, 2)
line = np.array((*line, line[0]))
ax.plot(line[:,0], line[:,1], 'k')
ax.text(*(es[0] - W*ep[0]), "bike is faster", ha='right', va='top', rotation=90)
ax.text(*(es[1] + 1.5*W*ep[1]), "car is faster", ha='left', va='bottom', rotation=30)
ax.text(*(es[2] - 1.5*W*ep[2]), "bus is faster", ha='right', va='bottom', rotation=-30)
fig.tight_layout()


# Illustration of special lines and points.
fig, axes = plt.subplots(2, 4, figsize=(10,7))
for i, ax in enumerate(axes.flatten()):
    if i == 0:
        cmap.colortriangle(ax, pix) #calculate once
        data = ax.get_images()[0].get_array()
    ax.imshow(data, origin='lower', extent=[-1, 1, -0.65, 1.05])
    ax.set_yticks([])
    ax.set_xticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    line = np.array([es[k] for k in [0,1,2,0]])
    ax.plot(line[:,0], line[:,1], 'k', linewidth=0.5)
axes[0, 0].add_patch(Polygon([[0,0],0.5*(es[0]+es[2]),es[0],0.5*(es[0]+es[1]),[0,0]], color='k', alpha=0.5))
axes[0, 0].set_title('a > b ∧\na > c  ')
axes[0, 1].plot(*np.array([0*es[0], 1*es[0]]).T, 'k')
axes[0, 1].set_title('a > b = c')
axes[0, 2].plot(0, 1, 'ko')
axes[0, 2].set_title('a = b + 1 ∧\na = c + 1  ')
axes[0, 3].add_patch(Polygon([[0,0],es[0],0.5*(es[0]+es[1]),[0,0]], color='k', alpha=0.5))
axes[0, 3].set_title('a > b > c')
axes[1, 0].add_patch(Polygon([[0,0],es[1],es[2],[0,0]], color='k', alpha=0.5))
axes[1, 0].set_title('a < b ∧\na < c  ')
axes[1, 1].plot(*np.array([0*es[0], -0.5*es[0]]).T, 'k')
axes[1, 1].set_title('a < b = c')
axes[1, 2].plot(0, -0.5, 'ko')
axes[1, 2].set_title('a = b - 0.5 ∧\na = c - 0.5  ')
axes[1, 3].add_patch(Polygon([[0,0],-0.5*es[0],es[2],[0,0]], color='k', alpha=0.5))
axes[1, 3].set_title('a < b < c')
fig.tight_layout()  


# %% Some other visualisations of the color map.

pix = 150
fig, axs = plt.subplots(3, 1, figsize=(6, 9))
for i, ax in enumerate(axs.flatten()):
    if i == 0:   fun = lambda a: f.color(a, 0, 0)
    elif i == 1: fun = lambda b: f.color(0, b, 0)
    else:        fun = lambda c: f.color(0, 0, c)
    img = np.zeros((1, pix, 4))
    img[0, :, :] = [fun(v) for v in np.linspace(-1, 1, pix)]
    img[0, img.shape[1]//2, :] = to_rgba('w')
    ax.imshow(img, interpolation='none', origin='lower', extent=[-1,1,-.1,.1])
    ax.set_title('abc'[(i+1)%3] + ' and ' + 'abc'[(i+2)%3] + ' are equal')
    ax.set_xticks([-1, 0, 1])
    ax.set_xticklabels(['abc'[i] + ' slowest', 'all same speed', 'abc'[i] + ' fastest'])
    ax.set_yticks([])
    
pix = 200
fig, axs = plt.subplots(3, 1, figsize=(6, 9))
for i, ax in enumerate(axs.flatten()):
    if i == 0:   fun = lambda a, bc: f.color(a, max(bc, 0), -min(bc, 0))
    elif i == 1: fun = lambda b, ac: f.color(-min(ac, 0), b, max(ac, 0))
    else:        fun = lambda c, ab: f.color(max(ab, 0), -min(ab, 0), c)
    img = np.ones((pix, pix*4, 4))
    img[:, :, :] = [[fun(x, y*x/(0.3+0.7*x**2)) if abs(y) < 0.3+0.7*x**2 else to_rgba('w')
                     for x in np.linspace(0, 1, pix*4)]
                    for y in np.linspace(-1, 1, pix)]
    ax.imshow(img, interpolation='none', origin='lower', extent=[0,10,-1,1])
    ax.set_title('abc'[i] + ' fastest')
    ax.yaxis.tick_right()
    ax.set_xticks([0])
    ax.set_xticklabels(['all same speed'])
    ax.set_yticks([-1, 1])
    ax.set_yticklabels(['abc'[(i+2)%3] + ' second fastest',
                        'abc'[(i+1)%3] + ' second fastest'])
    for s in ax.spines.values():
        s.set_visible(False)