# -*- coding: utf-8 -*-
"""
Tests and illustrations of the functions.
"""

from interpol import point2color as p2c
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import numpy as np

SR3 = np.sqrt(3)


# %% Interpolation functions in single polygon.

fig, axes = plt.subplots(2, 3, figsize=(10, 6))
fig.suptitle('interpolation with numbers (top) and colors (below)')
pix = 150
corners = 7

anchors = np.array([(np.cos(a), np.sin(a)) for a in np.pi*np.linspace(0.25, 2.25, corners, False)])
extremes = np.array([anchors.min(axis=0), anchors.max(axis=0)])
extremes = np.array([[1.3, -0.3], [-0.3, 1.3]]) @ extremes # add 30% margin
extent = extremes.T.flatten() #xmin xmax ymin ymax

for j in range(2):
    if j == 0:
        values = np.random.rand(len(anchors)) + 0.1
        f0 = p2c.triangles(anchors, values)
        f1 = p2c.polygon(anchors, values)
    else:
        values = np.random.rand(len(anchors), 3)
        f0 = p2c.triangles(anchors, values)
        f1 = p2c.polygon(anchors, values)

    data = np.array([[[v0:=f0((x,y)), v1:=f1((x,y)), np.abs(v1-v0)]
                      for x in np.linspace(extent[0], extent[1], pix)]
                     for y in np.linspace(extent[2], extent[3], pix)])
    for i, ax in enumerate(axes[j]):
        ax.imshow(data[:, :, i], origin='lower', extent=extent, cmap='gist_heat',
                  vmin=0, vmax=data.max())
        ax.set_xticks([])
        ax.set_yticks([])
        ax.plot(anchors[:,0], anchors[:,1], 'ko', markersize=4)
        ax.plot(anchors[:,0], anchors[:,1], 'wo', markersize=2)
    axes[j, 0].set_title('Standard, with triangle(s)')
    axes[j, 1].set_title('Generalised, with polygon')
    axes[j, 2].set_title('Delta')
fig.tight_layout()


# %% Interpolation functions in/around set of points.

fig, axes = plt.subplots(2, 3, figsize=(10, 6))
fig.suptitle('')
pix = 150
corners = 10

anchors = np.random.rand(corners*2).reshape(-1, 2)
extremes = np.array([anchors.min(axis=0), anchors.max(axis=0)])
extremes = np.array([[1.3, -0.3], [-0.3, 1.3]]) @ extremes # add 30% margin
extent = extremes.T.flatten() #xmin xmax ymin ymax

for j in range(2):
    if j == 0:
        values = np.random.rand(len(anchors)) + 0.1
        f0 = p2c.triangles(anchors, values)
        f1 = p2c.polygons(anchors, values)
    else:
        values = np.random.rand(len(anchors), 3)
        f0 = p2c.triangles(anchors, values)
        f1 = p2c.polygons(anchors, values)

    data = np.array([[[v0:=f0((x,y)), v1:=f1((x,y)), np.abs(v1-v0)]
                      for x in np.linspace(extent[0], extent[1], pix)]
                     for y in np.linspace(extent[2], extent[3], pix)])
    for i, ax in enumerate(axes[j]):
        ax.imshow(data[:, :, i], origin='lower', extent=extent, cmap='gist_heat',
                  vmin=0, vmax=data.max())
        ax.set_xticks([])
        ax.set_yticks([])
        ax.plot(anchors[:,0], anchors[:,1], 'ko', markersize=4)
        ax.plot(anchors[:,0], anchors[:,1], 'wo', markersize=2)
    axes[j, 0].set_title('Standard, with triangle(s)')
    axes[j, 1].set_title('Generalised, with polygon(s)')
    axes[j, 2].set_title('Delta')
fig.tight_layout()


# %% Interpolation on 1 axis.

pix = 250

# 'Normal' on 1 axis.
f = plt.get_cmap('terrain')
img = [[f(s) for s in np.linspace(0, 1, pix)]]
fig, ax = plt.subplots(1, 1)
ax.imshow(img, extent=[0, 1, -0.05, 0.05])
ax.set_yticks([])
xvals = [t[0] for t in f._segmentdata['red']]
ax.plot(xvals, [0]*len(xvals), 'ko', markersize=4)
ax.plot(xvals, [0]*len(xvals), 'wo', markersize=2)
fig.tight_layout()


# %% Interpolation on 2 axes.

pix = 250

# Diverging - income.
colorsA = [to_rgb('#40181B'), to_rgb('r')]
colorsB = [to_rgb('#0D5D00'), to_rgb('#159E00'), to_rgb('#3AD622')]
f = p2c.cmap2(colorsA, colorsB, to_rgb('#aaa'))
img = [[f(4, b) for b in np.linspace(3, 5, pix)]]
fig, ax = plt.subplots(1, 1)
ax.imshow(img, extent=[-1, 1, -.1, .1])
ax.set_yticks([])
ax.set_xticks([-1, 0, 1])
ax.set_xticklabels(['more expenses (a)', 'break even', 'more income (b)'])
xvals = np.array((*np.linspace(-1, 0, len(colorsA), False), 
                  *np.linspace(1, 0, len(colorsB), False), 0))
ax.plot(xvals, [0]*len(xvals), 'ko', markersize=4)
ax.plot(xvals, [0]*len(xvals), 'wo', markersize=2)
fig.tight_layout()

# Diverging - time.
colorsA = [to_rgb('#D87A11'), to_rgb('#F39935')]
colorsB = [to_rgb('#6C2FB1'), to_rgb('#8249C2')]
f = p2c.cmap2(colorsA, colorsB)
img = [[f(4, b) for b in np.linspace(3, 5, pix)]]
fig, ax = plt.subplots(1, 1)
ax.imshow(img, extent=[-1, 1, -.1, .1])
ax.set_yticks([])
ax.set_xticks([-1, -0.5, 0, 0.5, 1])
ax.set_xticklabels(['10:00\nbike is faster', '5:00', 'same-same', '5:00', '10:00\ncar is faster'])
xvals = np.array((*np.linspace(-1, 0, len(colorsA), False), 
                  *np.linspace(1, 0, len(colorsB), False)))
ax.plot(xvals, [0]*len(xvals), 'ko', markersize=4)
ax.plot(xvals, [0]*len(xvals), 'wo', markersize=2)
fig.tight_layout()


# %% Interpolation on 3 axes.

colorsA = [to_rgb('#D87A11'), to_rgb('#F39935')]
colorsB = [to_rgb('#6C2FB1'), to_rgb('#8249C2')]
colorsC = [to_rgb('#22751E'), to_rgb('#3C9537')]

# Partial and full image.
pix = 500
shape = (int(pix*0.8), pix)
f = p2c.cmap3(colorsA, colorsB, colorsC)
img = [np.ones((*shape, 3)), np.ones((*shape, 3))] # 2 white images
es = np.array([(0, 1), (-SR3/2, -0.5), (SR3/2, -0.5)]) #unit vectors in axes' directions
W = 0.1
ep = es @ np.array([[0, -1], [1, 0]]) #unit vectors perpendicular to axes
for j, x in enumerate(np.linspace(-1, 1, shape[1])):
    for i, y in enumerate(np.linspace(-0.65, 1.05, shape[0])):
        v = np.array((x, y))
        dot = es.dot(v) #projections onto the axes
        k = np.argmin(dot)
        if dot[k] > -0.5:
            img[1][i, j, :] = f(*dot) #full image
        k = np.argmax(dot)
        if not(0 < dot[k] < 1): continue
        if abs(ep[k].dot(v)) > W: continue
        img[0][i, j, :] = f(*np.roll([dot[k],0,0],k)) #partial image
for l in range(2):
    fig, ax = plt.subplots(1, 1)
    ax.imshow(img[l], extent=[-1, 1, -0.65, 1.05], origin='lower')
    ax.set_yticks([])
    ax.set_xticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if l == 0:
        line = np.array([[es[k] + W*ep[k], es[k] - W*ep[k], W*(-ep[k]+es[k]/SR3)] for k in range(3)]).reshape(-1, 2)
        line = np.array((*line, line[0]))
    else:
        line = np.array([es[k] for k in [0,1,2,0]])
    ax.plot(line[:,0], line[:,1], 'k')
    if l == 0:
        ax.text(*(es[0] - W*ep[0]), "bike is faster", ha='right', va='top', rotation=90)
        ax.text(*(es[1] + 1.5*W*ep[1]), "car is faster", ha='left', va='bottom', rotation=30)
        ax.text(*(es[2] - 1.5*W*ep[2]), "bus is faster", ha='right', va='bottom', rotation=-30)
    else:
        ax.text(*(es[0]+[0,0.01]), "bike is fastest", ha='center', va='bottom')
        ax.text(*(es[1]+[0,-0.03]), "car is fastest", ha='center', va='top')
        ax.text(*(es[2]+[0,-0.03]), "bus is fastest", ha='center', va='top')
    if l == 1:
        anchors = np.array([es[k] * np.linspace(1, 0, len(c), False)[:, np.newaxis] 
                            for k, c in enumerate([colorsA, colorsB, colorsC])]).reshape(-1, 2)
        ax.plot(*anchors.T, 'ko', markersize=4)
        ax.plot(*anchors.T, 'wo', markersize=2)
    fig.tight_layout()

# Special lines and points.
fig, axes = plt.subplots(2, 4, figsize=(10,7))
for ax in axes.flatten():
    ax.imshow(img[1], extent=[-1, 1, -0.65, 1.05], origin='lower')
    ax.set_yticks([])
    ax.set_xticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    line = np.array([es[k] for k in [0,1,2,0]])
    ax.plot(line[:,0], line[:,1], 'k', linewidth=0.5)
axes[0, 0].add_patch(mpl.patches.Polygon([[0,0],0.5*(es[0]+es[2]),es[0],0.5*(es[0]+es[1]),[0,0]], color='k', alpha=0.5))
axes[0, 0].set_title('a > b ∧\na > c  ')
axes[0, 1].plot(*np.array([0*es[0], 1*es[0]]).T, 'k')
axes[0, 1].set_title('a > b = c')
axes[0, 2].plot(0, 1, 'ko')
axes[0, 2].set_title('a = b + 1 ∧\na = c + 1  ')
axes[0, 3].add_patch(mpl.patches.Polygon([[0,0],es[0],0.5*(es[0]+es[1]),[0,0]], color='k', alpha=0.5))
axes[0, 3].set_title('a > b > c')
axes[1, 0].add_patch(mpl.patches.Polygon([[0,0],es[1],es[2],[0,0]], color='k', alpha=0.5))
axes[1, 0].set_title('a < b ∧\na < c  ')
axes[1, 1].plot(*np.array([0*es[0], -0.5*es[0]]).T, 'k')
axes[1, 1].set_title('a < b = c')
axes[1, 2].plot(0, -0.5, 'ko')
axes[1, 2].set_title('a = b - 0.5 ∧\na = c - 0.5  ')
axes[1, 3].add_patch(mpl.patches.Polygon([[0,0],-0.5*es[0],es[2],[0,0]], color='k', alpha=0.5))
axes[1, 3].set_title('a < b < c')
fig.tight_layout()  


# %% Some other visualisations of the color map.

pix = 150
fig, axs = plt.subplots(3, 1, figsize=(6, 9))
for i, ax in enumerate(axs.flatten()):
    if i == 0:   fun = lambda a: f(a, 0, 0)
    elif i == 1: fun = lambda b: f(0, b, 0)
    else:        fun = lambda c: f(0, 0, c)
    img = np.zeros((1, pix, 3))
    img[0, :, :] = [fun(v) for v in np.linspace(-1, 1, pix)]
    img[0, img.shape[1]//2, :] = [1., 1., 1.]
    ax.imshow(img, interpolation='none', origin='lower', extent=[-1,1,-.1,.1])
    ax.set_title('abc'[(i+1)%3] + ' and ' + 'abc'[(i+2)%3] + ' are equal')
    ax.set_xticks([-1, 0, 1])
    ax.set_xticklabels(['abc'[i] + ' slowest', 'all same speed', 'abc'[i] + ' fastest'])
    ax.set_yticks([])
    
pix = 200
fig, axs = plt.subplots(3, 1, figsize=(6, 9))
for i, ax in enumerate(axs.flatten()):
    if i == 0:   fun = lambda a, bc: f(a, max(bc, 0), -min(bc, 0))
    elif i == 1: fun = lambda b, ac: f(-min(ac, 0), b, max(ac, 0))
    else:        fun = lambda c, ab: f(max(ab, 0), -min(ab, 0), c)
    img = np.ones((pix, pix*4, 3))
    img[:, :, :] = [[fun(x, y*x/(0.3+0.7*x**2)) if abs(y) < 0.3+0.7*x**2 else (1., 1., 1.)
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