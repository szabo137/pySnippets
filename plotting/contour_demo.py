#!/usr/bin/env python
"""
Illustrate simple contour plotting, contours on an image with
a colorbar for the contours, and labelled contours.

See also contour_image.py.
"""
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from reader import *

#matplotlib.rcParams['xtick.direction'] = 'out'
#matplotlib.rcParams['ytick.direction'] = 'out'

tempData = dataReader('res.dat')
(x1,y1,z1)=(tempData[1],tempData[2],tempData[3])
X,Y,Z = buildMesh2(x1,y1,z1,15,21)
"""
delta = 0.025
x = np.arange(-3.0, 3.0, delta)
y = np.arange(-2.0, 2.0, delta)
X, Y = np.meshgrid(x, y)
Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# difference of Gaussians
Z = 10.0 * (Z2 - Z1)
print Z
"""
x=np.linspace(0,1)
y = np.linspace(0,1)
X,Y = np.meshgrid(x,y)
Z=X + Y
# Or you can use a colormap to specify the colors; the default
# colormap will be used for the contour lines
plt.figure()
print Z
im = plt.imshow(Z,cmap=cm.Reds,extent=(0.0,1.0,1.0,0.0))
#levels = np.arange(-1.2, 1.6, 0.2)
#CS = plt.contour(Z, levels,
#                origin='lower',
#                 linewidths=2,
#                 extent=(-3, 3, -2, 2))

# Thicken the zero contour.
#zc = CS.collections[6]
#plt.setp(zc, linewidth=4)

#plt.clabel(CS, levels[1::2],  # label every second level
#           inline=1,
#           fmt='%1.1f',
#           fontsize=14)

# make a colorbar for the contour lines
#CB = plt.colorbar(CS, shrink=0.8, extend='both')

plt.title('Lines with colorbar')
#plt.hot()  # Now change the colormap for the contour lines and colorbar
#plt.flag()

# We can still add a colorbar for the image, too.
#CBI = plt.colorbar(im, orientation='horizontal', shrink=0.8)

# This makes the original colorbar look a bit out of place,
# so let's improve its position.

#l, b, w, h = plt.gca().get_position().bounds
#ll, bb, ww, hh = CB.ax.get_position().bounds
#CB.ax.set_position([ll, b + 0.1*h, ww, h*0.8])


plt.show()
