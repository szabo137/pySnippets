"""
module to test the contour plotting
"""
import matplotlib
#matplotlib.use('Agg') #maybe in config
import matplotlib.pyplot as plt
import numpy as np
from reader import *
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
#matplotlib.rcParams['xtick.direction'] = 'out'
#matplotlib.rcParams['ytick.direction'] = 'out'
from scipy import interpolate


fig = plt.figure()
#ax = fig.add_subplot(111)
#get data
tempData = dataReader('res.dat')
(x1,y1,z1)=(tempData[1],tempData[2],tempData[3])
xTest = x1.reshape(15,21)[0]
yTest = np.array([el[0] for el in y1.reshape(15,21)])
#print xTest
#print yTest
#print np.meshgrid(xTest,yTest)
#build mesh
XX,YY,ZZ = buildMesh2(x1,y1,z1,15,21)
#f = interpolate.interp2d(XX, YY, ZZ, kind='cubic')
#print XX
#print YY
#print ZZ
#build contour
#im = plt.imshow(f(np.linspace(0.01,0.8),np.linspace(0.75,1.87)),cmap=cm.Reds,extent=(x1.min(),x1.max(),y1.max(),y1.min()))
#plt.ylim(0.0,0.9)
#plt.xlim(0.7,1.87)
im = plt.imshow(ZZ.T, interpolation='bilinear', origin='lower',cmap=cm.Reds,extent=(x1.min(),x1.max(),y1.max(),y1.min()))
#cset = plt.contourf(YY, XX, ZZ, interpolation='bilinear',cmap=cm.Reds)
#cset = ax.contourf(XX, YY, ZZ, zdir='z', offset=-0.04, cmap=cm.Reds)
xPlot = np.linspace(0.01,0.8)
yPlot = np.linspace(0.75,1.87)
XXp,YYp = np.meshgrid(xPlot,yPlot)
CS = plt.contour(XX,YY,ZZ,4,colors='k', origin='lower')
#CS = plt.contour(XXp,YYp,f(xPlot,yPlot),4,colors='r', origin='lower')
#plt.clabel(CS, inline=1, fontsize=6)
#ax.plot_surface(XX,YY,ZZ,cmap=cm.Reds)
#ax.scatter(x1, y1, z1, c='b', marker='.')
plt.show()
