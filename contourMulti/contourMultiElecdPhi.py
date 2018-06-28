import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, rc, ticker, colors
import scipy.ndimage as spn
from plotter.settings import *
from plotter.reader import *

def getData(job,obsName):
    tempData = dataReader(getMydataPath('res.dat',job,obsName))
    rows=[4,5,6]
    (x1,y1,z1)=(tempData[rows[0]],tempData[rows[1]],tempData[rows[2]])
    ppr1,ppr2 = [128,128]
    Xraw,Yraw,Zraw = buildMesh2(x1,y1,z1,ppr1,ppr2)
    XX,YY,ZZ=Xraw,Yraw,Zraw
    conData = spn.zoom(ZZ,4)
    conData2 = spn.zoom(XX,4)
    conData3 = spn.zoom(YY,4)
    return XX,YY,ZZ,conData,conData2,conData3


xLim=[ 0.5,1.9]
yLim=[0.01,0.8]

levels_conf=np.logspace(-10,-6,1200)
levels_bar = np.logspace(-10,-6,5)

lineColors = [str(round(el,1)) for el in np.linspace(1.0,0.0,5)]


#build figure and axes
plt.close('all')
fig, axes = plt.subplots(2, 2, sharex='col', sharey='row',figsize=(13,9))
((ax1, ax2), (ax3, ax4)) = axes



"""
plot left top
"""
XX1,YY1,ZZ1,conData1,conData12,conData13 = getData('sfTrident','dblDiff/elecDistFixEl/transInteg/phi25/a0m4/a')

CS = ax1.contourf(XX1,YY1,ZZ1,levels = levels_conf,vmin=levels_conf.min(), vmax=levels_conf.max(),norm=colors.LogNorm(),cmap='jet')

levels_conl1=np.logspace(-10,-7,5)
CSlines = ax1.contour(conData12,conData13,conData1,levels = levels_conl1, colors =lineColors)
ax1.clabel(CSlines, inline=1, fontsize=9,fmt='%1.2e')
ax1.set_xlim(xLim[0],xLim[1])
ax1.set_ylim(yLim[0],yLim[1])
#ax1.set_ylabel('$\\frac{p_T}{m_e}$')
ax1.set_title('$\\Delta \\varphi=25$',fontsize=20)

"""
plot right top
"""
XX2,YY2,ZZ2,conData2,conData22,conData23 = getData('sfTrident','dblDiff/elecDistFixEl/transInteg/phi50/a0m4/a')
CS = ax2.contourf(XX2,YY2,ZZ2,levels = levels_conf,vmin=levels_conf.min(), vmax=levels_conf.max(),norm=colors.LogNorm(),cmap='jet')

levels_conl2=np.logspace(-10,-6,5)
CSlines = ax2.contour(conData22,conData23,conData2,levels = levels_conl2, colors =lineColors)
ax2.clabel(CSlines, inline=1, fontsize=9,fmt='%1.2e')
ax2.set_xlim(xLim[0],xLim[1])
ax2.set_ylim(yLim[0],yLim[1])
ax2.set_title('$\\Delta \\varphi=50$',fontsize=20)

"""
plot left bottom
"""
XX3,YY3,ZZ3,conData3,conData32,conData33 = getData('sfTrident','dblDiff/elecDistFixEl/transInteg/phi250/a0m4/a')

CS = ax3.contourf(XX3,YY3,ZZ3,levels = levels_conf,vmin=levels_conf.min(), vmax=levels_conf.max(),norm=colors.LogNorm(),cmap='jet')

levels_conl3=np.logspace(-7,-6,5)
CSlines = ax3.contour(conData32,conData33,conData3,levels = levels_conl3, colors =lineColors)
ax3.clabel(CSlines, inline=1, fontsize=9,fmt='%1.2e')
ax3.set_xlim(xLim[0],xLim[1])
ax3.set_ylim(yLim[0],yLim[1])
ax3.set_ylabel('$\\frac{p_T}{m_e}$',fontsize=20)
ax3.set_xlabel('$y$',fontsize=20)
ax3.set_title('$\\Delta \\varphi=250$',fontsize=20)

"""
plt right bottom
"""
XX4,YY4,ZZ4,conData4,conData42,conData43 = getData('sfTrident','dblDiff/elecDistFixEl/transInteg/phi500/a0m4/a')

CS = ax4.contourf(XX4,YY4,ZZ4,levels = levels_conf,vmin=levels_conf.min(), vmax=levels_conf.max(),norm=colors.LogNorm(),cmap='jet')

levels_conl4=np.logspace(-7,-6,5)
CSlines = ax4.contour(conData42,conData43,conData4,levels = levels_conl4, colors =lineColors)
ax4.clabel(CSlines, inline=1, fontsize=9,fmt='%1.2e')
ax4.set_xlim(xLim[0],xLim[1])
ax4.set_ylim(yLim[0],yLim[1])
ax4.set_title('$\\Delta \\phi=500$',fontsize=20)


plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)


#put in colorbar (right from the plots)
cbar=fig.colorbar(CS,ticks=levels_bar,ax=axes.ravel().tolist(),shrink=0.95,format="%1.2e")

"""
f, axarr = plt.subplots(2, 2)
axarr[0, 0].plot(x, y)
axarr[0, 0].set_title('Axis [0,0]')
axarr[0, 1].scatter(x, y)
axarr[0, 1].set_title('Axis [0,1]')
axarr[1, 0].plot(x, y ** 2)
axarr[1, 0].set_title('Axis [1,0]')
#axarr[1, 1].scatter(x, y ** 2)
#axarr[1, 1].set_title('Axis [1,1]')
# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
#plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
#plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
"""
plt.show()
