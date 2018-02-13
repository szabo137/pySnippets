import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.gridspec as gridspec
from itertools import imap

# Simple data to display in various forms
def fkt(a,b):
    return a**2 + b**2


x = np.linspace(0,10, 100)
y = np.linspace(-10,10, 100)

XX,YY = np.meshgrid(x,y)
ZZ=fkt(XX,YY)


plt.close('all')


def flattenContour(Z,axis='row',func=np.sum):
    if axis=='row':
        return np.array(map(func,Z))
    elif axis=='col':
        return np.array(map(func,Z.T))
    else:
        raise ValueError("axis need to be 'row' or 'col'. ('%s' given)"%axis)



# Create 2x2 sub plots
gs = gridspec.GridSpec(2, 3,width_ratios=[1, 3,0.3],height_ratios=[3, 1])



fig = plt.figure()
#ax1.set_ybound(3.0)


ax2 = fig.add_subplot(gs[1, 1]) # row 0, col 1
axX_X=XX[0,:]
#axX_Y=np.array(map(np.sum,ZZ.T))
axX_Y=flattenContour(ZZ,axis='col')
ax2.plot(axX_X,axX_Y)
ax2.set_ylabel('Xy-label')
ax2.set_xlabel('Xx-label')
#ax2.set_title('axis [1,1]')

ax3 = fig.add_subplot(gs[0, 0]) # row 1, span all columns
axY_X=YY[:,0]
axY_Y=np.array(map(np.sum,ZZ))
axY_Y=flattenContour(ZZ,axis='row')
#ax3.plot(axY_Y,axY_X)
ax3.set_ylabel('Yx-label')
ax3.set_xlabel('Yy-label')
ax3.plot(axY_Y,axY_X)
ax3.xaxis.set_major_formatter(ScalarFormatter())
for el in dir(ax3):print el
#ax3.set_title('axis [0,0]')

ax0=fig.add_subplot(gs[0,2])
ax1 = fig.add_subplot(gs[0, 1],sharex=ax2,sharey=ax3) # row 0, col 0
CS = ax1.contourf(XX,YY,ZZ)
ax1.set_xlim([0,10])
ax1.set_ylim([-10,10])
ax1.set_ylabel('y-label')
ax1.set_xlabel('x-label')
#ax1.set_title('axis [0,1]')

barlines = np.linspace(ZZ.min(),ZZ.max(),6)
cbar=fig.colorbar(CS,ticks=barlines,format="%1.2e",cax=ax0)
fig.tight_layout()
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.show()
