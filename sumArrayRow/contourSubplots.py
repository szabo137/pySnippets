import numpy as np
import matplotlib.pyplot as plt

# Simple data to display in various forms
x = np.linspace(-2 * np.pi, 2 * np.pi, 400)
y = np.linspace(-2 * np.pi, 2 * np.pi, 400)
#y = np.sin(x ** 2)

plt.close('all')


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Create 2x2 sub plots
gs = gridspec.GridSpec(2, 2)


XX,YY = np.meshgrid(x,y)
ZZ=XX**2 + YY**2


fig = plt.figure()
ax1 = fig.add_subplot(gs[0, 0]) # row 0, col 0

#ax1.plot(x,y)
ax1.contourf(XX,YY,ZZ)
ax1.set_xlim([-2,2])
ax1.set_title('axis [0,0]')

ax2 = fig.add_subplot(gs[0, 1]) # row 0, col 1
ax2.scatter(x,y)
ax2.set_title('axis [0,1]')

ax3 = fig.add_subplot(gs[1, 0]) # row 1, span all columns
ax3.plot(x,y**2)
ax3.set_title('axis [1,0]')

plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
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
