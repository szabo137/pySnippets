import numpy as np
import wwqCalc.qft as qft
m  = 1.3
px = 2.23
py = 1.1
pz = -1.5
E  = np.sqrt(m**2 + px**2 + py**2 + pz**2)
P= qft.MinkowskiVector([E,px,py,pz])
