import numpy as np
import mks
import spinors as sp

p=5.0
m=3.0
    
E=np.sqrt(p**2 + m**2)
    
A=mks.MinkowskiVector([E,0,0,p])

u=sp.SpinorU((A,3.0),1)
print u

uBar=sp.SpinorUBar((A,3.0),1)
print uBar

print uBar*u
