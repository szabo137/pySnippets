# -*- coding: utf-8 -*-
import sys,os


from numpy import *
from spinors import *

### NME - propagators #################################################################
# v1.0:	2010-12-13	Daniel Seipt


class ScalarPropagator(object):

    def __init__(self,(P,m),epsilon = 1e-14):
	if isinstance(P,MinkowskiVector):
	    self.data	= 1./(P*P -m**2 + 1j*epsilon)
        elif (isinstance(P,ndarray) or isinstance(P,list) or isinstance(P,tuple)) and len(P) == 4:
	    self.data	= 1./(P[0]**2 - P[1]**2 - P[2]**2 - P[3]**2 - m**2 + 1j*epsilon)

    def __call__(self):
	return self.data
	

class FermiPropagator(object):

    def __init__(self,(p,m),epsilon = 1e-14):

	self.data	= (feyndagg(p)+m) * ScalarPropagator((p,m))()
	
    def __call__(self):
	return self.data	

if __name__ == '__main__':
	
    m  = 1.
    px = 2.
    py = 1.1
    pz = -1.5
    E  = 3.5
    
    
    P	= [E,px,py,pz]
    
    print ScalarPropagator((P,m)).__class__
    print ScalarPropagator((P,m))().__class__
    
    
    print FermiPropagator((P,m)).__class__
    print FermiPropagator((P,m))().__class__
    