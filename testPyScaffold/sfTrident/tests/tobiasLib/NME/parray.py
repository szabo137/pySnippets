# -*- coding: utf-8 -*-
import numpy as np
### NME - parray #################################################################
# v1.0:	2010-12-13	Daniel Seipt



class parray(np.ndarray):
    # replacement for numpy arrays that have different behaviour
    # under multiplication
    def __new__(cls, input_array, axis=None):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        #obj.info = info
	obj.axis = axis
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self,obj):
        # reset the attribute from passed original object
        self.axis = getattr(obj, 'axis', None)
        # We do not need to return anything

    def __repr__(self):
	return str(self.data)
    
    
    def __add__(self,other):
	
	if isinstance(other,parray):
		return np.asarray(self) + other
	else:
		return other + self
	
    def __sub__(self,other):
	
	if isinstance(other,parray):
		return np.asarray(self) - other
	else:
		return -1.0*other + self
    
    
    def __mul__(self,other):
	
	if isinstance(other,parray):
		return np.asarray(self) * other
	else:
		return other * self



if __name__ == '__main__':
	
	x = parray( 1j*np.ones((2,2)) )
	
	print x*3
	print 4*x
	print x*x
	print x.__class__
	
	a = np.array([[5,-3],[2,4]])
	b = parray([[5,-3],[2,4]])
	
	z = a+x
	z = x+a
	print z
	print z.__class__
	
	#t = b-x
	#print t
	#print t.__class__

