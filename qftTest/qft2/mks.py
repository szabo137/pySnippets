"""
class for minkowski vectors
"""

import numpy as np

def isscalar(x):
    return isinstance(x,int) or isinstance(x,float) or isinstance(x,complex)

class MinkowskiVector(object):
    def __init__(self,A):
        if (isinstance(A,list) or isinstance(A,tuple)) and len(A) == 4:
            self.arr=np.array(A)
        elif isinstance(A,np.ndarray):
            self.arr=A
        else:
            raise TypeError("<%s> is not convertable to MinkowskiVector!"%(type(A)))
        
        self.data = self.arr.data
        self.shape = self.arr.shape
    

    def __getItem__(self,item):
        return self.arr[item]
    
    def __mul__(self,other):
        if isinstance(other,MinkowskiVector):
            myVec=self()
            otherVec = other()
            #return self[0]*other[0] - self[1]*other[1] - self[2]*other[2] - self[3]*other[3]
            return myVec[0]*otherVec[0] - myVec[1]*otherVec[1] - myVec[2]*otherVec[2] - myVec[3]*otherVec[3]
        elif isscalar(other):
            return MinkowskiVector(self.arr*other)
        else:
            raise TypeError("Multiplication %s * %s is not defined."%(type(self,),type(other)))

    def __rmul__(self,other):
        return self*other

    def __div__(self,other):
        if isscalar(other):
            return self*(1.0/other)
    
    def __add__(self,other):
        return  MinkowskiVector(self.arr+other())
    
    def __sub__(self,other):
        return self+(-1.0)*other

    def __neg__(self):
        return (-1.0)*self
        
    def isonshell(self,mass,tol=1e-3):
        return abs(self*self-mass**2)<tol
    
    def __call__(self):
        return self.arr
    
    def __repr__(self):
        return self.arr.__repr__()
    
    def _0(self):
        return self.arr[0]
    
    def _1(self):
        return self.arr[1]
    
    def _2(self):
        return self.arr[2]
    
    def _3(self):
        return self.arr[3]
    
    
if __name__=='__main__':
    p=5.0
    m=3.0
    
    E=np.sqrt(p**2 + m**2)
    
    A=MinkowskiVector([E,0,0,p])
    
    print A
    
    p=6.0
    m=2.0
    
    E=np.sqrt(p**2 + m**2)
    
    B=MinkowskiVector([E,0,0,p])

    print B
    
    #add:
    print A+B
    
    #sub:
    print A-B

    #prod
    print A*B
    
    #prod2
    print B*A
    
    #prod with scalar
    print A*3.0
    print 3.0*A
    
    #neg
    print -A
    
    #isonshell
    print A.isonshell(3.0)
    
    #mass:
    print A*A
    print B*B


    import spinors as sp
    
   
