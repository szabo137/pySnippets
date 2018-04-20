"""
this module contains special matrices
"""

import numpy as np

gamma0	= np.array([ [ 1, 0, 0, 0],
                    [ 0, 1, 0, 0],
                    [ 0, 0,-1, 0],
                    [ 0, 0, 0,-1]])
	
gamma1	= np.array([ [ 0, 0, 0,-1],
                    [ 0, 0,-1, 0],
                    [ 0, 1, 0, 0],
                    [ 1, 0, 0, 0]])
			  
gamma2	= np.array([ [  0, 0 , 0 , 1j],
                    [  0, 0 ,-1j,  0],
                    [  0,-1j, 0 ,  0],
                    [ 1j, 0 , 0 ,  0]])
			  
gamma3	= np.array([ [ 0, 0,-1, 0],
                    [ 0, 0, 0, 1],
                    [ 1, 0, 0, 0],
                    [ 0,-1, 0, 0]])

#gamma matrices (low lorentz index -> contravariant)

gamma=np.array([gamma0,gamma1,gamma2,gamma3])



#metric
metric = np.diag([1,-1,-1,-1])
