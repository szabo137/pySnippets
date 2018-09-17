"""
module to compare plot the data and shifted data in the same setting
"""

import reader as r


tempData = r.dataReader('res.dat')
rows = [1,2,3]
xshift = -4.066674294159558


(x1,y1,z1)=(tempData[rows[0]],tempData[rows[1]],tempData[rows[2]])
(x1,y1,z1)=(tempData[rows[0]]+xshift,tempData[rows[1]],tempData[rows[2]])
