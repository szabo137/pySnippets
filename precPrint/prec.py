"""
print with full precition
"""
#from __future__ import print_function

import numpy as np
import mpmath as mp
aRaw = 10.00000000000000000000000000900000002
mp.mp.dps=50
a=mp.mp.mpf('9.00000000000000000000000000900000002')
b=mp.mp.mpf('1e26')
mp.nprint(a,50)
c=a*b

print c

d=np.sqrt(a)
print d*3

from decimal import Decimal as dec
a=dec("10.00000000000000000000000000900000002")
