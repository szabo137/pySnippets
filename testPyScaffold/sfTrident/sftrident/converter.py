"""
converts constants
"""

import math as m

def alpha2charge(alpha):
    return m.sqrt(4.0*m.pi*alpha)

def charge2alpha(charge):
    return charge**2/(4.0*m.pi)
