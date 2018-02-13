"""
test error handling in MCMISER
"""
from skmonaco import mcmiser

def test(x):
    raise TypeError,"blabla FEHLER blabla"
    return 0

try:
    print mcmiser(test,xl=[0.0],xu=[1.0],npoints=100)
except TypeError:
    print "schade"
