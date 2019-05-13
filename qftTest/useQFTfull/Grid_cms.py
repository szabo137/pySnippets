from numpy import *
from Konstanten import *


def Grid_rap(y,w_cms,x_min,x_max,x_genau,a0,p_lim=0):
    import warnings
    warnings.filterwarnings("ignore")

    d      = lambda l,y_z: l  * w_cms * exp(y_z)  / (l*exp(2.*y_z)  + 1.)
    c      = lambda l,y_z: m**2 * a0**2 / (2.*(l*exp(2.*y_z)  + 1.))

    def p_max (l,y_z):
        p = sqrt(4.*l**2*w_cms**2/(l*exp(y_z)+exp(-y_z))**2-m**2)
        if p_lim != 0 and p > p_lim :
            p = p_lim
        return nan_to_num(p)

    def p_min (l,y_z):
        p = sqrt((d(l,y_z) + sqrt(d(l,y_z)**2 - c(l,y_z)))**2-m**2)
        return nan_to_num(p)

    p_perp = transpose([linspace(p_min(x_min,y_z),p_max(x_max,y_z),x_genau) for y_z in y])    # array y
    warnings.filterwarnings("ignore")
    return p_perp
