""" example of wrapping cos function from math.h using cython"""

cdef extern from "dblcos.h":
    double dblcosFunc(double arg)

def cos_func(arg):
    return dblcosFunc(arg)

