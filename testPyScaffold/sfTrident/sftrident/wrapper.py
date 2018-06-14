"""
contains several wrapper functions

    -   coordinate trafo of functions and derivative
"""

def funcWrapper(func,trafo):
    """Transfomes a given function.

    Parameters
    ----------
    func : python callable
        the python function which will be transformed.

    trafo: python callable
        the transformation of new coordinates to old. len(input)=len(ouput)

    Returns
    -------
    python callable
        A transformed function of the transformed parameters
    """
    def newFunc(newArgs):
        oldArgs = trafo(newArgs)
        return func(oldArgs)
    return newFunc


def diffWrapper(func,trafo,jac):
    """Transformes the derivative with given jacobian.

    Parameters
    ----------
    func :  python callable
        the python function which describes the derivative that will be transformed.

    trafo: python callable
        the transformation of new coordinates to old. len(input)=len(ouput)

    jac: python callable
        the jacobian of 'trafo' as a function of the new parameters

    Returns
    -------
    python callable
        The transformed derivative as a function of the new parameters
    """
    def newDiff(newArgs):
        oldArgs = trafo(newArgs)
        return func(oldArgs)*jac(newArgs)
    return newDiff


if __name__=='__main__':
    import numpy as np
    import matplotlib.pylab as plt

    def func(x):
        return np.sqrt(1-x**2)

    eps = 1e-8

    args = np.linspace(-1,1,100)
    vals = func(args)
    plt.plot(args,vals,label="cartesian")


    def cTrafo(theta):
        return np.cos(theta)

    tFunc = funcWrapper(func,cTrafo)

    cArgs = np.linspace(0,2.0*np.pi,100)
    cVals = tFunc(cArgs)
    plt.plot(cTrafo(cArgs),cVals,'x',label="trigonom")

    plt.show()
