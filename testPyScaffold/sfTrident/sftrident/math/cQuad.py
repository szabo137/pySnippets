"""
routine to complex quadrature using scipy quad
"""
import scipy
from scipy.integrate import quad

def cQuad(func, a, b, **kwargs):
    r"""integration of a complex valued function

    Parameters
    ----------
    func : python_callable
        scalar function to integrate (may return complex scalar)

    a : float
        lower boundary of the integration

    b : float
        upper boundary of the integration

    Returns
    -------
    result : complex
        result of the integration as a complex number

    intError : list
        numerical errors of the integration (real part, imag part)

    Other Parameters
    ----------------
    kwargs : dict (optional)
        keyword arguments for quad (see scipy.integrate.quad)

    Notes
    -----
        -   uses scipy.integrate.quad internaly
        -   executes the <func> two times per point
    """
    def real_func(x):
        return scipy.real(func(x))
    def imag_func(x):
        return scipy.imag(func(x))
    real_integral = quad(real_func, a, b, **kwargs)
    imag_integral = quad(imag_func, a, b, **kwargs)
    return (real_integral[0] + 1j*imag_integral[0], real_integral[1:], imag_integral[1:])
