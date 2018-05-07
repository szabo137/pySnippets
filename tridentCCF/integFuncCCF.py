"""
module to compute the integrand in ccf

todo:
    -   insert into sfTrident
    -   connect to kinClass
    -   use Integrator (asimps)
    -   compare to myCalc with CCF (is it possible? -> boundarys)
"""

from phaseInt import B0ccf


def func(rStar,kinObj):
    photoNumC=rStar #insert translator
    photoNumBW = 1-photoNumC #insert BW
    alphasC = kinObj.getAlpha('c')
    alphasBW = kinObj.getAlpha('bw')
    return B0ccf(photoNumC,alphaC[0]/2.0,alphaC[2]/3.0)*B0ccf(photoNumBW,alphaBW[0]/2.0,alphaBW[2]/3.0)

def integFunc(rStar,kinObj):
    return (func(rStar,kinObj) - func(-rStar,kinObj))/rStar
