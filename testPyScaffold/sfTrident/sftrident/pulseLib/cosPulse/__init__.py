"""
translates the option in the correct pulse shape and internal integrals
"""
def getInternalInt(opt=[]):
    if 'analytic' in opt:
        import analyticInt
        return analyticInt.analyticInternalInt
    else:
        import integrands
        return integrands.getIntegrals

def getEnvelope(opt):
    import pyEnvelope
    return pyEnvelope.envelope
