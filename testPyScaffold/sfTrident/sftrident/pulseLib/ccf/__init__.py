"""
translates the option in the correct pulse shape and internal integrals
"""
def getInternalInt(opt=[]):
    import internalInt
    return internalInt.internalInt
    
def getEnvelope(opt):
    print "get ccf env"
    import pyEnvelope
    return pyEnvelope.envelope
