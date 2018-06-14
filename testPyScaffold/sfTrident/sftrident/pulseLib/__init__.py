"""
choses the pulse from and some optional args (like analytics)
"""

def getPulse(pulseName,options=[]):
    if pulseName=='cos':
        import cosPulse
        return (cosPulse.getEnvelope(options),cosPulse.getInternalInt(options))
    elif pulseName=='ccf':
        import ccf
        return (ccf.getEnvelope(options),ccf.getInternalInt(options))
    else:
        raise AttributeError("There is no pulse named: <%s>"%pulseName)
