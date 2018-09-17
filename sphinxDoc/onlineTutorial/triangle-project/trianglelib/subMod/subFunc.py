"""
this is the test function for a sub module
"""

def wrapper(func):
    """
    this function wraps 'func'
    """
    def outFunc(args):
        return func(*args)

    return outFunc
