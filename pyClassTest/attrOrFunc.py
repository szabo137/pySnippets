"""
module to test what is faster: Attribute or function call
"""

class testClass(object):
    def __init__(self,foo):
        self.__foo = foo
        self.outFoo = foo
    
    def getFoo(self):
        return self.__foo


if __name__=='__main__':
    import time as T
    
    A=testClass(143.3)
    
    start = T.time()
    for el in range(1000000):
        bar=A.outFoo
    end = T.time() - start
    print "attribute time: %1.2e"%(end)
    
    start = T.time()
    for el in range(1000000):
        bar=A.getFoo()
    end = T.time() - start
    print "func call time: %1.2e"%(end)
    
