"""
tests the inheritance of classes in python
"""


class A(object):
    def __init__(self,data):
        self.__data = data
        self.data=data

    def setData(self,data):
        self.__data = data

    def getData(self):
        return self.__data


class B(A):
    def __init__(self,data):
        A.__init__(self,data)
        print"test: %s"%self.data

    def transformData(self):
        self._A__data = self._A__data*2


if __name__=='__main__':
    testA = A("bla")
    print "data: %s"%(testA.getData())


    test = B("bla")
    print "data: %s"%(test.getData())
    test.transformData()
    print "transformed data: %s"%test.getData()

    print testA.getData()
