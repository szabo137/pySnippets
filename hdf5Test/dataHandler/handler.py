"""
dev playground of the data handler
"""

import h5py
import numpy as np

class handler(object):
    def __init__(self,filePath="untitled.hdf5",mode="r"):
        self.__filePath = filePath
        self.__mode = mode

    def __enter__(self):
        self.__file = h5py.File(self.__filePath,self.__mode)
        print "file is open"
        """
        here a function recieves self.__file and gets the data structure -> sets all attributes
        """

        return self.__file

    def __exit__(self, type, value, tb):
        self.__file.close()
        print "file is closed"

if __name__=='__main__':
    with handler(mode='w') as m:
        print "here it is open"
