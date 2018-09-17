"""
module to read data from file
"""
import numpy as np
import csv




def dataReader(filename):
    with open(str(filename),"rb") as f:
        data = csv.reader(f,delimiter = " ")
        tempData = np.array([np.array(t) for t in data])
        tempData = np.array([el[el!=''] for el in tempData])#avoiding whitespace
        temp=np.zeros(tempData.shape[::-1])
        for rowInd in np.arange(len(tempData)):
            row=tempData[rowInd]#[:-1]
            for index in np.arange(len(temp)):
                temp[index][rowInd] = row[index]
        return list(temp)


def dataReaderPaper(filename):
    with open(str(filename),"rb") as f:
        data = csv.reader(f,delimiter = " ")
        tempData = np.array([t for t in data])
        temp=np.zeros(tempData.shape[::-1])
        for rowInd in np.arange(len(tempData)):
            row=tempData[rowInd]
            for index in np.arange(len(temp)):
                temp[index][rowInd] = row[index]
        return list(temp)


def buildMesh2(xt,yt,zt,nx=0,ny=0):
    """
    rebuild for non-quadratic shapes
    """
    if nx:
        tempNx = nx
    else:
        tempNx = np.sqrt(len(x))
        if not tempNx.is_integer():
            raise "Automatic search of dimension fail! Give right dimension."
        else:
            tempNx = int(tempNx)
    if ny:
        tempNy = ny
    else:
        tempNy = np.sqrt(len(y))
        if not tempNy.is_integer():
            raise "Automatic search of dimension fail! Give right dimension."
        else:
            tempNy = int(tempNy)
    return xt.reshape(tempNx,tempNy).T,yt.reshape(tempNx,tempNy).T,zt.reshape(tempNx,tempNy).T


def buildMesh(xt,yt,zt,n=0):
    """
    rebuild for non-quadratic shapes
    """
    if n:
        tempN = n
    else:
        tempN = np.sqrt(len(xt))
        if not tempN.is_integer():
            raise "Automatic search of dimension fail! Give right dimension."
        else:
            tempN = int(tempN)
    return xt.reshape(tempN,tempN).T,yt.reshape(tempN,tempN).T,zt.reshape(tempN,tempN).T
