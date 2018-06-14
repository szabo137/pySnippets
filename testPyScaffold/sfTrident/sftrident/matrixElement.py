"""
calculation of the matrix element


format kinObj:

kinObj.physArea : bool
kinObj.momenta : k, p,p1,p2,p3
kinObj.photoNum : photoNumBW, photoNumC, photoNumBWx, photoNumCx
kinObj.spinors : u,ubar1,ubar2,u3
kinObj.alphas : [alphaBW],[alphaC],[alphaBWx],[alphaCx]
kinObj.photoNum: s_bw, r_c, s_bwx, r_cx <- onshell!
"""
import numpy as np
#import copy
from sfTrident.phaseIntegral import prodSum,gaussPoints # improve!!
#from sfTrident.current.currentClass import currentClass
from sfTrident.kinematics import kinClass,photoNumBW
from vertex import vertexFunctionCombine as vertexFunction

#from scipy.integrate import quad
#import sfTrident.src.qft as qft

def integCauchy(fkt,gauss=None,mode='quad'):
    """
    integrate fkt(x)/x over symetric interval [-bound,bound]
    """

    #tempFunc=lambda x: (fkt(x)-fkt(-x))/x
    tempFunc=lambda x: (fkt(1.0/(1+x))-fkt(-1.0/(1+x)))/(1+x)
    tempVals=np.array(map(tempFunc,gauss.points)).T
    outTemp=(gauss.weights*tempVals)
    return np.sum(outTemp,axis=1)
    #return np.sum(outTemp)




class matrixElement(object):
    def __init__(self,config):
        self.__kinObj = kinClass(config)
        self.__vertexObj = vertexFunction(config)



    def setKin(self,kinPara):
        """
        sets the kinematic for breit wheeler and compton vertex function
        and calculates the currents
        """
        self.__kinObj.evalKin(kinPara)
        if self.__kinObj.physArea:
            self.__vertexObj.setKin(self.__kinObj)

    """

    GENERAL MATRIX FUNCTIONS

    """


    def __setMatrixFromVertex(self,vertexA,vertexB):
        """
        builds the matrix element from vertex function

        -> source out
        """
        #return np.array([vertexA[item]*vertexB[item2] for item in np.arange(4) for item2 in np.arange(4)]).reshape(16,)
        return np.outer(vertexA,vertexB).flatten()

    def __setMatrixElementDirect(self,photoNums):
        """
        evaluates the offshell matrix element for given photonums

        returns Mdirect

        NEED TO BE IN PHYSAREA

        todo:
            -   check prefac
            -   improve prodsum
        """
        vertexBW = self.__vertexObj.getVertexBW(photoNums[0])
        vertexC = self.__vertexObj.getVertexC(photoNums[1])
        Mdirect = self.__setMatrixFromVertex(vertexBW,vertexC)
        return Mdirect


    def __setMatrixElementExchange(self,photoNums):
        """
        evaluates the offshell matrix element for given photonums

        returns Mexchange

        NEED TO BE IN PHYSAREA

        todo:
            -   check prefac
            -   improve prodsum
        """
        vertexBWx = self.__vertexObj.getVertexBWx(photoNums[0])
        vertexCx = self.__vertexObj.getVertexCx(photoNums[1])
        Mexchange =self.__setMatrixFromVertex(vertexBWx,vertexCx)
        return Mexchange


    """

    ONSHELL FUNCTIONS

    """
    def evalMdirONSHELL(self):
        """
        returns the onshell matrix element of the direct part

        (before sum)
        """
        if self.__kinObj.physArea:
            photoNumsOnshell = self.__kinObj.getPhotoNumONSHELL()

            Mdirect=self.__setMatrixElementDirect(photoNumsOnshell[[0,1]])
            #for item in np.arange(4):
            #    for item2 in np.arange(4):
            #        print "PhaseIntBW_%s: %s"%(item,self.__vertexObj.BBW[item])
            #        print "PhaseIntC_%s: %s"%(item2,self.__vertexObj.BC[item2])
            #        print "B_%s%s: %s"%(item,item2,self.__vertexObj.BBW[item]*self.__vertexObj.BC[item2])
            #print photoNumsOnshell
        return self.__kinObj.PreFac[0]*Mdirect

    def evalMexONSHELL(self):
        """
        returns the onshell matrix element of the exchange part

        (before sum)
        """
        if self.__kinObj.physArea:
            photoNumsOnshell = self.__kinObj.getPhotoNumONSHELL()
            Mexchange=self.__setMatrixElementExchange(photoNumsOnshell[[2,3]])
        return self.__kinObj.PreFac[1]*Mexchange


    def evalMsquareONSHELL(self):
        """
        returns the full onshell Msquare with sum over all spins

        -> remove this
        """
        if self.__kinObj.physArea:
            photoNumsOnshell = self.__kinObj.getPhotoNumONSHELL()
            #print "photoNum onshell: %s"%(photoNumsOnshell)
            Mdirect=self.__setMatrixElementDirect(photoNumsOnshell[[0,1]])
            #Mexchange=self.__setMatrixElementExchange(photoNumsOnshell[[2,3]])
            M=self.__kinObj.PreFac[0]*Mdirect# - self.__kinObj.PreFac[1]*Mexchange
            return self.__kinObj.PreFac[4]*prodSum(M,np.conjugate(M))
        else:
            return 0.0

    """

    OFFSHELL FUNCTIONS

    """
    def evalMdirOFFSHELL(self):
        """
        returns the offshell matrix element of the direct part

        (before sum)

        -> check preFac!!!
        """
        if self.__kinObj.physArea:
            Mdirect=self.__kinObj.PreFac[2]*self.__vertexObj.getVertexCombineDirect()
        return Mdirect

    def evalMexOFFSHELL(self):
        """
        returns the offshell matrix element of the exchange part

        (before sum)

        -> check preFac!!!
        """
        if self.__kinObj.physArea:
            Mexchange=self.__kinObj.PreFac[3]*self.__vertexObj.getVertexCombineExchange()
        return Mexchange

    """

    PARTIAL FUNCTIONS

    """
    def evalMdirPartBW(self):
        """
        returns the partial matix element with strong field direct bw part

        M_11 = J_0^C * (B_0J_0 + sum B_l*J_l)^BW

        """
        if self.__kinObj.physArea:
            photoNumPart = photoNumBW(0.0,self.__kinObj.getMomenta())
            vertexBW = self.__vertexObj.getVertexBW(photoNumPart)
            MdirBW=self.__setMatrixFromVertex(self.__vertexObj.J0c,vertexBW)
        return 2.0*np.pi*MdirBW

    def evalMdirPartC(self):
        """
        returns the partial matix element with strong field direct c part

        M_12 = J_0^BW * (B_0J_0 + sum B_l*J_l)^C

        """
        if self.__kinObj.physArea:
            photoNumPart = photoNumBW(0.0,self.__kinObj.getMomenta())
            vertexC = self.__vertexObj.getVertexC(photoNumPart)
            MdirC=self.__setMatrixFromVertex(self.__vertexObj.J0bw,vertexC)
        return 2.0*np.pi*MdirC

    def evalMexPartBW(self):
        """
        returns the partial matix element with strong field exchange bw part

        M_11x = J_0^Cx * (B_0J_0 + sum B_l*J_l)^BWx

        """
        if self.__kinObj.physArea:
            photoNumPart = photoNumBW(0.0,self.__kinObj.getMomenta())
            vertexBWx = self.__vertexObj.getVertexBWx(photoNumPart)
            MexBW=self.__setMatrixFromVertex(self.__vertexObj.J0cx,vertexBWx)
        return 2.0*np.pi*MexBW

    def evalMexPartC(self):
        """
        returns the partial matix element with strong field exchange bw part

        M_12x = J_0^BWx * (B_0J_0 + sum B_l*J_l)^Cx

        """
        if self.__kinObj.physArea:
            photoNumPart = photoNumBW(0.0,self.__kinObj.getMomenta())
            vertexCx = self.__vertexObj.getVertexCx(photoNumPart)
            MexC=self.__setMatrixFromVertex(self.__vertexObj.J0bwx,vertexCx)
        return 2.0*np.pi*MexC
    """
    part MATRIX ELEMENT
    """

    def evalMatrixElementpart(self):
        """
        evaluates the full matrix element

        TODO:
            -   benchmark: do the physArea check here and remove from the
                single Ms
        """
        Mpbw = self.evalMdirPartBW()
        Mpc = self.evalMdirPartC()
        Mpbwx = self.evalMexPartBW()
        Mpcx = self.evalMexPartC()
        return Mpbw+Mpc+Mpbwx+Mpcx

    def evalMsquarePart(self):
        if self.__kinObj.physArea:
            M = self.evalMatrixElementpart()
            return self.__kinObj.PreFac[4]*prodSum(M,np.conjugate(M))
        else:
            return 0.0


    """
    FULL MATRIX ELEMENT
    """

    def evalMatrixElement(self):
        """
        evaluates the full matrix element

        TODO:
            -   benchmark: do the physArea check here and remove from the
                single Ms
        """
        Mdon = self.evalMdirONSHELL()
        Mxon = self.evalMexONSHELL()
        Mdoff = self.evalMdirOFFSHELL()
        Mxoff = self.evalMexOFFSHELL()
        Mpbw = self.evalMdirPartBW()
        Mpc = self.evalMdirPartC()
        Mpbwx = self.evalMexPartBW()
        Mpcx = self.evalMexPartC()
        return Mdon+Mxon+Mdoff+Mxoff+Mpbw+Mpc+Mpbwx+Mpcx

    def evalMsquare(self):
        if self.__kinObj.physArea:
            M = self.evalMatrixElement()
            return self.__kinObj.PreFac[4]*prodSum(M,np.conjugate(M))
        else:
            return 0.0
