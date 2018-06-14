"""
this class provides the vertexFunctions
"""
import numpy as np
from sfTrident.phaseIntegral import phaseIntegral as phaseIntegral
from sfTrident.phaseIntegral import prodSum
from sfTrident.current import currentClass
from sfTrident.kinematics import photoNumBW, rFromRstar
from sfTrident.math.cQuad import cQuad


def ProdSingle(arr1,arr2):
    return np.multiply.outer(arr1,arr2).flatten()

def Prod(arr1,arr2):
    return np.array([ProdSingle(arr1[item],arr2[item2]) for item in np.arange(4) for item2 in np.arange(4) ])


class vertexFunction(object):
    """
    class to calculate the single VertexFunctions
    """
    def __init__(self,config):
        """
        evaluation of the vertex function

        need to be in physArea!
        """

        self.config = config #avoid after testing
        self.phaseIntObj = phaseIntegral(config,deg = config['deg']) #generalize the config! deg -> config
        self.__currentObj = currentClass(config)

    def setKin(self,kinObj):
        """
            kinObj needs to be evaluated!

            todo:
                -   proof concept, if kinObj is evaluated?!
        """
        self.kinObj = kinObj
        self.__setCurrents()
        #print self.J0c



    def __evalSingleCurrent(self,mode):
        """
        evaluation of currents on the given kinObj

        need to be in physArea
        """
        self.__currentObj.setKin(self.kinObj.getMomenta(mode)[0],self.kinObj.getMomenta(mode)[1],self.kinObj.getMomenta(mode)[2],self.kinObj.getSpinors(mode)[0],self.kinObj.getSpinors(mode)[1])
        return self.__currentObj.J0,self.__currentObj.J1,self.__currentObj.J2,self.__currentObj.J3


    def __setCurrents(self):
        """
        sets currents for every mode

        need to be in physArea!!!
        """
        self.J0bw,self.J1bw,self.J2bw,self.J3bw = self.__evalSingleCurrent('bw')
        self.J0c,self.J1c,self.J2c,self.J3c = self.__evalSingleCurrent('c')
        self.J0bwx,self.J1bwx,self.J2bwx,self.J3bwx = self.__evalSingleCurrent('bwx')
        self.J0cx,self.J1cx,self.J2cx,self.J3cx = self.__evalSingleCurrent('cx')
        test = ProdSingle(self.J1bw,self.J0c)
        #for el in test:
        #    print np.asarray(el)


    def evalPhaseIntOFFSHELL(self,alphas,photoNums):
        """
        calculates all phaseInts for given alpha and photoNum
        """
        self.phaseIntObj.setPhotoNum(photoNums)
        self.phaseIntObj.setKin(alphas)
        return [self.phaseIntObj.getB1(),self.phaseIntObj.getB2(),self.phaseIntObj.getB3(),self.phaseIntObj.getB0()]


    def getVertexBW(self,photoNum):
        """
        returns the single BW-Vertex for given photoNum
        """
        B1bw,B2bw,B3bw,B0bw=self.evalPhaseIntOFFSHELL(self.kinObj.getAlpha('bw'),photoNum)
        self.BBW = np.array([B0bw,B1bw,B2bw,B3bw])
        return B0bw*self.J0bw+B1bw*self.J1bw+B2bw*self.J2bw+B3bw*self.J3bw

    def getVertexC(self,photoNum):
        """
        returns the single C-Vertex for given photoNum
        """
        B1c,B2c,B3c,B0c=self.evalPhaseIntOFFSHELL(self.kinObj.getAlpha('c'),photoNum)
        self.BC = np.array([B0c,B1c,B2c,B3c])
        return B0c*self.J0c+B1c*self.J1c+B2c*self.J2c+B3c*self.J3c

    def getVertexBWx(self,photoNum):
        """
        returns the single BWx-Vertex for given photoNum
        """
        B1bwx,B2bwx,B3bwx,B0bwx=self.evalPhaseIntOFFSHELL(self.kinObj.getAlpha('bwx'),photoNum)
        return B0bwx*self.J0bwx+B1bwx*self.J1bwx+B2bwx*self.J2bwx+B3bwx*self.J3bwx

    def getVertexCx(self,photoNum):
        """
        returns the single Cx-Vertex for given photoNum
        """
        B1cx,B2cx,B3cx,B0cx=self.evalPhaseIntOFFSHELL(self.kinObj.getAlpha('cx'),photoNum)
        return B0cx*self.J0cx+B1cx*self.J1cx+B2cx*self.J2cx+B3cx*self.J3cx

class vertexFunctionCombine(vertexFunction):
    """
    This class models the combined vertexFunctions Delta_BW*Delta_C.
    """
    def __init__(self,config):
        vertexFunction.__init__(self,config)



    def getCurrentCombineDirect(self):
        """
        returns the combine funtion J_ij for any combination of spins, direct part

        need to evaluate the kinematic first
        """
        #self.Jdir = Prod([self.J0bw,self.J1bw,self.J2bw,self.J3bw],[self.J0c,self.J1c,self.J2c,self.J3c])
        self.JBW = [self.J0bw,self.J1bw,self.J2bw,self.J3bw]
        self.JC = [self.J0c,self.J1c,self.J2c,self.J3c]
        #self.Jdir = np.array([ProdSingle(self.JBW[item],self.JC[item2]) for item in np.arange(4) for item2 in np.arange(4)])
        #self.Jdir = np.array([ProdSingle(self.JBW[item],self.JC[item2]) for item in np.arange(4) for item2 in np.arange(4)])
        return 0

    def getCurrentCombineExchange(self):
        """
        returns the combine funtion J_ij for any combination of spins, exchange part

        need to evaluate the kinematic first
        """
        self.JBWx = np.array([self.J0bwx,self.J1bwx,self.J2bwx,self.J3bwx])
        self.JCx = np.array([self.J0cx,self.J1cx,self.J2cx,self.J3cx])
        #self.Jex = Prod(self.JBWx,self.JCx)
        return 0

    """

    Phase Integral Combine: direct part

    """

    def __integrandPhaseIntegralDirect(self,item1,item2,rStar):
        """
        returns B_ij for a given set of photoNums

        check order!!

        -> direct part!!!
        """
        #print "i1: %s"%item1
        #print "i2: %s"%item2
        pNumC=rFromRstar(rStar,self.kinObj.getMomenta('c'))
        #print "rstar: %s"%rStar
        #print "pNumC: %s"%pNumC
        self.phaseIntObj.setKin(self.kinObj.getAlpha('c'))
        phaseC=[self.phaseIntObj.B0func,self.phaseIntObj.B1func,self.phaseIntObj.B2func,self.phaseIntObj.B3func]
        phaseIntC = phaseC[item2](pNumC)
        #print"phaseIntC_%s: %s"%(item2,phaseIntC)

        self.phaseIntObj.setKin(self.kinObj.getAlpha('bw'))
        phaseBW=[self.phaseIntObj.B0func,self.phaseIntObj.B1func,self.phaseIntObj.B2func,self.phaseIntObj.B3func]
        pNumBW = photoNumBW(pNumC,self.kinObj.getMomenta())
        #print "pNumBW: %s"%pNumBW
        phaseIntBW = phaseBW[item1](pNumBW)
        #print"phaseIntBW_%s: %s"%(item1,phaseIntBW)
        return phaseIntBW*phaseIntC


    def integFuncDirect(self,item1,item2,rStar):
        """
        integrand after trafo r to rStar
        """
        #return self.__integrandPhaseIntegralDirect(item1,item2,rStar)
        return (self.__integrandPhaseIntegralDirect(item1,item2,rStar)-self.__integrandPhaseIntegralDirect(item1,item2,-rStar))/rStar

    def getPhaseIntCombineDirect(self,item1,item2):
        """
        returns B_item1_item2

        todo:
            -   move singular points calc to kinObj
            -   prefac!

        -> direct part!!!
        """
        #test
        #photoNumsOnshell = self.kinObj.getPhotoNumONSHELL('c')
        #print photoNumsOnshell
        #code
        func=lambda x: self.integFuncDirect(item1,item2,x)


        mom=self.kinObj.getMomenta('c')
        dP = -mom[2]-mom[1]
        #rStar = dP*dP + 2.0*dP*mom[0]*photoNumsOnshell
        #return self.kinObj.PreFac[0]*func(rStar)
        mom2=self.kinObj.getMomenta()
        pt = mom2[2] + mom2[3] + mom2[4]
        pkt =  (pt*pt - 1)*(dP*mom2[0])/((-mom2[0]*mom2[1])) + (dP*dP)
        a1=1e-10 #generalize!!
        a2=80.0#np.maximum(np.float(-dP*dP)+eps,np.float(pkt)+eps,np.array([15.0]))
        #print"sing1: %s"%(np.float(-dP*dP))
        #print"sing2: %s"%(np.float(pkt))
        return cQuad(func,a1, a2,points=(np.float(-dP*dP),np.float(pkt)))[0]


    def getVertexCombineDirect(self):
        """
        returns the combined vertexFunctions:

            M   =   Delta_BW*Delta_C
                =   sum_ij B_ij J_ij
        """
        self.getCurrentCombineDirect()
        #print self.Jdir.shape
        #for el in self.Jdir:
        #    print el.shape
        #Bdir = np.array([self.getPhaseIntCombineDirect(item,item2) for item in np.arange(4) for item2 in np.arange(4) ])
        #for item in np.arange(4):
        #    for item2 in np.arange(4):
        #        print "B_%s%s: %s"%(item,item2,self.getPhaseIntCombineDirect(item,item2))
        #print Bdir
        #for index in np.arange(16):
        #    print "index: %s"%index
        #    for el in self.Jdir[index]:
        #        print np.asarray(el)
        """
        TODO:
            -   check order of Jdir!!!
        """
        #for el in np.sum(np.array([Bdir[index]*self.Jdir[index] for index in np.arange(16)]),axis=1):
        #    print np.asarray(el)
        res = np.sum(np.array([self.getPhaseIntCombineDirect(item,item2)*ProdSingle(self.JBW[item],self.JC[item2]) for item in np.arange(4) for item2 in np.arange(4)]),axis=0)
        #print res.shape
        return res

    """

    Phase Integral Combine: exchange part

    """


    def __integrandPhaseIntegralExchange(self,item1,item2,rStar):
        """
        returns B_ij for a given set of photoNums

        check order!!

        -> exchange part!!!
        """
        pNumC=rFromRstar(rStar,self.kinObj.getMomenta('cx'))
        self.phaseIntObj.setKin(self.kinObj.getAlpha('cx'))
        phaseC=[self.phaseIntObj.B0func,self.phaseIntObj.B1func,self.phaseIntObj.B2func,self.phaseIntObj.B3func]
        phaseIntC = phaseC[item2](pNumC)
        #print phaseIntC

        self.phaseIntObj.setKin(self.kinObj.getAlpha('bwx'))
        phaseBW=[self.phaseIntObj.B0func,self.phaseIntObj.B1func,self.phaseIntObj.B2func,self.phaseIntObj.B3func]
        pNumBW = photoNumBW(pNumC,self.kinObj.getMomenta())

        phaseIntBW = phaseBW[item1](pNumBW)
        #print phaseIntBW
        return phaseIntBW*phaseIntC


    def integFuncExchange(self,item1,item2,rStar):
        """
        integrand after trafo r to rStar
        """
        #return self.__integrandPhaseIntegralExchange(item1,item2,rStar)
        return (self.__integrandPhaseIntegralExchange(item1,item2,rStar)-self.__integrandPhaseIntegralExchange(item1,item2,-rStar))/rStar

    def getPhaseIntCombineExchange(self,item1,item2):
        """
        returns B_item1_item2

        todo:
            -   move singular points calc to kinObj
            -   prefac!

        -> exchange part!!!
        """
        #photoNumsOnshell = self.kinObj.getPhotoNumONSHELL('cx')
        #print photoNumsOnshell
        #code
        #func=lambda x: self.integFuncExchange(item1,item2,x)


        #mom=self.kinObj.getMomenta('cx')
        #dP = -mom[2]-mom[1]
        #rStar = dP*dP + 2.0*dP*mom[0]*photoNumsOnshell
        #return self.kinObj.PreFac[1]*func(rStar)
        func=lambda x: self.integFuncExchange(item1,item2,x)
        mom=self.kinObj.getMomenta('cx')
        dP = -mom[2]-mom[1]
        mom2=self.kinObj.getMomenta()
        pt = mom2[2] + mom2[3] + mom2[4]
        pkt =  (pt*pt - 1)*(dP*mom2[0])/((-mom2[0]*mom2[1])) + (dP*dP)
        a1=1e-10 #generalize!!
        a2=80.0
        return cQuad(func,a1, a2,points=(np.float(-dP*dP),np.float(pkt)))[0]

    def getVertexCombineExchange(self):
        """
        returns the combined vertexFunctions:

            M   =   Delta_BW*Delta_C
                =   sum_ij B_ij J_ij
        """
        self.getCurrentCombineExchange()
        #res = np.sum(np.array([Bex[index]*self.Jex[index] for index in np.arange(16)]),axis=1)
        res = np.sum(np.array([self.getPhaseIntCombineExchange(item,item2)*ProdSingle(self.JBWx[item],self.JCx[item2]) for item in np.arange(4) for item2 in np.arange(4)]),axis=0)
        return res
