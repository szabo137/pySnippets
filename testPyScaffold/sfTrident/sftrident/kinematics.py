# -*- coding: utf-8 -*-
"""
The kinematics module
=====================

alphaClass
----------

contains class of kinematic prefactors alpha1,alpha2,alpha3

todo:
    -   scaling of polarisation! -> a_i!=pol_i

kinClass
--------

contains all kinematic for trident process

getKin Format:

[physArea,np.array(allMom, len=5),preFac]
momenta = k,-p,p1,p2,p3
spinors = u,ubar1,ubar2,v3

todo:
    -   improve the __evalAlphas function
    -   set momenta on (k,p,p1,p2,p3) -> check eval spinors
    -   set default for config
    -   build function: getMomenta, getSpinor, ... as functions of mode, which
        returns only the quantity of the given mode -> the status of kinObj will be unchanged

maybe: set -p as a momentum in kinClass or check where the - need to be!!!
"""

import numpy as np
import qft
import lightCone
from kinutility import *

modeDict={
    'bw': [[0,3,4],[2,3],0],#k,p2,p3; ubar2,v3
    'c': [[0,2,1],[1,0],1],#k,p1,-p; ubar1,u
    'bwx':[[0,2,4],[1,3],2],#k,p1,p3; ubar1,v3
    'cx':[[0,3,1],[2,0],3] # k,p2,-p; ubar2,u
}




pol1=BGpolarisationBase(1)
pol2=BGpolarisationBase(2)

class alpha(object):
    r"""Calculation of kinematic factors.

    Parameters
    ----------
    mom : array_like
        A 3-array which contains the momenta as qft.minkowskiVector.
        The format is mom = [k,pa,pb].
    config : dict
        configuration dictionary which contains at least a0 and mass.

    Returns
    -------
    python callable
        A function which provides the alphas (internal terms are saved).

    Notes
    -----
    Calculation of

    .. math:: \alpha _{1,2} = 2ma_0\left(\frac{p_a \epsilon_{1,2}}{2 k p_a} - \frac{p_b \epsilon_{1,2}}{2 k p_b}\right)

    and

    .. math:: \alpha_3 = m^2a_0^2\left(\frac{1}{2kp_a} + \frac{1}{2kp_b}\right).

    Signs of $p_a, p_b$ is chosen from Breit-Wheeler-Part.
    """
    def __init__(self,mom,config):
        #print "pol1: %s"%pol1
        #print "pol2: %s"%pol2
        self.__a0 = config['a0']
        self.__me = config['mass']
        self.__pak=mom[0]*mom[1]
        self.__pbk=mom[0]*mom[2]
        #print"pak: %s"%self.__pak
        #print"pbk: %s"%self.__pbk
        self.__pa=mom[1]
        #print "pa: %s"%self.__pa
        self.__pb=mom[2]
        #print "pb: %s"%self.__pb

        #print"term1 a1: %s"%((self.__pa*pol1)/(2.0*self.__pak))
        #print"term2 a1: %s"%((self.__pb*pol1)/(2.0*self.__pbk))

        #print"term1 a2: %s"%((self.__pa*pol2)/(2.0*self.__pak))
        #print"term2 a2: %s"%((self.__pb*pol2)/(2.0*self.__pbk))

    def __call__(self,item):
        r""" function call of alpha object

        Parameters
        ----------
        item : integer
            index of alpha

        Returns
        -------
        float
            the value of alpha(item)
        """
        if item==1:
            #print "calc alpha1"
            #print "a_1: %s"%(2.0*((self.__pa*pol1)/(2.0*self.__pak) - (self.__pb*pol1)/(2.0*self.__pbk))*self.__a0*self.__me)
            return 2.0*((self.__pa*pol1)/(2.0*self.__pak) - (self.__pb*pol1)/(2.0*self.__pbk))*self.__a0*self.__me
        elif item==2:
            #print "calc alpha2"
            #print "a_2: %s"%(2.0*((self.__pa*pol2)/(2.0*self.__pak) - (self.__pb*pol2)/(2.0*self.__pbk))*self.__a0*self.__me)
            return 2.0*((self.__pa*pol2)/(2.0*self.__pak) - (self.__pb*pol2)/(2.0*self.__pbk))*self.__a0*self.__me
        elif item==3:
            #print "calc alpha3"
            return self.__me**2*self.__a0**2*(1.0/(2.0*self.__pak) + 1.0/(2.0*self.__pak))
        else:
            raise ValueError("The index for alpha needs to be 1, 2 or 3! (<%s> given)"%item)


class kinClass(object):
    def __init__(self,config=None):
        self.__config = config
        self.__kinOutInit()
        self.__kinLib = lightCone.getKin #generalize!

    def __kinOutInit(self):
        self.physArea=False
        self.__allMomenta=None
        self.PreFac=None
        self.__allSpinors= None
        self.__allAlphas = None
        self.__allPhotoNumONSHELL = None


    def evalKin(self,kinPara):
        self.__calcMomenta(kinPara)
        if self.physArea:
            self.__calcSpinors()
            self.__evalAlphas()
            self.__evalPhotoNumONSHELL()
        else:
            self.__kinOutInit()

    def __calcMomenta(self, kinPara):
        """
        calculates all necessary momenta and set them to the spinor attr
        """
        self.physArea,self.__allMomenta,self.PreFac = self.__kinLib(kinPara)
        return 0

    def getMomenta(self,mode='full'):
        """
        returns selected momenta
        """
        if mode=='full':
            return self.__allMomenta
        else:
            return self.__allMomenta[modeDict[mode][0]]


    def __calcSpinors(self):
        """
        calculates all necessary spinors and set them to the spinor attr

        only evaluate if physArea==True
        """
        self.__allSpinors = self.__evalSpinors(self.__allMomenta)

    def __evalSpinors(self,momenta):
        """
        evals all necesarry spinors

        todo:
            -   move this func to src.qft


            CHECK THIS -> momentum change!


        """
        evalU=lambda mom: np.array([qft.SpinorU((mom,1.0),1),qft.SpinorU((mom,1.0),2)])
        evalUbar=lambda mom: np.array([qft.SpinorUBar((mom,1.0),1,eigenspinor = 'helicity'),qft.SpinorUBar((mom,1.0),2,eigenspinor = 'helicity')])
        evalV=lambda mom: np.array([qft.SpinorV((mom,1.0),1),qft.SpinorV((mom,1.0),2)])
        evalAll=lambda allMom: np.array([evalU(-allMom[0]),evalUbar(allMom[1]),evalUbar(allMom[2]),evalV(allMom[3])])
        return evalAll(momenta[1:])

    def getSpinors(self,mode='full'):
        """
        returns selected spinors
        """
        if mode=='full':
            return self.__allSpinors
        else:
            return self.__allSpinors[modeDict[mode][1]]

    def __evalAlphas(self):
        """
        evaluation of all alpha

        check order and sign at C/Cx !!!
        """
        #breit wheeler
        self.__alphaObjBW = alpha(self.getMomenta('bw'),self.__config)
        self.__alphaBW = [self.__alphaObjBW(index) for index in [1,2,3]]
        #compton
        self.__alphaObjC = alpha(self.getMomenta('c'),self.__config)
        self.__alphaC = [self.__alphaObjC(index) for index in [1,2,3]]
        #breit wheeler exchange
        self.__alphaObjBWx = alpha(self.getMomenta('bwx'),self.__config)
        self.__alphaBWx = [self.__alphaObjBWx(index) for index in [1,2,3]]
        #compton exchange
        self.__alphaObjCx = alpha(self.getMomenta('cx'),self.__config)
        self.__alphaCx = [self.__alphaObjCx(index) for index in [1,2,3]]
        self.__allAlphas = [self.__alphaBW,self.__alphaC,self.__alphaBWx,self.__alphaCx]

    def getAlpha(self,mode='full'):
        """
        return selected alpha
        """
        if mode=='full':
            return self.__allAlphas
        else:
            return self.__allAlphas[modeDict[mode][2]]

    def __evalPhotoNum(self,photonNumC,photonNumCx):
        """
        calc the photoNums for abritrary compton photonNum
        needs evalMomenta first!
        """
        tempphotoNumBW = photoNumBW(photonNumC,self.__allMomenta)
        tempphotoNumBWx = photoNumBW(photonNumCx,self.__allMomenta)
        return np.array([tempphotoNumBW,photonNumC,tempphotoNumBWx,photonNumCx])

    def __evalPhotoNumONSHELL(self):
        """
        calculation of photoNums with onhell intermediate photon

        needs evalMomenta first!
        specialized on two step process
        """
        photoNumC = photoNumConshell(self.getMomenta('c'))
        photoNumCx = photoNumConshell(self.getMomenta('cx'))

        self.__allPhotoNumONSHELL = self.__evalPhotoNum(photoNumC,photoNumCx)

    def getPhotoNumONSHELL(self,mode='full'):
        """
        returns selected photoNums
        """

        if mode=='full':
            return self.__allPhotoNumONSHELL
        else:
            return self.__allPhotoNumONSHELL[modeDict[mode][2]]

    def setMode(self,mode):
        raise AttributeError("There is no method setMode anymore")
