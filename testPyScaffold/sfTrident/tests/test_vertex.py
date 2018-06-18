"""
module to test the bw vertex against tobias results
"""
import unittest
import numpy as np
import time

import sftrident.qft as qft
from sftrident.vertex import vertexFunction
from sftrident.kinematics import alpha

from tobiasLib.Lambda import Lambda_2 as MatTob
m=1.0
a0=1.0
#xi=np.pi/4.0
theta_ph      = 0.
phi_ph        = 0.

def eps_photon(l):
    """
    single photon polariation
    """
    return (l==1)*qft.MinkowskiVector([0,np.cos(theta_ph)*np.cos(phi_ph),np.cos(theta_ph)*np.sin(phi_ph),-np.sin(theta_ph)])\
           +(l==2)*qft.MinkowskiVector( [0,-np.sin(phi_ph),np.cos(phi_ph),0])

def eps_laser(l):
    """
    laser polarisation
    """
    return (l==1)*qft.MinkowskiVector( [0,1,0,0] ) + (l==2)*qft.MinkowskiVector( [0,0,1,0] )



spinDict=[[0,0],[1,1],[0,1],[1,0]]

class testVertex(vertexFunction):
    def __init__(self,config):
        vertexFunction.__init__(self,config)

    def testGetVertex(self, momenta,spinors,alphas,photoNum,verba=False):
        """
        only for testing purpose
        """
        self.polPhoto = np.array([eps_photon(l) for l in [1,2]])
        self.xi = self.config['xi']
        self._vertexFunction__currentObj.setKin(momenta[0],momenta[1],momenta[2],spinors[0],spinors[1])
        testJ0,testJ1,testJ2,testJ3 = self._vertexFunction__currentObj.J0,self._vertexFunction__currentObj.J1,self._vertexFunction__currentObj.J2,self._vertexFunction__currentObj.J3
        testB1bw,testB2bw,testB3bw,testB0bw=self.evalPhaseIntOFFSHELL(alphas,photoNum)
        if verba:
            print "-------------myCurrents-------------"

            for item,el in enumerate(testJ0):
                print "J0 %s: %s"%(spinDict[item],el*self.polPhoto[0])


            testJplus = 0.5*(np.cos(self.xi)*testJ1 - 1j*np.sin(self.xi)*testJ2)
            for item,el in enumerate(testJplus):
                print "Jplus %s: %s"%(spinDict[item],el*self.polPhoto[0])


            testJminus = 0.5*(np.cos(self.xi)*testJ1 + 1j*np.sin(self.xi)*testJ2)
            for item,el in enumerate(testJminus):
                print "Jminus %s: %s"%(spinDict[item],el*self.polPhoto[0])

            testJ20 = 0.5*testJ3
            for item,el in enumerate(testJ20):
                print "J20 %s: %s"%(spinDict[item],el*self.polPhoto[0])



            print "------------- my phase integrals --------------"
            print "B0: %s"%testB0bw
            print "B1: %s"%testB1bw
            print "B2: %s"%testB2bw
            print "B3: %s"%testB3bw
        return (testB0bw*testJ0+testB1bw*testJ1+testB2bw*testJ2+testB3bw*testJ3,momenta,spinors,photoNum)



class mytest(unittest.TestCase):
    def setUp(self):
        #tobias
        start = time.time()
        self.tobRes = MatTob()
        end=time.time()-start
        print "tob time: %1.2e"%end
        self.MatTobRes = self.tobRes[0]
        #print self.MatTobRes
        #print "pPos tobias: %s \n\t(type: %s)"%(self.tobRes[1].data,type(self.tobRes[1].data))
        self.pPos = qft.MinkowskiVector(np.asarray(self.tobRes[1].data))

        #print "pPos: %s \n\t(type: %s)"%(self.pPos,type(self.pPos))
        #print qft.SpinorV((self.pPos,m),1)
        self.pEl = qft.MinkowskiVector(np.asarray(self.tobRes[2].data))
        #print "pEl: %s \n\t(type: %s)"%(self.pEl,type(self.pEl))
        self.pLaser = qft.MinkowskiVector(np.asarray(self.tobRes[3].data))
        #print "pLaser: %s \n\t(type: %s)"%(self.pLaser,type(self.pLaser))
        self.pPhoton = qft.MinkowskiVector(np.asarray(self.tobRes[4].data))
        #print "pPhoton: %s \n\t(type: %s)"%(self.pPhoton,type(self.pPhoton))

        self.bwMomenta = np.array([self.pLaser,self.pEl,self.pPos])

        #spinors
        self.V = np.array([qft.SpinorV ((self.pPos, m),s1) for s1 in [1,2]])
        self.Ubar  = np.array([qft.SpinorUBar((self.pEl,m),s1) for s1 in [1,2]])
        self.spinors = np.array([self.Ubar,self.V])
        #config
        self.config = {'a0':1.0,'mass':1.0,'xi':np.pi/4.0,'dPhi':25,'envelope':'cos','pulseOpt':['analytic'],'deg':350}


        #extra
        self.eps=1e-6
        self.polPhoto = np.array([eps_photon(l) for l in [1,2]])



    def test1(self):
        """
        tests vertex*pol mod squared and spin summed of breit wheeler against tobias result
        """
        start = time.time()
        self.vertexObj = testVertex(self.config)
        self.photoNum =  (self.pPos.minus() + self.pEl.minus() - self.pPhoton.minus())/(self.pLaser.minus())
        self.alphas = np.array([alpha(self.bwMomenta,self.config)(item) for item in [1,2,3]])



        self.myVertex,self.myMom,self.mySp,self.myPhotoNum = self.vertexObj.testGetVertex(self.bwMomenta,self.spinors,self.alphas,self.photoNum,verba=True)
        self.myMatrixEl0 = np.array([self.myVertex[item]*self.polPhoto[0] for item in np.arange(4)])/(self.pLaser.minus()*0.5)
        self.myMatrixEl1 = np.array([self.myVertex[item]*self.polPhoto[1] for item in np.arange(4)])/(self.pLaser.minus()*0.5)
        self.amp = 0.5*(np.dot(self.myMatrixEl0,np.conjugate(self.myMatrixEl0)) + np.dot(self.myMatrixEl1,np.conjugate(self.myMatrixEl1)))
        end=time.time()-start
        print "my time: %1.2e"%end
        self.assertLessEqual(np.abs(self.amp-self.MatTobRes),self.eps)




if __name__=='__main__':
    unittest.main()
