"""
contains the current handling

todo:
    -   think about inheritance: default class with the stuff from cleanSet
    -   check sign in J1 and J2 (polariation!)

spinor order: 00,11,01,10

"""
import numpy as np
import qft
from kinutility import laserPolarisation
gamma=qft.GammaMatrix()
#print "is there feyndagg: %s"%('feyndagg' in dir(qft))





class currentClass(object):
    def __init__(self,config=None):
        self.a0=config['a0']
        self.mass = config['mass']
        self.xi=config['xi']

        self.polBase =[laserPolarisation(item) for item in [1,2]]#-> right choice
        #self.polBase =[ laserPolarisation(1)*np.cos(self.xi)-1j*laserPolarisation(2)*np.sin(self.xi),laserPolarisation(1)*np.cos(self.xi)+1j*laserPolarisation(2)*np.sin(self.xi)]
        self.cleanSet()

    def cleanSet(self):
        self.J0 =None# np.array([qft.MinkowskiVector([0,0,0,0]) for t in range(4)])
        self.J1 =None# np.array([qft.MinkowskiVector([0,0,0,0]) for t in range(4)])
        self.J2 =None# np.array([qft.MinkowskiVector([0,0,0,0]) for t in range(4)])
        self.J3 = None#np.array([qft.MinkowskiVector([0,0,0,0]) for t in range(4)])
        self.momPhoto = None#qft.MinkowskiVector([0,0,0,0])
        self.momA = None #qft.MinkowskiVector([0,0,0,0])
        self.momB = None #qft.MinkowskiVector([0,0,0,0])
        self.spinorA = None
        self.spinorB = None


    def generalJ0(self,spinorA,spinorB):
        """
        calc of J0 for general spinors
        """
        current1 = np.zeros(4,dtype=np.complex)
        current2 = np.zeros(4,dtype=np.complex)
        current3 = np.zeros(4,dtype=np.complex)
        current4 = np.zeros(4,dtype=np.complex)
        for item in np.arange(4):
            diracMatrix = gamma[item]
            current1[item] = (spinorA[0]*(diracMatrix*spinorB[0]))
            current2[item] = (spinorA[1]*(diracMatrix*spinorB[1]))
            current3[item] = (spinorA[0]*(diracMatrix*spinorB[1]))
            current4[item] = (spinorA[1]*(diracMatrix*spinorB[0]))
        return np.array([qft.MinkowskiVector(current1),qft.MinkowskiVector(current2),qft.MinkowskiVector(current3),qft.MinkowskiVector(current4)])

    def generalJ12(self,momPhoto,momA,momB,spinorA,spinorB,pol):
        """
        calc of J1/2 for general spinors and momenta

        pol ... mks.Vector for background (general) polarisation
        """
        current1 = np.zeros(4,dtype=np.complex)
        current2 = np.zeros(4,dtype=np.complex)
        current3 = np.zeros(4,dtype=np.complex)
        current4 = np.zeros(4,dtype=np.complex)
        factor1 = self.a0*self.mass*qft.feyndagg(pol)*qft.feyndagg(momPhoto)/(2*momPhoto*momA)
        factor2 = self.a0*self.mass*qft.feyndagg(momPhoto)*qft.feyndagg(pol)/(2*momPhoto*momB)
        for item in np.arange(4):
            diracMatrix1 = factor1*gamma[item]
            diracMatrix2= gamma[item]*factor2
            current1[item] = spinorA[0]*(diracMatrix1*spinorB[0]) - spinorA[0]*(diracMatrix2*spinorB[0])
            current2[item] =spinorA[1]*(diracMatrix1*spinorB[1]) - spinorA[1]*(diracMatrix2*spinorB[1])
            current3[item] =spinorA[0]*(diracMatrix1*spinorB[1]) - spinorA[0]*(diracMatrix2*spinorB[1])
            current4[item] =spinorA[1]*(diracMatrix1*spinorB[0]) - spinorA[1]*(diracMatrix2*spinorB[0])
        return np.array([qft.MinkowskiVector(current1),qft.MinkowskiVector(current2),qft.MinkowskiVector(current3),qft.MinkowskiVector(current4)])

    def generalJ3(self,momPhoto,momA,momB,spinorA,spinorB):
        """
        calc of J3 for general spinors and momenta
        """
        prefac=-self.a0**2*self.mass**2/2.0
        diracMatrix = prefac * qft.feyndagg(momPhoto)/(momPhoto*momB)/(momPhoto*momA)
        current1 = spinorA[0]*(diracMatrix*spinorB[0])
        current2 = spinorA[1]*(diracMatrix*spinorB[1])
        current3 = spinorA[0]*(diracMatrix*spinorB[1])
        current4 = spinorA[1]*(diracMatrix*spinorB[0])
        return np.array([momPhoto*current1,momPhoto*current2,momPhoto*current3,momPhoto*current4])

    def setKin(self,momPhoto,momA,momB,spinorA,spinorB):
        self.momPhoto = momPhoto
        self.momA = momA
        self.momB = momB
        self.spinorA = spinorA
        self.spinorB = spinorB
        self.evalJ0()
        self.evalJ1()
        self.evalJ2()
        self.evalJ3()

    def evalJ0(self):
        self.J0 = self.generalJ0(self.spinorA,self.spinorB)
        return self.J0

    def evalJ1(self):
        self.J1 = self.generalJ12(self.momPhoto,self.momA,self.momB,self.spinorA,self.spinorB,self.polBase[0])
        return self.J1

    def evalJ2(self):
        self.J2 = self.generalJ12(self.momPhoto,self.momA,self.momB,self.spinorA,self.spinorB,self.polBase[1])
        return self.J2

    def evalJ3(self):
        self.J3 = self.generalJ3(self.momPhoto,self.momA,self.momB,self.spinorA,self.spinorB)
        return self.J3
