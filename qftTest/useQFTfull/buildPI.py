
# coding: utf-8

# Build Phase integrals on a grid
# ===================

# In[ ]:


import sftrident as sf
print "sftrident version: %s"%(str(sf.__version__))
import sftrident.qft as qft
import time
import numpy as np
import sftrident.pulseLib as pulseLib
from sftrident.mathUtil import isOdd
from filon import cos_integral,sin_integral

times = {}


# In[ ]:




pol1=sf.kinutility.BGpolarisationBase(1)
pol2=sf.kinutility.BGpolarisationBase(2)

print pol1
print pol2

class alpha(object):
    r"""Calculation of kinematic factors.

    Methodes
    --------
    __call__(item)
        represents the kinemtatic functions alpha

    Notes
    -----
    Calculation of

    .. math:: \alpha _{1,2} = 2ma_0\left(\frac{p_a \epsilon_{1,2}}{2 k p_a} - \frac{p_b \epsilon_{1,2}}{2 k p_b}\right)

    and

    .. math:: \alpha_3 = m^2a_0^2\left(\frac{1}{2kp_a} + \frac{1}{2kp_b}\right).

    Signs of $p_a, p_b$ is chosen from Breit-Wheeler-Part.
    """
    def __init__(self,mom,config):
        r""" Initializes the alpha class.

        Parameters
        ----------
        mom : array_like
            A 3-array which contains the momenta as qft.minkowskiVector.
            The format is mom = [k,pa,pb].
        config : dict
            configuration dictionary which contains at least a0 and mass.
        """
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

        Raises
        ------
        ValueError
            <item> not in (1,2,3)
        """
        if item==1:
            return 2.0*((self.__pa*pol1)/(2.0*self.__pak) - (self.__pb*pol1)/(2.0*self.__pbk))*self.__a0*self.__me
        elif item==2:
            return 2.0*((self.__pa*pol2)/(2.0*self.__pak) - (self.__pb*pol2)/(2.0*self.__pbk))*self.__a0*self.__me
        elif item==3:
            #print self.__pak
            #print self.__pbk
            return self.__me**2*self.__a0**2*(1.0/(2.0*self.__pak) + 1.0/(2.0*self.__pbk))
        else:
            raise ValueError("The index for alpha needs to be 1, 2 or 3! (<%s> given)"%item)


# In[ ]:


def BuildGrid(ss,p1x,p1y,p1m,p2x,p2y,p2m):
    SS, P1x, P1y, P1m, P2x, P2y, P2m = np.meshgrid(ss,p1x,p1y,p1m,p2x,p2y,p2m,indexing='ij')
    return [SS,P1x,P1y,P1m,P2x,P2y,P2m]

ssInit = np.array([3.353])
#p1_x = np.linspace(0.0,3.5,3)
p1_x = np.random.uniform(0.0,3.5,3)
#print "p1_x: %s"%p1_x
#p1_y = np.linspace(0.0,3.5,3)
p1_y = np.random.uniform(0.0,3.5,3)
#print "p1_y: %s"%p1_y
#p1_m = np.linspace(0.0,3.5,20)
p1_m = np.random.uniform(0.0,3.5,40)
#print "p1_m: %s"%p1_m

#p2_x = np.linspace(0.0,3.5,3)
p2_x = np.random.uniform(0.0,3.5,3)
#print "p2_x: %s"%p2_x
#p2_y = np.linspace(0.0,3.5,3)
p2_y = np.random.uniform(0.0,3.5,3)
#print "p2_y: %s"%p2_y
#p2_m = np.linspace(-1.0,3.5,20)
p2_m = np.random.uniform(-1.0,3.5,20)
#print "p2_m: %s"%p2_m

start = time.time()
kinGridRaw = BuildGrid(ssInit,p1_x,p1_y,p1_m,p2_x,p2_y,p2_m)
end = time.time() - start
times['mesh'] = end
for i,el in enumerate(kinGridRaw):
    print "coord %d: %s"%(i,str(el.shape))

def omegaFromSS(ss):
    #lab
    return (ss**2 - 1.0)/2.0

#initial particles

start = time.time()
omega = omegaFromSS(kinGridRaw[0])
end =  time.time() -start
times['om']=end
print "omega shape: %s"%(str(omega.shape))
#print "omega: %s"%omega

start =  time.time()
E = np.ones(omega.shape)
end =  time.time() - start
times['E'] = end
print "E shape: %s"%(str(E.shape))
#print "E: %s"%E

start =  time.time()
#P=qft.MinkowskiVector([E,np.zeros(E.shape),np.zeros(E.shape),np.zeros(E.shape)])
P=qft.MinkowskiVector([E,0.0,0.0,0.0])
end =  time.time() - start
times['Pfirst']=end

start =  time.time()
P3m = 0.5*(P._0() - P._3()) - kinGridRaw[3] - kinGridRaw[6]
P3x =  - kinGridRaw[1] - kinGridRaw[4]
P3y =  - kinGridRaw[2] - kinGridRaw[5]
end =  time.time() - start
times['p3coord'] = end

print "P3m shape: %s"%(str(P3m.shape))
print "P3x shape: %s"%(str(P3x.shape))
print "P3y shape: %s"%(str(P3y.shape))

def physArea(P1m,P2m,P3m):
    return (P1m>0) * (P2m>0) * (P3m>0)


def BuildFinalMom(px,py,pm):
    #print pm
    pp = (px**2 + py**2 + 1)/(4.0*pm)
    return qft.MinkowskiVector(qft.parray([pp+pm,px,py,pp-pm]))

def buildAllMom(om,e,p1x,p1y,p1m,p2x,p2y,p2m,p3x,p3y,p3m):
    physAreaArr = physArea(p1m,p2m,p3m)
    #print physAreaArr.all()==False
    p1xT = p1x[physAreaArr]
    p1yT = p1y[physAreaArr]
    p1mT = p1m[physAreaArr]
    #print "p1m: %s"%(str(p1m.shape))
    #print "p1mT: %s"%(str(p1mT.shape))
    p2xT = p2x[physAreaArr]
    p2yT = p2y[physAreaArr]
    p2mT = p2m[physAreaArr]
    #print "p2m: %s"%(str(p2m.shape))
    #print "p2mT: %s"%(str(p2mT.shape))
    p3xT = p3x[physAreaArr]
    p3yT = p3y[physAreaArr]
    p3mT = p3m[physAreaArr]
    #print "p3m: %s"%(str(p3m.shape))
    #print "p3mT: %s"%(str(p3mT.shape))
    P1 = BuildFinalMom(p1xT,p1yT,p1mT)
    P2 = BuildFinalMom(p2xT,p2yT,p2mT)
    P3 = BuildFinalMom(p3xT,p3yT,p3mT)

    omegaT = om[physAreaArr]
    ET = e[physAreaArr]
    #P=qft.MinkowskiVector([ET,np.zeros(ET.shape),np.zeros(ET.shape),np.zeros(ET.shape)])
    P=qft.MinkowskiVector([ET,0.0,0.0,0.0])
    K = qft.MinkowskiVector([omegaT,np.zeros(omegaT.shape),np.zeros(omegaT.shape),omegaT])

    return P1,P2,P3,K,P

start =  time.time()
Pa,Pb,Pc,K,P = buildAllMom(omega,E,kinGridRaw[1],kinGridRaw[2],kinGridRaw[3],kinGridRaw[4],kinGridRaw[5],kinGridRaw[6],P3x,P3y,P3m)
end =  time.time() - start
times['allMom'] = end


print "Pa: %s"%(str(Pa.shape))
print "Pb: %s"%(str(Pb.shape))
print "Pc: %s"%(str(Pc.shape))
print "P: %s"%(str(P.shape))
print "K: %s"%(str(K.shape))

#print Pa._1()


# In[ ]:


tempDPHI = 50.0
config = {'a0':0.01,'mass':1.0,'xi':0.0,'dPhi':tempDPHI,'psBounds':[-tempDPHI,tempDPHI],'envelope':'cos','pulseOpt':['analytic'],'deg':2500,'mode':'filon'}
alphaFKT = alpha([K,Pa,Pb],config)
alphas = np.array([alphaFKT(i) for i in (1,2,3)])
#print alphas[2]


# In[ ]:


class phaseIntegralFilon(object):
    def __init__(self,config = {'xi':0.0,'dPhi':10,'psBounds':[-10,10],'envelope':'cos','pulseOpt':[]},deg=5):
        self.__xi = config['xi']
        self.__dPhi = config['dPhi']
        self.__pulse = pulseLib.getPulse(config['envelope'],config['pulseOpt'])
        self.__envelope = self.__pulse[0]
        self.__internalInt = self.__pulse[1]
        if isOdd(deg):
            self.__deg = deg
        else:
            self.__deg = deg + 1

        tempPoints, self.__steps = np.linspace(-1,1,self.__deg,retstep=True,endpoint=True)
        self.__points = tempPoints[:,np.newaxis]
        self.__evalPoints(self.__points)

    def __evalPoints(self,points):
        print "shape points: %s"%(str(self.__points.shape))
        self.__envGeneral = self.__envelope(self.__points*self.__dPhi,self.__dPhi)
        print "shape env: %s"%(str(self.__envGeneral.shape))
        self.__env1 = self.__envGeneral
        self.__beta1 = self.__evalBeta1(points*self.__dPhi)
        print "shape beta1:  %s"%(str(self.__beta1.shape))
        self.__beta2 = self.__evalBeta2(points*self.__dPhi)
        self.__beta3 = self.__evalBeta3(points*self.__dPhi)

    def setKin(self,alphas):
        self.__alphas = [alphas[i][np.newaxis,:] for i in np.arange(3)]
        print "shape alpha: %s"%(str(self.__alphas[0].shape))
        print "shape beta:  %s"%(str(self.__beta1.shape))
        test1 = 1j*self.__alphas[0]*self.__beta1
        print "t1: %s"%(str(test1.shape))
        test2 = 1j*self.__alphas[1]*self.__beta2
        #print "t2: %s"%(str(test2.shape))
        test3 = 1j*self.__alphas[2]*self.__beta3
        #print "t3: %s"%(str(test3.shape))

        self.__expTerm = np.exp(1j*self.__alphas[0]*self.__beta1 + 1j*self.__alphas[1]*self.__beta2 + 1j*self.__alphas[2]*self.__beta3)
        self.__integrand1 = self.__env1*self.__expTerm
        return 0

    def setPhotoNum(self,photoNum):
        """
        NEEDS self.setKin first --> to avoid call issues: self.__setclean
        """
        #print "pnum pre: %s"%photoNum
        photoNum = np.require(photoNum)
        if photoNum.ndim==0:
            photoNum = np.array([photoNum])
        self.__photoNum = photoNum

        #I_1(s+1)
        self.__sArgs11 = self.__dPhi*(self.__photoNum + np.ones(self.__photoNum.shape))
        resC11 = cos_integral(self.__integrand1,self.__steps,self.__sArgs11,-1,axis=0)
        resS11 = sin_integral(self.__integrand1,self.__steps,self.__sArgs11,-1,axis=0)
        self.__resFilon11 = resC11 + 1j*resS11
        #print "pnum: %s"%photoNum
        #print "resFilon: %s"%np.complex64(self.__resFilon11)

        #I_1(s-1)
        self.__sArgs12 = self.__dPhi*(self.__photoNum - np.ones(self.__photoNum.shape))
        resC12 = cos_integral(self.__integrand1,self.__steps,self.__sArgs12,-1,axis=0)
        resS12 = sin_integral(self.__integrand1,self.__steps,self.__sArgs12,-1,axis=0)
        self.__resFilon12 = resC12 + 1j*resS12



    def __evalBeta1(self,phiArray):
        """
        returns the evaluation of beta1 at the given gauss points
        """
        print "shape phiArr: %s"%(str(phiArray.shape))
        print "shape dphi: %s"%(str(self.__dPhi))
        res = np.zeros(phiArray.shape)
        if np.cos(self.__xi)==0:
            return res
        else:
            print "internalInt: %s"%(str(self.__internalInt(0)(phiArray,self.__dPhi).shape))
            #print self.__internalInt(0)(phiArray,self.__dPhi)
            #print self.__dPhi
            return np.cos(self.__xi)*(self.__internalInt(0)(phiArray,self.__dPhi))

    def __evalBeta2(self,phiArray):
        """
        returns the evaluation of beta2 at the given gauss points
        """
        res = np.zeros(self.__points.shape)
        if np.sin(self.__xi)==0:
            return res
        else:
            return np.sin(self.__xi)*(self.__internalInt(1)(phiArray,self.__dPhi))

    def __evalBeta3(self,phiArray):
        """
        returns the evaluation of beta3 at the given gauss points
        """
        res = np.zeros(self.__points.shape)
        if np.sin(self.__xi)==0:
            #print "internalInt: %s"%(str(self.__internalInt(2)(phiArray,self.__dPhi).shape))
            #print self.__internalInt(0)(phiArray,self.__dPhi)
            #print self.__dPhi
            return (np.cos(self.__xi))**2*(self.__internalInt(2)(phiArray,self.__dPhi))
        elif np.cos(self.__xi)==0:
            return (np.sin(self.__xi))**2*(self.__internalInt(3)(phiArray,self.__dPhi))
        else:
            return (np.cos(self.__xi))**2*(self.__internalInt(2)(phiArray,self.__dPhi)) + (np.sin(self.__xi))**2*(self.__internalInt(3)(phiArray,self.__dPhi))


    def getB1(self):
        res1 = np.cos(self.__xi)*np.complex64(self.__dPhi/2.0*(self.__resFilon11 + self.__resFilon12))
        if res1.shape[-1]==1:
            self.B1 = np.asscalar(res1)
        else:
            self.B1 = res1
        return self.B1



    def B1func(self,photoNum):
        if np.isscalar(photoNum):
            photoNum = np.array([photoNum])

        #I_1(s+1)
        sArgs11 = self.__dPhi*(photoNum + np.ones(photoNum.shape))
        resC11 = cos_integral(self.__integrand1,self.__steps,sArgs11,-1,axis=0)
        resS11 = sin_integral(self.__integrand1,self.__steps,sArgs11,-1,axis=0)
        resFilon11 = resC11 + 1j*resS11

        #I_1(s-1)
        sArgs12 = self.__dPhi*(photoNum - np.ones(photoNum.shape))
        resC12 = cos_integral(self.__integrand1,self.__steps,sArgs12,-1,axis=0)
        resS12 = sin_integral(self.__integrand1,self.__steps,sArgs12,-1,axis=0)
        resFilon12 = resC12 + 1j*resS12
        #print "integrand shape: %s"%(str(self.__integrand1.shape))
        #print "sargs shape: %s"%(str(sArgs12.shape))
        #print "res shape: %s"%(str(resFilon12.shape))
        return np.cos(self.__xi)*np.complex64(self.__dPhi/2.0*(resFilon11 + resFilon12))



# In[ ]:


start = time.time()
testPI = phaseIntegralFilon(config,config['deg'])
end = time.time() - start
times['PIinit'] = end


# In[ ]:



start = time.time()
testPI.setKin(alphas)
end = time.time() - start
times['PIkin'] = end


# In[ ]:


sTest = np.linspace(0.1,5.0,243)

start = time.time()
testPI.setPhotoNum(sTest)
end = time.time() - start
times['PIpnum'] = end




start = time.time()
testPIvals = testPI.getB1()
end = time.time() - start
times['PIvals'] = end
print testPIvals.shape

#print testPIvals


# In[ ]:


Npoints =K.shape[0]
print "eval points: %d (x %d = %d)"%(Npoints,len(sTest),Npoints*len(sTest))
print '-'*25
print "TIMES:"
print '-'*25
for el in times.iterkeys():
    print "%s: %1.4e"%(el,times[el])
    print "%s: %1.4e (avg.)"%(el,times[el]/Npoints)
    print ''
print ""
print '-'*25
print "TIMES (cummulative)"
print '-'*25
fullTime = 0.0
for el in times.iterkeys():
    fullTime+=times[el]
    print "%s: %1.4e"%(el,fullTime)
    print "%s: %1.4e (avg.)"%(el,fullTime/Npoints)
    print ''
print '-'*25
print "TIMES (properties)"
print '-'*25
maxVal = ['',0.0]

for el in times.iterkeys():
    if times[el]>=maxVal[1]:
        maxVal = [el,times[el]]
print "max: %s (%1.4e)"%(tuple(maxVal))
minVal = maxVal
for el in times.iterkeys():
    if times[el]<=minVal[1]:
        minVal = [el,times[el]]
print "min: %s (%1.4e)"%(tuple(minVal))



# Compare with serial evaluation
# -------------

# In[ ]:


start = time.time()
SQtestPI = sf.phaseIntegral(config,config['deg'])
end = time.time() - start

maxRelErrR = [0,0,0.0]
maxAbsErrR = [0,0,0.0]
maxRelErrI = [0,0,0.0]
maxAbsErrI = [0,0,0.0]

SQtimes = end


for i, als in enumerate([alphas[:,i] for i in np.arange(len(alphas[0]))]):
    #print '-'*25
    #print "%d: %s"%(i,als)
    start = time.time()
    SQtestPI.setKin(als)
    end = time.time() - start
    SQtimes += end
    for indS,sArg in enumerate(sTest):
        start = time.time()
        SQpiVal = SQtestPI.B1func(sArg)[0]
        end = time.time() - start
        SQtimes += end
        piVal = testPIvals[indS,i]

        #error estimation
        absErrR = np.abs(piVal.real - SQpiVal.real)
        #print absErrR
        if absErrR>=maxAbsErrR[-1]:
            maxAbsErrR = [i,indS,absErrR]

        relErrR = absErrR/(np.abs(piVal.real + SQpiVal.real))
        if relErrR>=maxRelErrR[-1]:
            maxRelErrR = [i,indS,absErrR,relErrR]

        absErrI = np.abs(piVal.imag - SQpiVal.imag)
        if absErrI>=maxAbsErrI[-1]:
            maxAbsErrI = [i,indS,absErrI]

        relErrI = absErrI/(np.abs(piVal.imag + SQpiVal.imag))
        if relErrI>=maxRelErrI[-1]:
            maxRelErrI = [i,indS,absErrI,relErrI]


print "----------- ERROR -----------"
print "max abs. real err.: %s (alpha: %s , pNum: %s)"%(maxAbsErrR[-1],alphas[:,maxAbsErrR[0]],sTest[maxAbsErrR[1]])
print "max rel. real err.: %s (alpha: %s , pNum: %s, absErr: %1.4e)"%(maxRelErrR[-1],alphas[:,maxRelErrR[0]],sTest[maxRelErrR[1]],maxRelErrR[2])

print "max abs. imag err.: %s (alpha: %s , pNum: %s)"%(maxAbsErrI[-1],alphas[:,maxAbsErrI[0]],sTest[maxAbsErrI[1]])
print "max rel. imag err.: %s (alpha: %s , pNum: %s, absErr: %1.4e)"%(maxRelErrI[-1],alphas[:,maxRelErrI[0]],sTest[maxRelErrI[1]],maxRelErrI[2])

print "----------- TIMES -----------"
print "full time (ser.): %1.4e"%SQtimes
print "avg. time (ser.): %1.4e"%(SQtimes/(float(len(sTest)*len(alphas[0]))))
PIwords = ['PIinit','PIkin','PIpnum','PIvals']
piTime = np.sum([times[el] for el in PIwords])

print "full time (mesh): %1.4e"%piTime
print "avg. time (mesh): %1.4e"%(piTime/(float(len(sTest)*len(alphas[0]))))



# In[ ]:


print "%1.4e"%(1.1602e-02/70.0)


# In[ ]:


3e-4*5e5*50.0*50.0/60.0/60.0/24.0


# In[ ]:


int(1e5**(1/4.0))
