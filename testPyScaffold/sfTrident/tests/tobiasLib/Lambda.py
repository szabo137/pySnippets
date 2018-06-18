from Config_paar import *
from A_Funktion import *
from Kinematik import *
import numpy
from zerhacker import *
from NME import *
from asimps import *

m=1.0

def fd(x):
    return feyndagg(x)

def Lambda_2():
    p_pos,p_el,k_laser,k_photon,q_pos,eps_m,eps_p = kinematik()


    r=2
    rp=2
    lp=2
    V     = array([SpinorV((p_pos, m),s1) for s1 in [1,2]])


    Ubar  = array([SpinorUBar((p_el,m),s1) for s1 in [1,2]])
    #print "ubar: %s"%(Ubar)
    #print "test ubar: %s"%(SpinorUBar((p_el,m),r))
    #print "ubar(%s): %s"%(r,Ubar[r-1])
    #print "v(%s): %s"%(rp,V[rp-1])
    #print "eps(%s): %s"%(lp,eps_photon(lp))

    d_el  = m*a0/(4.*k_laser*p_el)
    d_pos = m*a0/(4.*k_laser*p_pos)

    J_00  = lambda l: fd(eps_photon(l))
    #print "j0: %s"%(Ubar[r-1]*J_00(lp)*V[rp-1])
    J_11  = lambda l: d_el * fd(eps_m)* fd(k_laser)   * fd(eps_photon(l)) - d_pos * fd(eps_photon(l))*fd(k_laser)* fd(eps_m)
    J_1_1 = lambda l: d_el * fd(eps_p)* fd(k_laser)   * fd(eps_photon(l)) - d_pos * fd(eps_photon(l))*fd(k_laser)* fd(eps_p)
    J_20  = lambda l: -4.  *  k_laser * eps_photon(l) * d_pos * d_el  * fd(k_laser)
    J_22  = lambda l: 0.5  *  J_20(l) * (cos(ksi)**2  - sin(ksi)**2)
    J_2_2 = lambda l:         J_22(l)

    if Envelope == 'cos^2':
        a = sigma/w_laser
    elif Envelope == 'Gauss':
        a = sigma*5.
    elif Envelope == 'Box':
        a = sigma+1


    #@zerhacker( slicearguments = (0,1) , slicelength = 250 , axis = 1)
    #@zerhacker( slicearguments = (0,1) , slicelength = 250 , axis = 0)
    def integrate(p_el,p_pos,k_photon,k_laser,M,N):
        return asimps( lambda x: A_m_n_nSVEA(M,N,x,p_el,p_pos,k_photon,k_laser),- a, a, Nx=421 , errorabs = 1e-4, maxrecur = 100 )


    A_11  = integrate(p_el,p_pos,k_photon,k_laser,1,1)
    #print "A11: %s"%(str(A_11))
    A_1_1 = integrate(p_el,p_pos,k_photon,k_laser,1,-1)
    #print "A1-1: %s"%(str(A_1_1))
    A_20  = integrate(p_el,p_pos,k_photon,k_laser,2,0)
    #print "A20: %s"%(str(A_20))

    A_22  = integrate(p_el,p_pos,k_photon,k_laser,2,2)
    #print "A22: %s"%(str(A_22))
    A_2_2  = integrate(p_el,p_pos,k_photon,k_laser,2,-2)

    A_00  = A_0_0(A_11[0],A_1_1[0],A_20[0],A_22[0],A_2_2[0] )
    #print "A00: %s"%(str(A_00))


    Lambda_1  = array([ J_00(l)*parray(A_00) + ( J_11(l)*parray(A_11[0])+J_1_1(l)*parray(A_1_1[0])+\
                J_20(l)*parray(A_20[0]) ) for l in [1,2]])

    #print "lambda: %s"%(Lambda_1)
    #print "lambda info: %s"%(Lambda_1[newaxis,:,newaxis])

    Lambda    = (Ubar[newaxis,newaxis,:]*Lambda_1[newaxis,:,newaxis]*V[:,newaxis,newaxis])


    Lambda2  = 0.5 * sum(abs(Lambda)**2)

    #print "pPos:: %s"%(p_pos)
    return (Lambda2,p_pos,p_el,k_laser,k_photon,q_pos,eps_m,eps_p)

def LambdaSingle():
    p_pos,p_el,k_laser,k_photon,q_pos,eps_m,eps_p = kinematik()


    r=2
    rp=2
    lp=2
    V     = array([SpinorV((p_pos, m),s1) for s1 in [1,2]])


    Ubar  = array([SpinorUBar((p_el,m),s1) for s1 in [1,2]])
    #print "ubar: %s"%(Ubar)
    #print "test ubar: %s"%(SpinorUBar((p_el,m),r))
    #print "ubar(%s): %s"%(r,Ubar[r-1])
    #print "v(%s): %s"%(rp,V[rp-1])
    #print "eps(%s): %s"%(lp,eps_photon(lp))

    d_el  = m*a0/(4.*k_laser*p_el)
    d_pos = m*a0/(4.*k_laser*p_pos)

    J_00  = lambda l: fd(eps_photon(l))
    #print "j0: %s"%(Ubar[r-1]*J_00(lp)*V[rp-1])
    J_11  = lambda l: d_el * fd(eps_m)* fd(k_laser)   * fd(eps_photon(l)) - d_pos * fd(eps_photon(l))*fd(k_laser)* fd(eps_m)
    J_1_1 = lambda l: d_el * fd(eps_p)* fd(k_laser)   * fd(eps_photon(l)) - d_pos * fd(eps_photon(l))*fd(k_laser)* fd(eps_p)
    J_20  = lambda l: -4.  *  k_laser * eps_photon(l) * d_pos * d_el  * fd(k_laser)
    J_22  = lambda l: 0.5  *  J_20(l) * (cos(ksi)**2  - sin(ksi)**2)
    J_2_2 = lambda l:         J_22(l)

    print "-----Tobias Currents ------"
    testJ00dir = np.array([J_00(item) for item in [1,2]])
    testJ00=(Ubar[newaxis,newaxis,:]*testJ00dir[newaxis,:,newaxis]*V[:,newaxis,newaxis])
    for item in [0,1]:
        for item2 in [0,1]:
            print "J00 (%s,%s): %s"%(item,item2,testJ00[item][0][item2])

    testJ11dir = np.array([J_11(item) for item in [1,2]])
    testJ11=(Ubar[newaxis,newaxis,:]*testJ11dir[newaxis,:,newaxis]*V[:,newaxis,newaxis])
    for item in [0,1]:
        for item2 in [0,1]:
            print "J11 (%s,%s): %s"%(item,item2,testJ11[item][0][item2])

    testJ1_1dir = np.array([J_1_1(item) for item in [1,2]])
    testJ1_1=(Ubar[newaxis,newaxis,:]*testJ1_1dir[newaxis,:,newaxis]*V[:,newaxis,newaxis])
    for item in [0,1]:
        for item2 in [0,1]:
            print "J1_1 (%s,%s): %s"%(item,item2,testJ1_1[item][0][item2])

    testJ20dir = np.array([J_20(item) for item in [1,2]])
    testJ20=(Ubar[newaxis,newaxis,:]*testJ20dir[newaxis,:,newaxis]*V[:,newaxis,newaxis])
    for item in [0,1]:
        for item2 in [0,1]:
            print "J20 (%s,%s): %s"%(item,item2,testJ20[item][0][item2])


    if Envelope == 'cos^2':
        a = sigma/w_laser
    elif Envelope == 'Gauss':
        a = sigma*5.
    elif Envelope == 'Box':
        a = sigma+1


    #@zerhacker( slicearguments = (0,1) , slicelength = 250 , axis = 1)
    #@zerhacker( slicearguments = (0,1) , slicelength = 250 , axis = 0)
    def integrate(p_el,p_pos,k_photon,k_laser,M,N):
        return asimps( lambda x: A_m_n_nSVEA(M,N,x,p_el,p_pos,k_photon,k_laser),- a, a, Nx=421 , errorabs = 1e-4, maxrecur = 100 )


    A_11  = integrate(p_el,p_pos,k_photon,k_laser,1,1)
    #print "A11: %s"%(str(A_11))
    A_1_1 = integrate(p_el,p_pos,k_photon,k_laser,1,-1)
    #print "A1-1: %s"%(str(A_1_1))
    A_20  = integrate(p_el,p_pos,k_photon,k_laser,2,0)
    #print "A20: %s"%(str(A_20))
    A_00  = A_0_0(A_11[0],A_1_1[0],A_20[0],A_11[0]*0,A_11[0]*0)
    #print "A00: %s"%(str(A_00))

    #not contained in original?!
    A_22  = integrate(p_el,p_pos,k_photon,k_laser,2,2)
    #print "A22: %s"%(str(A_22))
    A_2_2  = integrate(p_el,p_pos,k_photon,k_laser,2,-2)



    print "------------ tobias phase integral --------------"
    pKlaser = k_laser
    B1 = cos(ksi)/2.0*pKlaser.minus()/2.0*(A_11[0] + A_1_1[0])
    B2 = np.sin(ksi)/(2.0*1j)*pKlaser.minus()/2.0*(A_11[0] - A_1_1[0])
    B3 = pKlaser.minus()/4.0*A_20[0] + pKlaser.minus()/8.0*np.cos(2.0*ksi)*(A_22[0] + A_2_2[0])
    B0 = pKlaser.minus()/2.0*A_00

    print "B0: %s"%(B0)
    print "B1: %s"%(B1)
    print "B2: %s"%(B2)
    print "B3: %s"%(B3)


    Lambda_1  = array([ J_00(l)*parray(A_00) + ( J_11(l)*parray(A_11[0])+J_1_1(l)*parray(A_1_1[0])+\
                J_20(l)*parray(A_20[0]) ) for l in [1,2]])

    #print "lambda: %s"%(Lambda_1)
    #print "lambda info: %s"%(Lambda_1[newaxis,:,newaxis])

    Lambda    = (Ubar[newaxis,newaxis,:]*Lambda_1[newaxis,:,newaxis]*V[:,newaxis,newaxis])


    Lambda2  = 0.5 * sum(abs(Lambda)**2)

    return (Lambda,p_pos,p_el,k_laser,k_photon,q_pos,eps_m,eps_p)


if __name__ == '__main__':
    import time
    t1=time.time()
    res =  Lambda_2()
    t2=time.time()
    print "res: %s (time: %s)"%(res[0],t2-t1)
    print "pPos: %s"%res[1]
    print "pEl: %s"%res[2]
    print "kLaser: %s"%res[3]
    print "kPhoton: %s"%res[4]
    print "qPos: %s"%res[5]
    print "epsM: %s"%res[6]
    print "epsP: %s"%res[7]
