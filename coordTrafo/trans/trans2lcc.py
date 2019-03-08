"""
transverse to light cone coordinates
"""
import numpy as np

def pm(y,pt,phi,m=1.0):
    mt = np.sqrt(pt**2 + m**2)
    return mt*np.exp(-y)/2.0

def px(y,pt,phi,m=1.0):
    return pt*np.cos(phi)

def py(y,pt,phi,m=1.0):
    return pt*np.sin(phi)

def y(pm,px,py,m=1.0):
    #print "-"*20
    #print "pm: %s"%pm
    #print "px: %s"%px
    #print "py: %s"%py
    pperp = px**2 + py**2
    return .5*np.log((pperp + m**2)/(4*(pm**2)))

def pt(pm,px,py,m=1.0):
    return np.sqrt(px**2 + py**2)

def phi(pm,px,py,m=1.0):
    return np.arctan(py/px)



if __name__=='__main__':

    def trans2LCCAndBack(Y,PT,PHI):
        tempPM = pm(Y,PT,PHI)
        tempPX = px(Y,PT,PHI)
        tempPY = py(Y,PT,PHI)
        #print "-"*20
        #print tempE
        #print tempC
        #print Y
        #print PT
        tempY = y(tempPM,tempPX,tempPY)
        tempPT = pt(tempPM,tempPX,tempPY)
        tempPHI = phi(tempPM,tempPX,tempPY)
        return tempY,tempPT,tempPHI

    testY = np.linspace(-4.0,4.0,20)
    testPT = np.linspace(1e-4,3.0,20)
    testPHI = np.linspace(-np.pi/2.0,np.pi/2.0,10)
    eps=1e-6
    #print c(0.6315789473684212,1e-7)
    #print np.sqrt(1e-8**2 +1)
    #print np.arccos(1.0000000000000002)
    for elY in testY:
        for elPT in testPT:
            for elPHI in testPHI:
                yb,ptb,phib = trans2LCCAndBack(elY,elPT,elPHI)
                if (np.abs(elY-yb)>eps) or (np.abs(elPT-ptb)>eps) or (np.abs(elPHI-phib)>eps):

                    print "y:    %s"%elY
                    print "pt:   %s"%elPT
                    print "phi:  %s"%elPHI
                    print "yb:   %s"%yb
                    print "ptb:  %s"%ptb
                    print "phib: %s"%phib

    #print c(1.4,2.033333)
