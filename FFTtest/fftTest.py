"""
module to test FFT for fourier transforms
"""
import numpy as np


def func(x):
    return np.exp(-(x+5)**2)
    
def func2(x,N):
    xOLD = x*N/(2.0*np.pi)
    return func(x)*N/(2.0*np.pi)


def exact(w):
    return np.exp(-w**2/4.0-1j*5*w)/(np.sqrt(2))



import numpy as np
import matplotlib.pyplot as pl

#Consider function f(t)=1/(t^2+1)
#We want to compute the Fourier transform g(w)

#Discretize time t
t0=-100.
dt=0.001
t=np.arange(t0,-t0,dt)
#Define function
f=func(t)

#Compute Fourier transform by numpy's FFT function
g=np.fft.ifft(f,norm='ortho')
#frequency normalization factor is 2*np.pi/dt
w = np.fft.fftfreq(f.size)*2*np.pi/dt


#In order to get a discretisation of the continuous Fourier transform
#we need to multiply g by a phase factor
g*=dt*np.exp(-complex(0,1)*w*t0)/(np.sqrt(2*np.pi))*np.sqrt(len(t))

#Plot Result
pl.plot(w,np.imag(g),"r.")
#For comparison we plot the analytical solution
pl.plot(w,np.imag(exact(w)),color="k")

pl.gca().set_xlim(-10,10)
pl.show()
pl.close()
