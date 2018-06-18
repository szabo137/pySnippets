from scipy.special import erf
from numpy import *


def g(phi,sigma,Envelope):

    if Envelope == 'cos^2':
        return cos(pi*phi/2./sigma)**2*(phi>-sigma)*(phi<sigma)
    
    elif Envelope == 'cos^4':
        return cos(pi*phi/2./sigma)**4*(phi>-sigma)*(phi<sigma)
    
    elif Envelope == 'Gauss':  
        return exp(-phi**2/2./sigma**2)
        
    elif Envelope == 'Box':  
        return 1.*(phi>=-sigma)*(phi<=sigma)
        
    elif Envelope == 'cosh':
        return 1./cosh(phi/sigma)
        
            
def Int_g_2(phi,sigma,Envelope):
    
    if   Envelope == 'cos^2':
        a = pi/2./sigma
        Wert = (  0.25  * cos(a*phi)**3* sin(a*phi) \
                + 3./8. * cos(a*phi)   * sin(a*phi) \
                + 3./8. * (phi)  *a ) /a         * (phi>=-sigma)*(phi<=sigma)\
                + 3./8. *  sigma *(phi>sigma)

        return Wert 
      
    elif   Envelope == 'cos^4':
        a = pi/2./sigma
        Wert = (  1./8.   * cos(a*phi)**7 * sin(a*phi) + 7. /48.  * cos(a*phi)**5 * sin(a*phi)\
                +35./192. * cos(a*phi)**3 * sin(a*phi) + 35./128. * cos(a*phi)    * sin(a*phi)\
                +35./128. * a * ( phi+sigma ) ) / a                                           *(phi>=-sigma)*(phi<sigma)\
                +35./64.  * sigma                                                             *(phi>=sigma)
        return Wert
        
    elif Envelope == 'Gauss':
        return sqrt(pi)/2.*sigma*erf(phi/sigma)
        
    elif Envelope == 'Box'  :
        return (phi+sigma)*(phi>=-sigma)*(phi<=sigma)+2.*sigma *(phi>sigma) 
     
    elif Envelope == 'cosh':
        return sigma* tanh(phi/sigma)


