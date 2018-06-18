from Config_paar import *
from Envelopefkt import *
from Kinematik import *


def A_m_n(M,N,x_plus,p_el,p_pos,k_photon,k_laser):

    def f1(p):
        return -(m*a0)/(pk(p)) * g(phi,sigma,Envelope) *( pe(1,p) * cos(ksi) * cos(phi) + pe(2,p) * sin(ksi) * sin(phi) )
        
    def f2(p):
        return -(m*a0)**2/(2.*pk(p))*g(phi,sigma,Envelope)**2*((cos(ksi)*cos(phi))**2+(sin(ksi)*sin(phi))**2)          
    
    def f(p):
        return f1(p)+f2(p)
    
    def f1_SVEA(p):
        return -(m*a0)/(pk(p))*g(phi,sigma,Envelope)*(pe(1,p)*cos(ksi)*sin(phi)-pe(2,p)*sin(ksi)*cos(phi))

    def f2_SVEA(p):
        return -(m*a0)**2/(4.*pk(p))*(Int_g_2(phi,sigma,Envelope)+g(phi,sigma,Envelope)**2*cos(phi)*sin(phi)*(cos(ksi)**2-sin(ksi)**2))

    def f_SVEA(p):
        return f1_SVEA(p)+f2_SVEA(p)

    pk  =  lambda imp:    (imp * k_laser)
    pe  =  lambda l,imp:  (imp * eps_laser(l))
    
    P_  =  p_pos.minus() + p_el.minus() - k_photon.minus()    
    s   =  P_/k_laser.minus()

    phi =  w_laser       * x_plus
    
    H_plus =  s*phi - f_SVEA(p_el) + f_SVEA(-p_pos)

    if M == 0:        
        A       = -1./s *  (f(-p_pos) - f(p_el)) * exp(1j  * H_plus)
    
    else:
        A       =  g(phi,sigma,Envelope)**M *exp(  1j* ( H_plus + N*phi))
    
    return A  


def A_m_n_nSVEA(M,N,x_plus,p_el,p_pos,k_photon,k_laser):
        
    def f1(p):
        
        if Envelope == 'cos^2':
            
            fakt_a  = sigma/(sigma-pi)
            fakt_b  = sigma/(sigma+pi)
            
            Int_sin = -0.25 *( fakt_a * cos( phi/fakt_a ) + fakt_b * cos( phi/fakt_b ) +2.*cos(phi) )
            Int_cos =  0.25 *( fakt_a * sin( phi/fakt_a ) + fakt_b * sin( phi/fakt_b ) +2.*sin(phi) )
            
            return -(m*a0)/(pk(p)) *( pe(1,p) * cos(ksi) * Int_cos + pe(2,p) * sin(ksi) * Int_sin )
        
        
        elif Envelope == 'cos^4':
            
            fakt_a  = lambda n: (  1. + n*pi/sigma )
            fakt_b  = lambda n: ( -1. + n*pi/sigma )
            
            Int_sin =  0.25 *( ( - cos( fakt_a(2.)*phi ) / fakt_a(2.) + cos( fakt_b(2.)*phi ) / fakt_b(2.) ) * 0.25 \
                                 - cos( fakt_a(1.)*phi ) / fakt_a(1.) + cos( fakt_b(1.)*phi ) / fakt_b(1.)   - 3./2. * cos(phi) )
            
            Int_cos =  0.25 *( (   sin( fakt_a(2.)*phi ) / fakt_a(2.) + sin( fakt_b(2.)*phi ) / fakt_b(2.) ) * 0.25 \
                                 + sin( fakt_a(1.)*phi ) / fakt_a(1.) + sin( fakt_b(1.)*phi ) / fakt_b(1.)   - 3./2. * sin(phi) )
        
            return -(m*a0)/(pk(p)) * ( pe(1,p) * cos(ksi) * Int_cos + pe(2,p) * sin(ksi) * Int_sin )
            
        
        elif Envelope == 'cosh':
            raise IOError,'cosh noch nicht implementiert -> benutze SEVA'
        
        else:
            raise IOError,'Nicht analytisch loesbar -> benutze SEVA'
    
    
    
    def f2(p):
        if Envelope == 'cos^2':
            
            a       =  pi/sigma/2.
            F       = lambda l,n: (  l + n*a )
            
            
            Int_cos = 1./8. *( 1.5*phi + 0.75*sin(2.*phi)  + sin(F(0,4.)*phi)/F(0,8.) + sin(F(0,2.)*phi)/a  \
                            + sin(F(-2.,4.)*phi)/F(-2.,4.)/4. + sin(F(2.,4.)*phi)/F(2.,4.)/4.               \
                            + sin(F(-2.,2.)*phi)/F(-2.,2.)    + sin(F(2.,2.)*phi)/F(2.,2.)                  )
                            
            Int_sin = 1./8. *( 1.5*phi - 0.75*sin(2.*phi)  + sin(F(0,4.)*phi)/F(0,8.) + sin(F(0,2.)*phi)/a  \
                            - sin(F(-2.,4.)*phi)/F(-2.,4.)/4. - sin(F(2.,4.)*phi)/F(2.,4.)/4.               \
                            - sin(F(-2.,2.)*phi)/F(-2.,2.)    - sin(F(2.,2.)*phi)/F(2.,2.)                  )
            
            return -( m*a0 )**2 / (2.*pk(p)) * ( cos(ksi)**2 * Int_cos +  sin(ksi)**2 * Int_sin )          
        
        elif Envelope == 'cos^4':
            
            
            Faktor  = lambda l,n: (  l + n*pi/sigma )
            
            Int_sin =  1./64. *( (- sin( Faktor(-2.,4.)*phi ) / Faktor(-2.,4.) - sin( Faktor(2.,4.)*phi ) / Faktor(2.,4.) ) / 8.  \
                                 -  sin( Faktor(-2.,3.)*phi ) / Faktor(-2.,3.) - sin( Faktor(2.,3.)*phi ) / Faktor(2.,3.)         \
                                 -( sin( Faktor(-2.,2.)*phi ) / Faktor(-2.,2.) + sin( Faktor(2.,2.)*phi ) / Faktor(2.,2.) ) * 3.5 \
                                 -( sin( Faktor(-2.,1.)*phi ) / Faktor(-2.,1.) + sin( Faktor(2.,1.)*phi ) / Faktor(2.,1.) ) * 7.  \
                                 +  sin( Faktor( 0.,4.)*phi ) / Faktor(0.,16.) + sin( Faktor(0.,3.)*phi ) / Faktor(0.,1.5)        \
                                 +( sin( Faktor( 0.,2.)*phi ) / Faktor(0.,2.)  + sin( Faktor(0.,1.)*phi ) / Faktor(0.,0.5)) * 7.  \
                                 +   35./4. * phi - 35./8. * sin( 2*phi ) ) 
            
            Int_cos =  1./64. *( (  sin( Faktor(-2.,4.)*phi ) / Faktor(-2.,4.) + sin( Faktor(2.,4.)*phi ) / Faktor(2.,4.) ) / 8.  \
                                 +  sin( Faktor(-2.,3.)*phi ) / Faktor(-2.,3.) + sin( Faktor(2.,3.)*phi ) / Faktor(2.,3.)         \
                                 +( sin( Faktor(-2.,2.)*phi ) / Faktor(-2.,2.) + sin( Faktor(2.,2.)*phi ) / Faktor(2.,2.) ) * 3.5 \
                                 +( sin( Faktor(-2.,1.)*phi ) / Faktor(-2.,1.) + sin( Faktor(2.,1.)*phi ) / Faktor(2.,1.) ) * 7.  \
                                 +  sin( Faktor( 0.,4.)*phi ) / Faktor(0.,16.) + sin( Faktor(0.,3.)*phi ) / Faktor(0.,1.5)        \
                                 +( sin( Faktor( 0.,2.)*phi ) / Faktor(0.,2.)  + sin( Faktor(0.,1.)*phi ) / Faktor(0.,0.5)) * 7.  \
                                 +   35./4. * phi - 35./8. * sin( 2*phi )  )
            
            return -( m*a0 )**2 / (2.*pk(p)) * ( cos(ksi)**2 * Int_cos +  sin(ksi)**2 * Int_sin ) 
            
         
        elif Envelope == 'cosh':
            raise IOError,'cosh noch nicht implementiert -> benutze SEVA'
        
        else:
            raise IOError,'Nicht analytisch loesbar -> benutze SEVA'
    
    def f(p):
        return f1(p)+f2(p)
    

    pk  =  lambda imp:    (imp * k_laser)
    pe  =  lambda l,imp:  (imp * eps_laser(l))
    
    P_  =  p_pos.minus() + p_el.minus() - k_photon.minus()    
    s   =  P_/k_laser.minus()
    
    
    
    phi =  w_laser       * x_plus
    
    H_plus =  s*phi - f(p_el) + f(-p_pos)

    A       =  g(phi,sigma,Envelope)**M *exp(  1j* ( H_plus + N*phi))
    
    return A 
    
    
def A_0_0  (A11,A1_1,A20,A22,A2_2):
    

    
    p_pos,p_el,k_laser,k_photon,q_pos,eps_m,eps_p = kinematik()
    
    pk  =  lambda p:    (p * k_laser)
    d_p = lambda p: m*a0 / ( 4.* pk(p) )    
        
    P_  =  p_pos.minus() + p_el.minus() - k_photon.minus()    
    s   =  P_/k_laser.minus()
        
    Wert = 2./s * (  ( d_p(p_pos)*p_pos*eps_m - d_p(p_el)*p_el*eps_m ) * A11  \
                   + ( d_p(p_pos)*p_pos*eps_p - d_p(p_el)*p_el*eps_p ) * A1_1 \
                   -  k_laser*k_photon*d_p(p_pos)*d_p(p_el)                   \
                   * ( 2.*A20 + (cos(ksi)**2 - sin(ksi)**2) * (A22 + A2_2) )  )
    
    return Wert


