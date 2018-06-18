from numpy import pi

m            = 0.511
m_elec       = 0.511                                   # Masse Elektron in MeV
m_elec_SI    = 9.10938291e-31                          # in kg
m_heinzel    = 0.5

alpha_fein   = 1./137.036
e2           = 4.*pi*alpha_fein                         # e^2 in natuerlichen Einheiten
e_elec_SI    = -1.602176565e-19                         # in C 
r_0          = alpha_fein / m_elec                      # klassischer e- Radius

c_licht      = 299792458                                # Lichtgeschw. in SI m/s
h_quer       = 6.58211928e-16                           # Planck. Wirk. in SI eV*s
epsilon_0    = 8.85418787162e-12                        # Elektr. Feldkonstante in A*s/V/m 

U_MeV_2_mb   = (c_licht * h_quer * 1e9)**2 * 10         # Umrechnung von MeV^-2 zu mb
 

if __name__ == '__main__': 

    from numpy import *
    
    faktor  = h_quer * c_licht * 1e9                    #   MeV * fm 
    scale   = 1e-1                                      # 1 mb in fm^2
    
    l       = 1                                         # lange in mb
    energie = 1                                         # in MeV^-2    
    
    #print  energie / scale * faktor**2,'mb'
    #print l * scale / faktor**2,'MeV^-2'
    
    E_feld  =   1e17
    lambda_ =   0.05e-9
    omega   =   2*pi*h_quer*c_licht/lambda_
    nu      =   1e4/0.511e6
    
    E_c     =    m_elec_SI**2 * c_licht**3 / h_quer/ e_elec_SI**2
    
    #print  'omega  = ',omega
    #print 0.5 * epsilon_0 * c_licht 
        
    #print 'a0 = ',E_feld/E_c/nu
    #print 
    #print e_elec_SI/c_licht**2
    
    #print 'E_c = ',E_c
    #print m_elec**2*1e12/4./pi/alpha_fein/c_licht**3/h_quer**3*1e-18
    

    
    