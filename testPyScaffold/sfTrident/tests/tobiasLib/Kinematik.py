from Config_paar import *

m=1.0
def kinematik():

    k_laser       = MinkowskiVector([w_laser,0,0,-w_laser])
    n1            = MinkowskiVector([1, sin(theta_ph)*cos(phi_ph), sin(theta_ph)*sin(phi_ph), cos(theta_ph)])
    k_photon      = w_photon*n1

    p_pos  = MinkowskiVector([ m_perp*cosh(y_z),
                            p_perp*cos(phi_pos),
                            p_perp*sin(phi_pos),
                            m_perp*sinh(y_z)] )

    q_pos      = p_pos + (m**2*a0**2)/(4.*(k_laser*p_pos))*k_laser # eff Mom Pos


    p_el_plus  =     k_photon.plus_component()     -     p_pos.plus_component()
    p_el_1     =     k_photon._1()                 -     p_pos._1()
    p_el_2     =     k_photon._2()                 -     p_pos._2()

    p_el_minus = (m**2+p_el_1**2+p_el_2**2)/p_el_plus

    p_el_0     = 0.5*(p_el_plus + p_el_minus)
    p_el_3     = 0.5*(p_el_plus - p_el_minus)
    p_el       = MinkowskiVector([ p_el_0, p_el_1 ,p_el_2, p_el_3 ])


    eps_m      = eps_laser(1)*cos(ksi)-1j*eps_laser(2)*sin(ksi)
    eps_p      = eps_laser(1)*cos(ksi)+1j*eps_laser(2)*sin(ksi)


    #print "kp: %s"%(k_laser)
    #print "k: %s"%(k_photon)
    #print "pp: %s"%(p_pos)
    #print "p: %s"%(p_el)


    return (p_pos,p_el,k_laser,k_photon,q_pos,eps_m,eps_p)
