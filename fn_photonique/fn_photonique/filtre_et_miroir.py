"""
Classe de fonction 
"""

import numpy as np

def cos_theta_i(n, n0, theta0):
    
    z1 = np.sqrt(1-(n0/n*np.sin(theta0))^2)    
    resultat = z1
    ncos = n*z1
    
    if np.imag(ncos) < 0: 
        resultat = -z1
    
    return resultat

def P(n, d, n0, theta0, longueur_donde):
    
    phi = 2*np.pi*n*d*cos_theta_i(n, n0, theta0)/longueur_donde

    return [[np.exp(-1j*phi), 0], [0, np.exp(1j*phi)]]

def rp(n1, n2, n0, theta0):
    costheta1 = cos_theta_i(n1,n0,theta0)
    costheta2 = cos_theta_i(n2,n0,theta0)
    
    return (n2*costheta1-n1*costheta2)/(n1*costheta2+n2*costheta1)

def rs(n1, n2, n0, theta0):
    costheta1 = cos_theta_i(n1,n0,theta0)
    costheta2 = cos_theta_i(n2,n0,theta0)
    
    return (n1*costheta1-n2*costheta2)/(n1*costheta1+n2*costheta2)

def ts(n1, n2, n0, theta0):
    costheta1 = cos_theta_i(n1,n0,theta0)
    costheta2 = cos_theta_i(n2,n0,theta0)
    
    return (2*n1*costheta1)/(n1*costheta1+n2*costheta2)

def tp(n1, n2, n0, theta0):
    costheta1 = cos_theta_i(n1,n0,theta0)
    costheta2 = cos_theta_i(n2,n0,theta0)
    
    return (2*n1*costheta1)/(n1*costheta2+n2*costheta1)

def T(pol, n1, n2, n0, theta0):
    
    t12s = ts(n1,n2,n0,theta0)
    r12s = rs(n1,n2,n0,theta0)
    t12p = tp(n1,n2,n0,theta0)
    r12p = rp(n1,n2,n0,theta0)

    if pol == 1:
        return 1/t12s*[[1, r12s],[r12s, 1]]

    else:
        return 1/t12p*[[1, r12p],[r12p, 1]]


