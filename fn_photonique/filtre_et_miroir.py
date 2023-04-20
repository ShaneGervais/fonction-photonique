#calcul basé sur les équations de Frenel
import numpy as np

"""
theta0 -> angle d'incidence
n0 -> milieu ambiant
n -> couche_n 
"""
def cos_theta_i(n, n0, theta0):
    
    z1 = np.sqrt(1 - ((n0**2)/(n**2))*(np.sin(theta0))**2)
    resultat = z1
    ncos = n*z1
    
    if np.imag(ncos) < 0: 
        resultat = -z1
    
    return resultat


"""
lambdatemp -> longeur d'onde avec Max de réf.
theta0 -> angle d'incidence
d -> epaisseur du bi-couche (empilement quart-d'onde)
n0 -> milieu ambiant
n -> couche_n 
"""
def P(n, d, n0, theta0, longueur_donde):
    
    phi = 2*np.pi*n*d*cos_theta_i(n, n0, theta0)/longueur_donde

    return np.array([[np.exp(-1j*phi), 0], [0, np.exp(1j*phi)]])

"""
Reflection en P

theta0 -> angle d'incidence
n0 -> milieu ambiant
n1 -> couche_1
n2 -> couche_2
"""
def rp(n1, n2, n0, theta0):
    costheta1 = cos_theta_i(n1,n0,theta0)
    costheta2 = cos_theta_i(n2,n0,theta0)
    
    return (n2*costheta1-n1*costheta2)/(n1*costheta2+n2*costheta1)

"""
Reflection en S

theta0 -> angle d'incidence
n0 -> milieu ambiant
n1 -> couche_1
n2 -> couche_2
"""
def rs(n1, n2, n0, theta0):
    costheta1 = cos_theta_i(n1,n0,theta0)
    costheta2 = cos_theta_i(n2,n0,theta0)
    
    return (n1*costheta1-n2*costheta2)/(n1*costheta1+n2*costheta2)

"""
Transmittance en S

theta0 -> angle d'incidence
n0 -> milieu ambiant
n1 -> couche_1
n2 -> couche_2
"""
def ts(n1, n2, n0, theta0):
    costheta1 = cos_theta_i(n1,n0,theta0)
    costheta2 = cos_theta_i(n2,n0,theta0)
    
    return (2*n1*costheta1)/(n1*costheta1+n2*costheta2)

"""
Transmittance en P

theta0 -> angle d'incidence
n0 -> milieu ambiant
n1 -> couche_1
n2 -> couche_2
"""
def tp(n1, n2, n0, theta0):
    costheta1 = cos_theta_i(n1,n0,theta0)
    costheta2 = cos_theta_i(n2,n0,theta0)
    
    return (2*n1*costheta1)/(n1*costheta2+n2*costheta1)

"""
theta0 -> angle d'incidence
pol -> 1=s  2=p
n0 -> milieu ambiant
n1 -> couche_1
n2 -> couche_2
"""
def T(pol, n1, n2, n0, theta0):
    
    t12s = ts(n1,n2,n0,theta0)
    r12s = rs(n1,n2,n0,theta0)
    t12p = tp(n1,n2,n0,theta0)
    r12p = rp(n1,n2,n0,theta0)

    if pol == 1:
        return np.array([[(1/(t12s))*1, (1/(t12s))*r12s],[(1/(t12s))*r12s, (1/(t12s))*1]])

    else:
        return np.array([[1/(t12p)*1, 1/(t12p)*r12p],[(1/(t12p))*r12p, 1/(t12p)*1]])
