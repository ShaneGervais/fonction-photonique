import matplotlib.pyplot as plt
import numpy as np
import scipy as sci
import fonction_polarisation as fp

def fonction_gaussien(sigma, t, t_0, N):
    A_t = np.zeros(N)

    for i in range(len(A_t)):
        A_t[i] = (np.sqrt(1/((np.sqrt(2*np.pi))*sigma)))*np.exp(-((t[i] -t_0)**2)/(4*sigma**2))

    return A_t

def impulsion_gaussien(theta, phi_x, phi_y, w, sigma, t, t_0, N):

    A = fonction_gaussien(sigma, t, t_0, N)
    #champ électrique
    E = np.zeros(len(t))
    
    for i in range(len(t)):
        E[i] = A[i]*fp.etat_polarisation(theta, phi_x, phi_y)*np.exp(1j*w*t[i])
    return E

def impulsion_gaussien_2D(theta, phi_x, phi_y, w, sigma, t, t_0, N):

    #vecteur de jones
    J = fp.j_theta_phi_xy(theta, phi_x, phi_y)
    A = fonction_gaussien(sigma, t, t_0, N)
    #champ électrique
    E = np.zeros([len(J), len(t)])
    
    for i in range(len(J)):
        for j in range(len(t)):
            E[i,j] = A[j]*J[i]*np.exp(1j*w*t[j])

    return E

def intensity(fonction):
    I = np.zeros(shape=(len(fonction)))
    
    for i in range(len(fonction)):
        I[i] = np.abs(np.dot(np.conjugate(fonction[i]), fonction[i]))
    
    return I

def intensity_2D(fonction, n, m):
    I = np.zeros(m)

    for i in range(m):
        I[i] = np.abs(np.dot(np.conjugate(fonction[0, i]), fonction[0, i]))

    return I