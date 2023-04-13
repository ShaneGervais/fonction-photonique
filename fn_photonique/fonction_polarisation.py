import numpy as np
import matplotlib.pyplot as plt
import math

def graph_j(J):
    #normalise J -> un vecteur 
    J =J/np.linalg.norm(J, 2)

    phi = np.linspace(0, 2*np.pi)

    #variable x, y
    x = J[0]*np.exp(-1j*phi)
    y = J[1]*np.exp(-1j*phi)

    #plt.show()

    dr = [i for i in range(2)]
    dr = [x[1] -x[0] , y[1] - y[0]]
    dr = (dr/np.linalg .norm(dr, 2))*0.5

    plt.plot(np.real(x), np.real(y))
    plt.quiver(np.real(x[1]), np.real(y[1]), dr[0], dr[1])
    plt.show()



#fonction qui trouve les vecteurs de jones a partir des parametres theta et epsilon
def j_theta_epsilon(theta, epsilon):
    #vecteur de jones
    J = [i for i in range(2)]
    J[0] = np.cos(theta)*np.cos(epsilon)+1j*np.sin(theta)*np.sin(epsilon)
    J[1] = np.sin(theta)*np.cos(epsilon)-1j*np.cos(theta)*np.sin(epsilon)

    phi = np.angle(J[0])

    return [J[i]*np.exp(-1j*phi) for i in range(2)]

def polariseur(phi):
    M = np.zeros((2,2))
    M[0][0] = np.cos(phi)**2
    M[0][1] = np.cos(phi)*np.sin(phi)
    M[1][0] = M[0][1]
    M[1][1] = np.sin(phi)**2
    return M


def psidelta_j(J):
    E_0 = np.dot(J, J)
    E_0x = np.linalg.norm(J[0])
    E_0y = np.linalg.norm(J[1])
    psi = np.pi/2

    if E_0x == 0:
        psi = math.atan(E_0y/E_0x)

    delta = np.angle(J[1]/J[0])
    

    return psi, delta

def compare_jones(J_0, J_1):
    J_0 = J_0/np.linalg.norm(J_0)
    J_1 = J_1/np.linalg.norm(J_1)

    return np.abs(np.sum(np.conj(J_0)*J_1))

def retardeur(phi_x, phi_y, phi):
    a = np.exp(1j*phi_x)*np.cos(phi)**2 + np.exp(1j*phi_y)*np.sin(phi)**2
    b = (np.exp(1j*phi_x) - np.exp(1j*phi_y))*np.cos(phi)*np.sin(phi)
    c = b
    d = np.exp(1j*phi_x)*np.sin(phi)**2 + np.exp(1j*phi_y)*np.cos(phi)**2

    return [[a, b], [c, d]]

def theta_epsilon_sens_j(J):
    
    [psi, delta] = psidelta_j(J)
    
    theta = (1/2)*np.arctan2(np.sin(2*psi)*np.cos(delta), np.cos(2*psi))
    epsilon= (-1/2)*np.arcsin(np.sin(2*psi)*np.sin(delta));
    sens = np.sign(epsilon)
    return theta, epsilon, sens

def trace_points_sphere_poincare(r):
    N = np.size(r, 1)

    for i in range(N):
        plt.scatter(r[i, 0], r[i, 1], r[i, 2])
    
    plt.xlabel('X');
    plt.ylabel('Y');
    plt.zlabel('Z');

def XYZ_poincare_j(J):

    theta, epsilon, sens = theta_epsilon_sens_j(J)

    r = [i for i in range(2)]
    r[0] = np.cos(2*theta)*np.cos(2*epsilon)
    r[1] = np.sin(2*theta)*np.cos(2*epsilon)
    r[2] = np.sin(2*epsilon)

    return r

def XYZ_poincare_theta_epsi(theta, epsilon):
    
    r = [i for i in range(2)]
    r[0] = np.cos(2*theta)*np.cos(2*epsilon)
    r[1] = np.sin(2*theta)*np.cos(2*epsilon)
    r[2] = np.sin(2*epsilon)

    return r
    