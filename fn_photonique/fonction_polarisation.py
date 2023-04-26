import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Graph le vecteur de Jones donné
def graph_j(J):
    #normalise J -> un vecteur 
    J =J/np.linalg.norm(J, 2)

    phi = np.linspace(0, 2*np.pi)

    #variable x, y
    x = J[0]*np.exp(-1j*phi)
    y = J[1]*np.exp(-1j*phi)


    """dr = [i for i in range(2)]
    dr = [x[1] -x[0] , y[1] - y[0]]
    dr = (dr/np.linalg .norm(dr, 2))*0.5
    """
    dr = np.array([x[1], y[1]]) - np.array([x[0], y[0]])
    dr = dr / np.linalg.norm(dr) * 0.5

    plt.plot(np.real(x), np.real(y))
    plt.quiver((x[0]), (y[0]), dr[0], dr[1])
    plt.axis([-1, 1, -1, 1])
    plt.show()

#fonction qui trouve les vecteurs de jones à partir des paramètres theta et epsilon
#-pi/4 < epsilone pi/4 polarisation circulaire epsi = (+/-)pi/4; linéaire = 0 
def j_theta_epsilon(theta, epsilon):
    #vecteur de jones
    J = [i for i in range(2)]
    J[0] = np.cos(theta)*np.cos(epsilon)+1j*np.sin(theta)*np.sin(epsilon)
    J[1] = np.sin(theta)*np.cos(epsilon)-1j*np.cos(theta)*np.sin(epsilon)

    phi = np.angle(J[0])

    return [J[i]*np.exp(-1j*phi) for i in range(2)]

#un polariseur orienté a un angle phi
#Polarisation linéaire: différence de phase = 0
#Polarisation circulaire: différence de phase = pi/2
#Polarisation élliptique: différence de phase 0 < delta < pi/2
def polariseur(phi):
    M = np.zeros((2,2))
    M[0][0] = np.cos(phi)**2
    M[0][1] = np.cos(phi)*np.sin(phi)
    M[1][0] = M[0][1]
    M[1][1] = np.sin(phi)**2
    return M

#trouve le psi(l'angle) et delta (différence de phase) d'un vecteur de Jones
def psidelta_j(J):

    E_0 = np.dot(J, J)
    E_0x = np.linalg.norm(J[0])
    E_0y = np.linalg.norm(J[1])
    #psi = np.pi/2

    #if E_0x == 0:
    
    psi = np.arctan2(E_0y,E_0x)

    delta = np.angle(J[1]/J[0])

    return psi, delta

#Comparaison de deux vecteur de Jones
def compare_jones(J_0, J_1):
    J_0 = J_0/np.linalg.norm(J_0)
    J_1 = J_1/np.linalg.norm(J_1)

    return np.abs(np.sum(np.conj(J_0)*J_1))

#Un retardeur orienté avec un angle phi
#Pour un lame quart-onde la difference de phase = +/-(pi/2)
#Pour un lame demi-onde la difference de phase = +/-(pi)
def retardeur(phi_x, phi_y, phi):
    a = np.exp(1j*phi_x)*np.cos(phi)**2 + np.exp(1j*phi_y)*np.sin(phi)**2
    b = (np.exp(1j*phi_x) - np.exp(1j*phi_y))*np.cos(phi)*np.sin(phi)
    c = b
    d = np.exp(1j*phi_x)*np.sin(phi)**2 + np.exp(1j*phi_y)*np.cos(phi)**2

    return [[a, b], [c, d]]

#Trouve le sens de rotation, élliptisité et l'angle
def theta_epsilon_sens_j(J):
    
    [psi, delta] = psidelta_j(J)
    
    theta = (1/2)*np.arctan2(np.sin(2*psi)*np.cos(delta), np.cos(2*psi))
    epsilon= (-1/2)*np.arcsin(np.sin(2*psi)*np.sin(delta));
    sens = np.sign(epsilon)
    return theta, epsilon, sens

#trace sur la sphère poincaré
def trace_points_sphere_poincare(r):
    # Create a 3D sphere plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')

    # Create a sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)

    # Plot the sphere
    #ax.plot_surface(x, y, z, color='b')
    
    # Get the number of rows in the matrix r
    N = r.shape[0]
    
    # Plot each point in r as a red dot
    for i in range(N):
        ax.scatter(r[i, 0], r[i, 1], r[i, 2], s=50, c='r', marker='o')
    
    # Add labels to the axes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    # Set the axes limits to ensure the sphere is plotted as a sphere, not an ellipsoid
    max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max()
    x_mean = x.mean()
    y_mean = y.mean()
    z_mean = z.mean()
    ax.set_xlim(x_mean - max_range/2, x_mean + max_range/2)
    ax.set_ylim(y_mean - max_range/2, y_mean + max_range/2)
    ax.set_zlim(z_mean - max_range/2, z_mean + max_range/2)
    

#Trouve X, Y et Z de la sphère
def XYZ_poincare_j(J):

    theta, epsilon, sens = theta_epsilon_sens_j(J)

    # Initialize an array r with three elements
    r = np.zeros(3)
    
    # Calculate the x, y, and z coordinates of the Poincare sphere
    r[0] = np.cos(2*theta)*np.cos(2*epsilon)
    r[1] = np.sin(2*theta)*np.cos(2*epsilon)
    r[2] = np.sin(2*epsilon)
    
    return r

#Trouve l'angle et l'éllipticité de la sphère
def XYZ_poincare_theta_epsi(theta, epsilon):
    
    # Initialize an array r with three elements
    r = np.zeros(3)
    
    # Calculate the x, y, and z coordinates of the Poincare sphere
    r[0] = np.cos(2*theta)*np.cos(2*epsilon)
    r[1] = np.sin(2*theta)*np.cos(2*epsilon)
    r[2] = np.sin(2*epsilon)
    
    return r
