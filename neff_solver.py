#Author: Alexander Greenwood
#The following scripts have been used for my BASc thesis on a "Segmented Microring Modulator"


#Based on the following sources:
#Silicon Photonics Design
#By Professor Lukas Chrostowski, Professor Michael Hochberg
#Photonics
#By Prof. Amnon Yariv, Prof. Pochi Yeh


#This is a work in progress :)

import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#First define some physical constants
n_i = 1.1*10**10 #intrinsic carrier concentration [cm^-3]
k_b = 1.38*10**(-23) #Boltzmann's constant [JK^-1]
epsilon_o = 8.85*10**(-14) #permittivity of free space [F/m] 
q = 1.602*10**(-19) # charge of an electron [C]
temp = 300 # [K]
pi = 3.14
c = 2.98*10**(8)
mu=4*pi*1e-7
#############################################################
def get_neff(d, wavelength, n1, n2):
#Obtain effective index profile for our waveguide
    k0 = 2*pi / wavelength
    k = n2*k0
    beta = np.linspace(n1*k0,n2*k0, num=1000)# Our propagation constant is somewhere between these two spatial freq's 
    beta = beta[:-1]
    h = np.sqrt(k**2 - beta**2) 
    q = np.sqrt(beta**2 - (k0*n1)**2)
    te_rel = np.tan(h*d)*(h**2 - q**2) - 2*h*q # The zeros of this relation correspond to our modes!
    zeros = list() 
    intervals = (te_rel>=0).astype(int) - (te_rel<0).astype(int)
    ############################################
    #find indexes of zeros!
    for i in range(1,te_rel.size):
        if(intervals[i]-intervals[i-1]>0):
            zeros.append(i)

    #zeros = zeros[:-1]
    neff = beta[zeros] / k0
    h_m = h[zeros]
    q_m = q[zeros]
    print(neff)

    return [neff,h_m,q_m,beta[zeros],k0]


############################################################# 
def get_coupling_coeff(coupling_length,spacing,d, wavelength, n1, n2): 
    print("This is a placeholder")

def get_field_profile(d, wavelength, n1, n2): 
    #In the Yariv text, we think of the x axis as normal to the substrate. 
    
    #Much of the following code has been taken from the Matlab scripts provided in Silicon Photonics Design by Chrostowski and Hochberg
    
    pts = 100
    M = 4 #Our simulation window is bounded between -4d and 4d
    x1=np.linspace( -M*d, -d/2, pts) 
    x2=np.linspace( -d/2, d/2, pts)
    x3=np.linspace( d/2, M*d, pts)
    x=[x1, x2, x3]
    nx=np.concatenate((n1*np.ones(pts), n2*np.ones(pts), n1*np.ones(pts)),axis=None)
     
    print(nx)
    mode_results = get_neff(220e-9,1.55e-6,1.444,3.47) #mode_results = [neff,h_m,q_m,beta,k]

    neff = mode_results[0]
    h_m = mode_results[1]
    q_m = mode_results[2]
    beta = mode_results[3]
    k0 = mode_results[4]
    omega = c / n2 * k0
    
    
    #Define normalization coefficient C:
    #More information on this can be found in Chapter 3.2 of Photonics by Yariv and Yeh
    C = 2 * h_m * np.power((omega * mu / (beta * (d + 1/q_m + 1/q_m)*(h_m**2 + q_m**2))),0.5)
    print(C)
    E = np.zeros(pts*3)
    for i in range(0,len(C)):
        #region 1 E field
        E1 = C[i]*np.exp(q_m*(x1+d/2))
        #region 2 E field
        E2 = C[i]*(np.cos(h_m[i]*(x2+d/2))-q/h_m[i]*np.sin(h_m[i]*(x2+d/2)))

        #region 3 E field. We're assuming that the refractive indexes of the BOX and cladding are identical.
        E3 = C[i]*(np.cos(h_m[i]*d)+q_m[i] / h_m[i] * np.sin(h_m[i]*d))*np.exp(-h_m[i]*(x3-d/2))

        E += np.concatenate((E1,E2,E3),axis=None)

    x = np.concatenate((x1,x2,x3), axis=None)
    plt.plot(x,E)
    plt.show()


    





get_field_profile(220e-9,1.55e-6,1.444,3.47)