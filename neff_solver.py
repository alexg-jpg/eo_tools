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
#First define some physical constants
n_i = 1.1*10**10 #intrinsic carrier concentration [cm^-3]
k_b = 1.38*10**(-23) #Boltzmann's constant [JK^-1]
epsilon_o = 8.85*10**(-14) #permittivity of free space [F/m] 
q = 1.602*10**(-19) # charge of an electron [C]
temp = 300 # [K]
pi = 3.14
#############################################################
def get_neff(d, wavelength, n1, n2):
#Obtain effective index profile for our waveguide
    k0 = 2*pi / wavelength
    k = n2*k0
    beta_0 = np.linspace(n1*k0,n2*k0, num=1000)# Our propagation constant is somewhere between these two spatial freq's 
    h = np.sqrt(k**2 - beta_0**2) 
    q = np.sqrt(beta_0**2 - (k0*n1)**2)
    te_rel = np.tan(h*d)*(h**2 - q**2) - 2*h*q # The zeros of this relation correspond to our modes!
    zeros = list() 
    intervals = (te_rel>=0).astype(int) - (te_rel<0).astype(int)
    ############################################
    #find indexes of zeros!
    for i in range(1,te_rel.size):
        if(intervals[i]-intervals[i-1]>0):
            zeros.append(i)

    neff = beta_0[zeros] / k0
    print(neff)

############################################################# 
def get_coupling_coeff(coupling_length,spacing,d, wavelength, n1, n2): 
    print("This is a placeholder")



get_neff(220e-9,1.55e-6,1.44,3.47)