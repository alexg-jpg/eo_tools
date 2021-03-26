#Author: Alexander Greenwood
#The following scripts have been used for my BASc thesis on a "Segmented Microring Modulator"


#Based on the following sources:
#Silicon Photonics Design
#By Professor Lukas cTEhrostowski, Professor Michael Hochberg
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
q = 1.602*10**(-19) # charge of an electron [cTE]
temp = 300 # [K]
pi = 3.14
c = 3*10**(8)
mu=4*pi*1e-7
#############################################################
def GetNeffTE(d, wavelength, n1, n2):
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

def GetNeffTM(d, wavelength, n1, n2):
#Obtain effective index profile for our waveguide
    k0 = 2*pi / wavelength
    k = n2*k0
    beta = np.linspace(n1*k0,n2*k0, num=1000)# Our propagation constant is somewhere between these two spatial freq's 
    beta = beta[:-1]
    h = np.sqrt(k**2 - beta**2) 
    q = np.sqrt(beta**2 - (k0*n1)**2)

    qBar=(n2**2)/(n1**2)*q 

    tm_rel = np.tan(h*d)*(h**2 - qBar**2) - 2*h*qBar # The zeros of this relation correspond to our modes!
    zeros = list() 
    intervals = (tm_rel>=0).astype(int) - (tm_rel<0).astype(int)
    ############################################
    #find indexes of zeros!
    for i in range(1,tm_rel.size):
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
    
    #Much of the following code has been taken from the Matlab scripts provided in Silicon Photonics Design by cTEhrostowski and Hochberg
    
    
    eta = np.sqrt(mu/epsilon_o)# The "characteristic impedance of the medium"
    pts = 100
    M = 4 #Our simulation window is bounded between -4d and 4d
    x1=np.linspace( -M*d, -d/2, pts) 
    x2=np.linspace( -d/2, d/2, pts)
    x3=np.linspace( d/2, M*d, pts)
    x=[x1, x2, x3]
    nx=np.concatenate((n1*np.ones(pts), n2*np.ones(pts), n1*np.ones(pts)),axis=None)
     
    print(nx)
    modeResultsTE = GetNeffTE(220e-9,1.55e-6,1.444,3.47) #modeResultsTE = [neff,h_m,q_m,beta,k]
    modeResultsTM = GetNeffTM(220e-9,1.55e-6,1.444,3.47) #modeResultsTE = [neff,h_m,q_m,beta,k]

    neffTE = modeResultsTE[0]
    hTE = modeResultsTE[1]
    qTE = modeResultsTE[2]
    qBarTE=(n2**2)/(n1**2)*qTE #Syntax taken from Matlab script by Lukas cTEhrostowski, 2012
    betaTE = modeResultsTE[3]
    k0 = modeResultsTE[4]
    omega0 = c * k0

    neffTM = modeResultsTM[0]
    hTM = modeResultsTM[1]
    qTM = modeResultsTM[2]
    qBarTM=(n2**2)/(n1**2)*qTM #Syntax taken from Matlab script by Lukas cTEhrostowski, 2012
    betaTM = modeResultsTM[3]
    
    
    #Define normalization coefficient cTE:
    #Each mode with field E_m has power flow of 1 W.
    #More information on this can be found in cTEhapter 3.2 of Photonics by Yariv and Yeh
    cTE = 2 * hTE * np.sqrt((omega0 * mu / (betaTE * (d + 1/qTE + 1/qTE)*(hTE**2 + qTE**2))))
    
    tEff = (qBarTM**2 + hTM**2) / (qBarTM**2)*(d/(n2**2) + (qTM**2 + hTM**2)/(qBarTM**2 + hTM**2)*((1/(n1**2*qTM)+(1/(n2**2*qTM)))))
    cTM = 2 * np.sqrt(omega0*epsilon_o / (np.abs(betaTM)*tEff))
    print(cTE)
    eTE = np.zeros(pts*3)
    hFieldTM = np.zeros(pts*3)
    for i in range(0,len(cTE)):
        #region 1 E field
        eTE1 = cTE[i]*hTE[i]/qBarTE[i]*np.exp(qTE[i]*(x1+d/2))
        #region 2 E field
        eTE2 = cTE[i]*(hTE[i]/qBarTE[i]*np.cos(hTE[i]*(x2+d/2))+np.sin(hTE[i]*(x2+d/2)))

        #region 3 E field. We're assuming that the refractive indexes of the BOX and cladding are identical.
        eTE3 = cTE[i]*(hTE[i]/qBarTE[i]*np.cos(hTE[i]*d)+np.sin(hTE[i]*d))*np.exp(-qTE[i]*(x3-d/2))

        eTE += np.concatenate((eTE1,eTE2,eTE3),axis=None)

    for i in range(0,len(cTM)):
        #region 1 E field
        hTM1 = cTM[i]*hTM[i]/qBarTM[i]*np.exp(qTM[i]*(x1+d/2))
        #region 2 E field
        hTM2 = cTM[i]*(hTM[i]/qBarTM[i]*np.cos(hTM[i]*(x2+d/2))+np.sin(hTM[i]*(x2+d/2)))

        #region 3 E field. We're assuming that the refractive indexes of the BOX and cladding are identical.
        hTM3 = cTM[i]*(hTM[i]/qBarTM[i]*np.cos(hTM[i]*d)+np.sin(hTM[i]*d))*np.exp(-qTM[i]*(x3-d/2))

        hFieldTM += np.concatenate((hTM1,hTM2,hTM3),axis=None)
    
    eFieldTM = hFieldTM / nx * eta
    x = np.concatenate((x1,x2,x3), axis=None)
    #plt.plot(x,eTE)
    #plt.show()

    plt.plot(x,hFieldTM)
    plt.show()
    return [x,eTE]

''''
def get_neffV(n_a,n_d,v_a,w_pn, w_rib,w_d, w_dn, w_dp, x_offset,d, wavelength, n1, n2, v_a):

    #begin by getting 1D transverse field profile of our WG.
    [x,E] = get_field_profile(d,wavelength,n1,n2)

    #Obtain 1D transverse carrier profile of our WG.
    [n_profile, p_profile] = get_carrier_concentrations(n_a,n_d,v_a, w_pn, w_rib,w_d, w_dn, w_dp, x_offset)

    #cTEompute overlap integrals...
    
    return
'''





get_field_profile(220e-9,1.55e-6,1.444,3.47)