#Author: Alexander Greenwood
#The following scripts have been used for my BASc thesis on a "Segmented Microring Modulator"


#Based on the following sources:
#Silicon Photonics Design
#By Professor Lukas cTEhrostowski, Professor Michael Hochberg
#Photonics
#By Prof. Amnon Yariv, Prof. Pochi Yeh


#This is a work in progress :)
from electrical_characteristics import *
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
def get_neff_te(d, wavelength, n1, n2,n3):
#Obtain effective index profile for our waveguide
    k0 = 2*pi / wavelength
    k = n2*k0
    beta = np.linspace(n1*k0,n2*k0, num=1000)# Our propagation constant is somewhere between these two spatial freq's 
    beta = beta[:-1]
    h = np.sqrt(k**2 - beta**2) 
    q = np.sqrt(beta**2 - (k0*n1)**2)
    p = np.sqrt(beta**2 - (k0*n3)**2)
    
    te_rel = np.tan(h*d) - ((p+q)/h)/(1-((p*q)/(h**2)))
    #te_rel = np.tan(h*d)*(h**2 - p*q) - (p+q)*h # The zeros of this relation correspond to our modes!
    zeros = list() 
    intervals = (te_rel>=0).astype(int) - (te_rel<0).astype(int)
    ############################################
    #find indexes of zeros!
    for i in range(1,te_rel.size):
        if(intervals[i]-intervals[i-1]<0):
            zeros.append(i)

    #zeros = zeros[:-1]
    neff = beta[zeros] / k0
    h_m = h[zeros]
    q_m = q[zeros]
    

    return [neff,h_m,q_m,beta[zeros],k0]

def get_neff_tm(d, wavelength, n1, n2,n3):
#Obtain effective index profile for our waveguide
    k0 = 2*pi / wavelength
    k = n2*k0
    beta = np.linspace(n1*k0,n2*k0, num=1000)# Our propagation constant is somewhere between these two spatial freq's 
    beta = beta[:-1]
    h = np.sqrt(k**2 - beta**2) 
    q = np.sqrt(beta**2 - (k0*n1)**2)
    p = np.sqrt(beta**2 - (k0*n3)**2)

    qBar=(n2**2)/(n1**2)*q
    pBar=(n2**2)/(n3**2)*q

    tm_rel = np.tan(h*d) - h*(pBar + qBar) / (h**2 - pBar*qBar) # The zeros of this relation correspond to our modes!
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

def get_field_profile(d, w_rib, h_slab, wavelength, n1, n2): 
    #In the Yariv text, we think of the x axis as normal to the substrate. 
    
    #Much of the following code has been taken from the Matlab scripts provided in Silicon Photonics Design by cTEhrostowski and Hochberg
    
    
    eta = np.sqrt(mu/epsilon_o)# The "characteristic impedance of the medium"
    pts = 500
    M = 4 #Our simulation window is bounded between -4d and 4d
    x1=np.linspace( -M*d, -d/2, pts) 
    x2=np.linspace( -d/2, d/2, pts)
    x3=np.linspace( d/2, M*d, pts)

    y1=np.linspace( -M*w_rib, -w_rib/2, pts) 
    y2=np.linspace( -w_rib/2, w_rib/2, pts)
    y3=np.linspace( w_rib/2, M*w_rib, pts)
     
    modeResultsTECentre = get_neff_te(d,1.55e-6,1.444,3.47,1) #modeResultsTE = [neff,h_m,q_m,beta,k]
    modeResultsTEEdge = get_neff_te(h_slab,1.55e-6,1.444,3.47,1) #modeResultsTE = [neff,h_m,q_m,beta,k]

    modeResultsTE = get_neff_te(w_rib,1.55e-6,modeResultsTEEdge[0],modeResultsTECentre[0],modeResultsTEEdge[0]) 
    modeResultsTM = get_neff_tm(w_rib,1.55e-6,modeResultsTEEdge[0],modeResultsTECentre[0],modeResultsTEEdge[0]) 
    '''
    neffTE = modeResultsTE[0]
    hTE = modeResultsTE[1]
    qTE = modeResultsTE[2]
    qBarTE=(n2**2)/(n1**2)*qTE #Syntax taken from Matlab script by Lukas cTEhrostowski, 2012
    betaTE = modeResultsTE[3]
    k0 = modeResultsTE[4]
    '''

    k0 = modeResultsTM[4]
    omega0 = c * k0
    

    neffTM = modeResultsTM[0]
    hTM = modeResultsTM[1]
    qTM = modeResultsTM[2]
    qBarTM=(modeResultsTECentre[0]**2)/(modeResultsTEEdge[0]**2)*qTM #Syntax taken from Matlab script by Lukas Chrostowski, 2012
    betaTM = modeResultsTM[3]
    
    
    #Define normalization coefficient cTE:
    #Each mode with field E_m has power flow of 1 W.
    #More information on this can be found in cTEhapter 3.2 of Photonics by Yariv and Yeh
    
    #cTE = 2 * hTE * np.sqrt((omega0 * mu / (betaTE * (d + 1/qTE + 1/qTE)*(hTE**2 + qTE**2))))
    
    tEff = (qBarTM**2 + hTM**2) / (qBarTM**2)*(w_rib/(modeResultsTECentre[0]**2) + (qTM**2 + hTM**2)/(qBarTM**2 + hTM**2)*((1/(modeResultsTEEdge[0]**2*qTM)+(1/(modeResultsTECentre[0]**2*qTM)))))
    cTM = 2 * np.sqrt(omega0*epsilon_o / (np.abs(betaTM)*tEff))
    #eTE = np.zeros(pts*3)
    hFieldTM = np.zeros(pts*3)

    '''
    for i in range(0,len(cTE)):
        #region 1 E field
        eTE1 = cTE[i]*hTE[i]/qBarTE[i]*np.exp(qTE[i]*(x1+d/2))
        #region 2 E field
        eTE2 = cTE[i]*(hTE[i]/qBarTE[i]*np.cos(hTE[i]*(x2+d/2))+np.sin(hTE[i]*(x2+d/2)))

        #region 3 E field. We're assuming that the refractive indexes of the BOX and cladding are identical.
        eTE3 = cTE[i]*(hTE[i]/qBarTE[i]*np.cos(hTE[i]*d)+np.sin(hTE[i]*d))*np.exp(-qTE[i]*(x3-d/2))

        eTE += np.concatenate((eTE1,eTE2,eTE3),axis=None)

    '''

    for i in range(0,len(cTM)):
        #region 1 E field
        hTM1 = cTM[i]*hTM[i]/qBarTM[i]*np.exp(qTM[i]*(y1+w_rib/2))
        #region 2 E field
        hTM2 = cTM[i]*(hTM[i]/qBarTM[i]*np.cos(hTM[i]*(y2+w_rib/2))+np.sin(hTM[i]*(y2+w_rib/2)))

        #region 3 E field. We're assuming that the refractive indexes of the BOX and cladding are identical.
        hTM3 = cTM[i]*(hTM[i]/qBarTM[i]*np.cos(hTM[i]*w_rib)+np.sin(hTM[i]*w_rib))*np.exp(-qTM[i]*(y3-w_rib/2))

        hFieldTM += np.concatenate((hTM1,hTM2,hTM3),axis=None)
    
    ny=np.concatenate((modeResultsTEEdge[0]*np.ones(pts), modeResultsTECentre[0]*np.ones(pts), modeResultsTEEdge[0]*np.ones(pts)),axis=None)
    
    eFieldTM = hFieldTM / ny * eta
    x = np.concatenate((x1,x2,x3), axis=None)

    y = np.concatenate((y1,y2,y3), axis=None)
    #plt.plot(x,eTE)
    #plt.show()

    plt.plot(y,eFieldTM)
    plt.show()
    return [y,eFieldTM]


def get_neff_v(n_a,n_d,w_pn, w_rib,w_d, w_dn, w_dp, x_offset,d,h_slab, wavelength, n1, n2, v_a):

    #begin by getting 1D transverse field profile of our WG.
    [y,E] = get_field_profile(d,w_rib, h_slab,wavelength,n1,n2)

    #Obtain 1D transverse carrier profile of our WG.
    [n_profile0, p_profile0] = get_carrier_concentrations(n_a,n_d,0, w_pn, w_rib,w_d, w_dn, w_dp, x_offset)

    [n_profile, p_profile] = get_carrier_concentrations(n_a,n_d,v_a, w_pn, w_rib,w_d, w_dn, w_dp, x_offset)

    delNProfile = n_profile - n_profile0
    delPProfile = p_profile - p_profile0

    delRefrProfileN = -3.64 * 10**(-10) * wavelength**2 * delNProfile
    delRefrProfileP = -3.51 * 10**(-6) * wavelength**2 * delPProfile

    #compute overlap integrals...
    d_y = np.diff(y)

    del_neff_v = np.sum(E[:-1]**2 * delRefrProfileN[:-1] * d_y) / np.sum(E[:-1]**2 * d_y) + np.sum(E[:-1]**2 * delRefrProfileP[:-1] * d_y) / np.sum(E[:-1]**2 * d_y)

    print(del_neff_v)

    return

h_slab = 90e-9  
n1=1.444
n2=3.47
x_offset = 120e-9
v_a = -1
V_bi = get_vbi(1e17,1e17)
w_d = get_depletion_width(1e17,1e17,V_bi,v_a)
get_neff_v(1e17,1e17,1.2e-6, 5e-7,w_d, w_dn, w_dp, x_offset,220e-9,h_slab, 1.55e-6, n1, n2, v_a)



