#Author: Alexander Greenwood
#The following scripts have been used for my BASc thesis on a "Segmented Microring Modulator"


#Based on the following papers:
#Analytical Modeling of Silicon Microring and Microdisk Modulators With Electrical and Optical Dynamics
#By Dubé-Demers et al.
#Ultrafast pulse-amplitude modulation with a femtojoule silicon photonic modulator
#By Dubé-Demers et al.
#Analytical Model and Fringing-Field Parasitics of Carrier-Depletion Silicon-on-Insulator Optical Modulation Diodes
#By Jayatilleka et al.


#This is a work in progress :)

import matplotlib.pyplot as mpl
import math
import numpy as np
import scipy as scipy
from scipy import special
#First define some physical constants
n_i = 1.1e10 #intrinsic carrier concentration [cm^-3]
k_b = 1.38*10**(-23) #Boltzmann's constant [JK^-1]
epsilon_0 = 8.85*10**(-14) #permittivity of free space [F/cm] 
epsilon_s = 11.7 # relative permittivity of silicon
epsilon_sio =3.9
q = 1.602*10**(-19) # charge of an electron [C]
temp = 300 # [K]
pi = 3.14

mu_n = 364
mu_p = 205

#d_n = k*T / q * mu_n
#d_p = k*T / q * mu_p

#l_p = np.sqrt(d_p*tau_p) # diffusion length for holes
#l_n = np.sqrt(d_n*tau_n) # diffusion length for electrons


NPOINTS = 1500
#############################################################

def get_carrier_concentrations(N_A,N_D,V,w_slab,w_rib,w_d, w_dn, w_dp, x_offset):
    p_0n = n_i**2 / N_D #minority carriers on n side.
    n_0p = n_i**2 / N_A #minority carriers on p side.

    x = np.linspace(-w_slab/2,w_slab/2,NPOINTS)
    
    iwdp = (abs(x + w_dp-x_offset)).argmin()
    iwdn = (abs(x - w_dn-x_offset)).argmin()
    
    w_dp *= -1
    w_dp += x_offset

    w_dn += x_offset
    xp = x[0:iwdp]
    xn = x[iwdn:len(x)]
    #Begin with minority carrier concentration in p region.
    n_xp = n_0p*((np.exp(q*V/(k_b*temp))-1) * (1-abs((w_dp-xp)/(w_dp - w_slab))))
    p_xn = p_0n*((np.exp(q*V/(k_b*temp))-1) * (1-abs((w_dn-xn)/(w_slab - w_dn))))

    #Depletion region shall be empty.
    n_xd = np.zeros(iwdn-iwdp)
    p_xd = n_xd

    n_xn = np.full(NPOINTS - iwdn,N_D) # is this correct?????
    p_xp = np.full(iwdp,N_A) 

    n_profile = np.concatenate((n_xp,n_xd,n_xn), axis=None)

    p_profile = np.concatenate((p_xp,p_xd,p_xn), axis=None)



    #carrier_profile = np.append(n_0_profile,p_0_profile)
    
    #For long base assumption, we are assuming that diffusion length is around same size or smaller
    #than base length.
    #Use einstein's relation for obtaining diffusion const
    return [n_profile,p_profile]
    
    

def get_depletion_width(N_A,N_D,V_bi,V):
    inner = ((2*epsilon_0 *epsilon_s)/q)*((N_A+N_D)/(N_A*N_D))*(V_bi-V-(2*k_b*temp)/q)
    w_d = np.sqrt(inner)

    w_dp = w_d / (1 + (N_A/N_D))
    w_dn = w_d / (1 + (N_D/N_A))
    
    inner = ((2*epsilon_0 *epsilon_sio)/q)*((N_A+N_D)/(N_A*N_D))*(V_bi-V-(2*k_b*temp)/q)
    w_d_tb = np.sqrt(inner)

    w_dp_tb = w_d_tb / (1 + (N_A/N_D))
    w_dn_tb = w_d_tb / (1 + (N_D/N_A))


    #Express as m instead of cm.
    w_d *=0.01
    w_dp *=0.01
    w_dn *=0.01
    w_dp_tb *=0.01
    w_dn_tb *= 0.01

    return [w_d,w_dp,w_dn,w_dp_tb,w_dn_tb]

def get_vbi(N_A,N_D):
    p_0 = n_i**2 / N_A
    n_0 = n_i**2 / N_D
    V_bi = k_b*temp / (2*q) * math.log(N_A*N_D / (p_0*n_0))
    return V_bi

def get_c_parallel(h_rib,w_d):
    #usually h_rib is just 200 nm
    c_parallel = (epsilon_0*epsilon_s*h_rib)/w_d
    print(c_parallel)
    return c_parallel

def get_c_fringe(h_rib,w_d):
    c_fringe = epsilon_0*(2*epsilon_sio / (2*np.pi))*math.log(2*np.pi * h_rib / w_d)
    print(c_fringe)
    return c_fringe
def get_c_vertical(w_d,t_p,t_n):
    k = np.sqrt(w_d*(w_d + t_p)*(w_d + t_n + t_p)/ (w_d + t_n))
    m = k**2
    k_prime = np.sqrt(1 - k**2)
    m_prime = k_prime**2

    c_top = epsilon_0 * epsilon_sio * scipy.special.ellipk(m_prime)/ scipy.special.ellipk(m)

    return c_top

def get_t(w_dp,w_dn,w_dp_tb,w_dn_tb,w_rib):
    #Given N_a |x_p| = N_d |x_n
    t_n = 2*w_dn_tb / np.pi
    t_p = 2*w_dp_tb / np.pi
   
    return [t_p,t_n]

def get_xpn(w_dp,w_dn,w_dp_tb,w_dn_tb,w_rib):
    v = np.linspace(-np.pi/2,0,NPOINTS)
    x_p = (2*w_dp_tb / np.pi) + w_dp - (4*w_dp_tb / np.pi * np.sin(v/2)**2)
    x_n = (4*w_dn_tb / np.pi * np.sin(v/2)**2) - ((2*w_dn_tb / np.pi) + w_dn)
    y = (2*w_dn_tb / np.pi)*(np.log(np.tan(np.pi/4 + v/2))-np.sin(v))
    
    
    return [x_p,x_n]

def get_r_np(w_pn,w_rib,N_A,N_D,h_slab,L,w_dn,w_dp,x_offset,h_rib):
    w_pn *=100
    w_rib *=100
    L *=100
    x_offset *=100
    w_dp *=100
    w_dn *=100
    h_rib *=100

    r_p = (w_pn - w_rib)/(2*q*N_A*mu_p*h_slab*L) + (w_rib/2 + x_offset - w_dp)/(q*N_A*mu_p*h_rib*L)
    r_n = (w_pn - w_rib)/(2*q*N_D*mu_n*h_slab*L) + (w_rib/2 - x_offset - w_dn)/(q*N_D*mu_n*h_rib*L)
    
    print("L = " + str(L))
    print("Rp = " + str(r_p))
    print("Rn= " + str(r_n))
    return [r_p, r_n]

############################################################
'''
V_a = -1
L = 0.75 * 2 * np.pi * 5e-6 # assuming radius is 5 um.
x_offset = 120e-9
w_pn = 1.2e-6
h_slab = 90e-9  
w_rib = 5e-7
h_rib = 200e-9
#Values below taken from Analytical Model and Fringing-Field Parasitics of Carrier-Depletion Silicon-on-Insulator Optical Modulation Diodes
N_a = 5e17
N_d = 5e17
#get_carrier_concentrations(N_a,N_d,-1,1.5e-6)
V_bi = get_vbi(N_a,N_d)
w_d = get_depletion_width(N_a,N_d,V_bi,V_a)

w_dp = w_d[1]
w_dn = w_d[2]
w_dp_tb = w_d[3]
w_dn_tb = w_d[4]
w_d = w_d[0]

t= get_t(w_dp,w_dn,w_dp_tb,w_dn_tb,w_rib)
t_p = t[0]
t_n = t[1]
c_vertical = get_c_vertical(w_d,t_p,t_n)
c_p = get_c_parallel(h_rib,w_d)
c_f = get_c_fringe(h_rib,w_d)
r_np = get_r_np(w_pn,w_rib,N_a,N_d,h_slab,L,w_dn,w_dp,x_offset,h_rib)
#bw = (2*np.pi*(r_np[0]+r_np[1])*(c_vertical+c_f+c_p)*L)**(-1)
A = (q/w_pn)*(mu_p*mu_n*N_a*N_d/(mu_p*N_a + mu_n*N_d))
bw = A/(2*np.pi)*(h_slab/h_rib)*(epsilon_0*epsilon_s / w_d + (c_f + c_vertical)/h_rib)**(-1)
print(str(bw/(1e9)) + " GHz")

get_carrier_concentrations(N_a,N_d,V_a,w_pn, w_rib,w_d, w_dn, w_dp, 0)
'''