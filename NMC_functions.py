
"""
Cosmological functions
Author: Arta Khosravi
Started: Jul. 2022
Last update: Dec. 2023
"""

import numpy as np
import scipy as sp
from scipy.integrate import quad
import scipy.integrate as integrate
import math
from Source_NMC.NMC_constants import *

#-----#For SNIa:#----#
def H2(z,om,omneg,ol,w): # Hubble function, s^-2
    return (H0**2)*(((om)*((1 + z)**3))+(omneg)*((1 + z)**(3*(1+w)))+ (ol)+(1-om-ol-omneg)*((1+z)**2))
def Xlintegrat(z,om,omneg,ol,w): #Comoving Distance to the object
    res=integrate.quad((lambda x : 1/(np.sqrt(H2(x,om,omneg,ol,w)))),0,z)
    return res[0]
def Skneg(z,om,omneg,ol,w): #S_k(chi) for Omega_curvature > 0
    return (1/(np.sqrt(abs(1-om-ol-omneg))))*(c/H0)*np.sinh(np.sqrt(abs(1-om-ol-omneg))*H0*Xlintegrat(z,om,omneg,ol,w))
def Skpos(z,om,omneg,ol,w): #S_k(chi) for Omega_curvature < 0
    return (1/(np.sqrt(abs(1-om-ol-omneg))))*(c/H0)*np.sin(np.sqrt(abs(1-om-ol-omneg))*H0*Xlintegrat(z,om,omneg,ol,w))
def Skflat(z,om,omneg,ol,w): #S_k(chi) for Omega_curvature = 0
    return (c)*Xlintegrat(z,om,omneg,ol,w)
def DL(z,om,omneg,ol,w): #Luminosity distance for Omega_curvature < 0 Or Close Universe #k=+1: sin
    return Skpos(z,om,omneg,ol,w)*(1+z)
def DLneg(z,om,omneg,ol,w): #Luminosity distance for Omega_curvature > 0 Or Open Universe #k=-1: inh 
    return Skneg(z,om,omneg,ol,w)*(1+z)
def DLflat(z,om,omneg,ol,w): #Luminosity distance for Omega_curvature = 0 Or Flat Universe #k=0 
    return Skflat(z,om,omneg,ol,w)*(1+z)
def Mp(dl): #Distance Modulus or m-M
    return (5*((np.log10(dl))))+25
def X2(mp,m0,sigma0): #Chi Square using distance modulus
    return (((mp-m0)**2)/(((sigma0)**2)+((Sigmav)**2)))

#----For BAO#----#

def Dv(z,om,omneg,ol,w): 
    return (((Skpos(z,om,omneg,ol,w))**2)*c*z/np.sqrt(H2(z,om,omneg,ol,w)))**(1/3)
def Dvflat(z,om,omneg,ol,w):
    return (((Skflat(z,om,omneg,ol,w))**2)*c*z/np.sqrt(H2(z,om,omneg,ol,w)))**(1/3)
def Dvneg(z,om,omneg,ol,w):
    return (((Skneg(z,om,omneg,ol,w))**2)*c*z/np.sqrt(H2(z,om,omneg,ol,w)))**(1/3)
def rhob(z):
    return OB*rho_cr*((1+z)**(3))
def rhoph(z):
    return OR*rho_cr*((1+z)**(4))
def cs(z): #https://arxiv.org/pdf/1411.1074.pdf
    return np.power(3,-1/2)*c*np.power(1+(3/4)*(rhob(z)/rhoph(z)),(-1/2))
def rdint(z,om,omneg,ol,w):
    res=integrate.quad((lambda x : cs(x)/(np.sqrt(H2(x,om,omneg,ol,w)))),z,9607)
    return res[0]
def Dvrdp(z,om,omneg,ol,w):
    return Dv(z,om,omneg,ol,w)/rdint(zdrag,om,omneg,ol,w)
def Dvrdflat(z,om,omneg,ol,w):
    return Dvflat(z,om,omneg,ol,w)/rdint(zdrag,om,omneg,ol,w)
def Dvrdneg(z,om,omneg,ol,w):
    return Dvneg(z,om,omneg,ol,w)/rdint(zdrag,om,omneg,ol,w)
def X2dv(Dvtrd,Dvobtrd,sigma0trd):
    return ((Dvtrd-Dvobtrd)**2)/((sigma0trd)**2)

#----General----#
def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))