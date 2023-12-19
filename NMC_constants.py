""""
Physical and cosmological constants
Author: Arta Khosravi
Started: Jul. 2022
Last update: Dec. 2023
"""
import numpy as np

# Misc
c = 300000 #km/s
#-------#
# It's worth noting that because of multiple fiducial models in various datasets which we investigated,
# We needed to check a different H0 for example BAO than SNIa'14. 
#-------#
H0=68.26 #planck no CMB lensing #100h * km s-1 MPC-1 #Our main choice for the SNIa
OB=0.02242/(0.6826)**2 #Baryonic matter density
#-------#
# H0=73.8-2.4 #riess2014 #or +2.4
# H0=73.48 #R18
# H0=74.03 #riess2019
# H0=68.14 #WMAP-BAO
# H0=65.2 #riess1998
# H0 = (67.32) #plancklcdm
#-------#Our main choice for the BAO#-------#
# H0=67.66 
# OB=0.02242/(0.6766)**2
# h=0.6766
#-------#Other Constants such as densities:#-------#

# Cosmo
OR= 2.4e-5 #Radiation density
a0 = 1 #Scale factor at t=t0
OMP =  0.3111 #Cold dark matter density according to Planck'18
OKP=0.001 #Curvature density according to Planck'18
OLP=0.6889 #Dark energy density according to Planck'18
OLdw= 0.6847 - 0.0073 #Dark energy density with negative error according to Planck'18
OLuw= 0.6847 + 0.0073 #Dark energy density with positive error according to Planck'18
G=6.67*1e-20 #km3 kg-1 s-2 #Newton's gravity constant
rho_cr=3*(H0**2)/(8*np.pi*G) #Mpc-2 kg km-1 #critical radius
zdrag= 1060.01 #redshift at the drag epoch
rdrag=147.21 #radius at the drag epoch
Sigmav=0.17 #Sigma_\nu defined in Chi Square for peculiar velocities