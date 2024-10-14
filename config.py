##########################
### MODULES & PACKAGES ###
##########################

import numpy as np
import matplotlib.pyplot as plt

import scipy as sp
from scipy import constants 
import scipy.constants as const
from scipy.special import gamma as gamma_func
from scipy.integrate import quad, simpson

import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import z_at_value


#################
### CONSTANTS ###
#################

# Introducing useful constants from scipy
# Units [SI]

h = const.h  # Planck constant
hbar = const.hbar 
c = const.c  # Speed of light
pi = const.pi 
k = const.k  # Boltzmann constant
G = const.G  # Newton gravitational constant
M_s = 1.99 * 10**30  # Solar mass 
c_SB = const.Stefan_Boltzmann  # Stefan-Boltzmann constant
b_wien_freq = 0.05878925757646824946e12  # Wien's proportionality constant [Hz/K]


# Cosmological model
HUBBLE_CONST = 70
OMEGA_M = 0.3
cosmo = FlatLambdaCDM(H0=HUBBLE_CONST, Om0=OMEGA_M)


# Unit conversions
M_SOLAR_2_GRAMS = 1.989e33   # Multiply to make grams from Mo
E_JOULE_2_eV = 6.242e18      # Divide to make eV to Joules
YEARS_2_SEC = 3.15*1e7     # Multiply to make sec from year
KM_2_MPC = 3.24e-20 
HUBBLE_TIME = HUBBLE_CONST*KM_2_MPC*YEARS_2_SEC  # In 1/years
J_CRIT_2_SI = 1e-3           # Multiply the cgs value to get the SI value