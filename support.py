from config import *


##########################################
# Useful functions for the main notebook #
##########################################

def normalise(lista):
    """ 
    Function that normalises a list of data 
    """
    
    return lista/np.sum(lista)



def Schwarzschild_radius(M_bh):
    """
    Function that returns the Schwarzschild radius
    of a BH in [m] for a M_bh in [Mo]
    """
    
    M_kg = M_bh*M_s
    
    return 2*G*M_kg/c**2



def theoretical_lum(Radius, Temp):
    """
    This function returns the luminosity of a black 
    hole as predicted by the Stefan-Boltzmann Law.
    
    Parameters
    ----------
    Radius : Object radius [m]
    
    Temp: Temperature [K]
    """
    
    return 4*pi*Radius**2*c_SB*Temp**4



##############################
### Cosmological functions ###
##############################


def age_at_redshift(z):
    """
    Function that calculates the age of the Universe in Gyrs
    for a given redshift
    
    Example: z=0 will give the age of the Universe today.
    """
    
    factor = (1/HUBBLE_TIME)*(2/(3*np.sqrt(1-OMEGA_M)))
    argument = (np.sqrt((1-OMEGA_M)*(1+z)**(-3))+np.sqrt((1-OMEGA_M)*(1+z)**(-3)+OMEGA_M))/np.sqrt(OMEGA_M)
    
    return factor*np.log(argument)/1e9


def H(z):
    """
    Hubble function
    """
    return np.sqrt(HUBBLE_CONST**2*(OMEGA_M*(1+z)**3+(1-OMEGA_M)))


def Omega_m_redshift(z):
    num = OMEGA_M*(1+z)**3
    denom = OMEGA_M*(1+z)**3+(1-OMEGA_M)
    
    return num/denom


def Delta_c(z):
    """
    Matter overdensity in haloes
    """
    
    y = Omega_m_redshift(z)-1
    
    return 18*pi**2+82*y-39*y**2


def virial_radius(z, Mvir):
    """
    Calculate virial radius [in units of m].
    
    Mvir = halo mass [in Mo]
    """
    
    Mvir_SI = Mvir*M_s
    
    num = 2*G*Mvir_SI
    
    denom = Delta_c(z)*(H(z)*KM_2_MPC)**2  # To make it in SI units
    
    rvir = (num/denom)**(1/3)
    
    return rvir


def virial_density(z):
    """
    Calculate virial density [in units of kg/m^3].
    """    
    
    rho_crit = 2*(H(z)*KM_2_MPC)**2/(8*pi*G)
    
    return rho_crit*Delta_c(z)



#######################################
### Hawking Temperature and Spectra ###
####################################### 

def Hawking_temperature_from_mass(M_bh, units='solar'):
    """
    Returns the Hawking temperature in K,
    for a BH mass in Mo
    """
    
    if units == 'solar':
        T_HR = hbar*c**3/(8*pi*k*G*M_bh*M_s)
    elif units == 'cgs':
        # For BH masses provided in grams
        T_HR = hbar*c**3/(8*pi*k*G*(M_bh/M_SOLAR_2_GRAMS)*M_s)
    
    return T_HR


def BH_mass_from_Hawking_temperature(T_HR, units='solar'):
    """
    Returns BH mass [in Mo, or gr] for Hawking temperature in K.
    """
    
    if units == 'solar':
        M_bh = hbar*c**3/(8*pi*k*G*T_HR*M_s)
    elif units == 'cgs':
        # For BH masses in grams
        M_bh = hbar*c**3/(8*pi*k*G*(1/M_SOLAR_2_GRAMS)*M_s*T_HR)
    
    return M_bh



def Spectrum_freq_temperature(nu, Temp, f_eff = 0.2):
    """
    Blackbody spectrum with frequency for a given temperature
    
    
    Parameters
    ----------
    nu : frequency [Hz]
    
    Temp : temperature [K]
    
    Returns
    -------
    Blackbody spectrum
    """
    
    exp_denominator = np.exp(h*nu/(k*Temp))-1
    return f_eff*2*h*nu**3/(c**2*exp_denominator)



def Spectrum_freq_mass(nu, Mbh, units='solar', f_eff=0.2, f_grey=0.24):
    """
    Blackbody spectrum with frequency for a given BH mass.
    
    Parameters
    ----------
    nu : frequency [Hz]
    
    Mbh : BH mass in solar or cgs units
    
    Returns
    -------
    Blackbody spectrum
    """
    
    Temp = Hawking_temperature_from_mass(Mbh, units)
    
    exp_denominator = np.exp(h*nu/(k*Temp))-1
    return f_grey*f_eff*2*h*nu**3/(c**2*exp_denominator)




def total_blackbody_spectrum(v, list_of_masses, units='solar', f_eff=0.2, f_grey=0.24):
    """
    This function takes a list of masses and calculates the total HR blackbody spectrum.
    
    Parameters:
    -----------
    
    ν: Frequency [Hz]
    
    list_of_masses [Solar masses]
    """
    
    ## First, a "holder variable" is called so we can collect all the contributions
    total_spectrum = 0
    
    ## Then, we need to loop through all the masses in list_of_masses
    for mass in list_of_masses:
        # To calculate each the blackbody spectrum of each individual mass, we can just 
        # call the function mass_blackbody_spectrum (we which defined previously).
        total_spectrum += Spectrum_freq_mass(v, mass, units=units, f_eff=f_eff, f_grey=f_grey)
        
    return total_spectrum



def energy_from_photons(nu = None, wavelength=None):
    """
    Function that returns energy in [eV] for provided frequency or wavelength.
    
    Parameters
    ----------
    nu : frequency [Hz]
    
    wavelegth : wavelength [nm]
    """
    
    if nu is None and wavelength is None:
        raise ValueError("Either `frequency` or `wavelength` must be specified")
          
    # Calculate frequency if wavelegth is provided
    nu = c/(wavelength*1e-9) if nu is None else nu
    
    # Calculate energy from frequency
    Eg = h*nu
    
    return Eg*E_JOULE_2_eV



def photons_from_energy(E_eV):
    """
    Function that returns frequency [Hz] and wavelength [nm], 
    for given energy in [eV].
    
    Parameters
    ----------
    E_eV : photon energy [eV]
    """
    
    nu = (E_eV/E_JOULE_2_eV)/h
    
    wavelength = c/nu*1e9  # Transform to [nm]
          
    return nu, wavelength



def peak_HR_energy_from_mass(Mbh):
    """
    Function that return peak energy [in eV] of blackbody spectrum, 
    for a given mass in [Mo]
    
    Parameters
    ----------
    Mbh : black hole mass [Mo]
    """
    
    T_HR = Hawking_temperature_from_mass(Mbh)
    
    # Wien's law (different b when using frequency/wavelength)
    freq_peak = b_wien_freq*T_HR
    
    energy_peak = energy_from_photons(nu=freq_peak)
    
    return energy_peak



def mass_from_bb_spectrum(nu, Bv):
    """
    Function that returns the mass of a BH [in Mo], based on its emission 
    of Bv at frequency v
    """
    
    Ao = 2*h*nu**3/c**2
    factor = c**3/(16*np.pi**2*G*nu)
    
    mass = factor*np.log(1+Ao/Bv)
    
    return mass/M_s


def initial_mass_with_Tcrit_and_z(z, T_crit):
    """
    Function that returns the initial mass value that would have
    a critical HR temperature at a given redshift.
    
    Parameters
    ----------
    T_crit : critical temperature
    
    z : redshift of interest
    
    Returns
    -------
    initial_mass : initial BH mass [in Mo]    
    """
    
    ## Temperature evolution 
    C_temperature = hbar*c**3/(8*pi*k*G)
    ## This has units of mass, so we need to rescale
    T_ratio = C_temperature/T_crit/M_s
    
    ## Mass evolution
    time_at_redshift = age_at_redshift(z)

    ## C in notes, in [kg^3/s]
    evaporation_const = 3*c_SB*hbar**4*c**6/(256*pi**3*k**4*G**2)
    ## transform to [Mo^3/Gyr]
    evaporation_const *= YEARS_2_SEC*1e9/M_s**3
    
    ## Time evolution
    time_factor = time_at_redshift*evaporation_const
    
    ## Initial mass at [Mo]
    initial_mass = (time_factor + T_ratio**3)**(1/3)
    
    return initial_mass
    



################################
### PBHs evolution due to HR ###
################################ 


def evaporation_time_from_mass(M_bh):
    """
    Function that returns the evaporation time [in Gyears] of a 
    black hole as a function of its mass [Mo].
    
    Parameters
    ----------
    M_bh : black hole mass [Mo]
    """
    
    evaporation_const = c_SB*hbar**4*c**8/(256*pi**3*k**4*G**2)
    return c**2*(M_s*M_bh)**3/(3*evaporation_const*YEARS_2_SEC)/1e9


def evaporation_redshift_from_mass(M_bh):
    """
    Function that returns the evaporation redshift of a 
    black hole as a function of its mass [Mo].
    
    Parameters
    ----------
    M_bh : black hole mass [Mo]
    """
    
    evaporation_const = c_SB*hbar**4*c**8/(256*pi**3*k**4*G**2)
    
    ## Evaporation time [Gyr]
    t_evap = c**2*(M_s*M_bh)**3/(3*evaporation_const*YEARS_2_SEC)/1e9
    
    ## Get redshift, assuming cosmology
    z_evap = z_at_value(cosmo.age, t_evap * u.Gyr, zmin=0, zmax=1e8).value
    
    return z_evap


def mass_evolution_half_time(M_bh):
    """
    Function that returns the half-time time [in Gyears] for
    a BH of initial mass Mi [Mo] to evolve to Mi.
    
    Parameters
    ----------
    M_bh : black hole mass [Mo]
    """
    
    ## C in notes, in [kg^3/s]
    evaporation_const = 3*c_SB*hbar**4*c**6/(256*pi**3*k**4*G**2)
    
    ## transform to [Mo^3/Gyr]
    evaporation_const *= YEARS_2_SEC*1e9/M_s**3
    
    t_half_const = 7/8/evaporation_const
    
    return t_half_const*M_bh**3


def mass_from_evaporation_time(t_evap):
    """
    Function that takes the evaporation time [in Gyears] and  
    returns the black hole mass [Mo] that would evaporate by then.
    
    Parameters
    ----------
    t_evap : evaporation time [Gyr]
    """
    
    evaporation_const = c_SB*hbar**4*c**8/(256*pi**3*k**4*G**2)
    
    M_bh = (t_evap/(c**2/(3*evaporation_const*YEARS_2_SEC)/1e9))**(1/3)/M_s
    
    return M_bh



def mass_from_evaporation_redshift(z_evap):
    """
    Function that takes the evaporation redshift and  
    returns the black hole mass [Mo] that would evaporate by then.
    
    Parameters
    ----------
    t_evap : evaporation time [Gyr]
    """
    
    evaporation_const = c_SB*hbar**4*c**8/(256*pi**3*k**4*G**2)
    t_evap = age_at_redshift(z_evap)
    
    M_bh = (t_evap/(c**2/(3*evaporation_const*YEARS_2_SEC)/1e9))**(1/3)/M_s
    
    return M_bh


def effective_evaporation_time_from_mass(M_bh):
    """
    Function that returns the evaporation time [in Gyears] of a 
    black hole as a function of its mass [Mo].
    
    Parameters
    ----------
    M_bh : black hole mass [Mo]
    """
    
    ## If BH mass in grams is smaller than 1e18
    if M_bh*M_s <= 1e18:
        c1 = -0.3015
        c2 = 0.3113
        p = -8e-4
        a_eff = c1 + c2*(M_bh*M_s)**p
    else:
        # a_eff = 2.011*1e-4
        a_eff = 1/(15360*np.pi)
    
    evaporation_const = G**2/(3*a_eff*c**4*hbar)
    return evaporation_const*(M_s*M_bh)**3/YEARS_2_SEC/1e9



def mass_evolution_from_HR(M_bh, t0, t, mass_units='solar', time_units='redshift'):
    """
    Function that returns the evolved BH mass due to HR, 
    given the initial mass and time difference.
    
    
    Parameters
    ----------
    M_bh : initial black hole mass
    
    t0 : initial time parameter (time in units of Gyrs)
    
    t : final time parameter (time in units of Gyrs)
    
    Returns
    -------
    M_evolved : evolved BH mass [in Mo]
    """
    
    ## Evaporation constant in SI
    evaporation_constant = 3*c_SB*hbar**4*c**6/(256*pi**3*k**4*G**2)
       
    ## Time difference in SI
    if time_units == 'redshift':
        t_init = age_at_redshift(t0)
        t_final = age_at_redshift(t)
        delta_t = (t_final - t_init)*1e9*YEARS_2_SEC
    elif time_units == 'Gyrs':
        delta_t = (t - t0)*1e9*YEARS_2_SEC
    else:
        print("Wrong time units!")
        
    ## Masses in SI
    if mass_units == 'solar':
        M_initial = M_bh*M_s
    elif mass_units == 'cgs':
        M_initial = M_bh/1000
    else:
        print("Wrong mass units!")
            
    M_evolved_cubed = M_initial**3 - evaporation_constant*delta_t
    
    if M_evolved_cubed < 0:
        M_evolved = 0
    else:
        M_evolved = M_evolved_cubed**(1/3)
    
    ## Return evolved BH mass in Mo
    return M_evolved/M_s



########################################
### PBHs Mass Functions  & Formation ###
########################################


def formation_mass_at_time(t_form):
    """
    Function that gives an estimate of the PBH mass [Mo]
    at formation time (at specific t_age).
    
    Parameters
    ----------
    t_form : time at formation (Gyr)
    
    Returns
    -------
    M_bh : PBH mass at formation (Mo)
    """
    
    factor = 1e15
    time_scale = t_form*1e9*YEARS_2_SEC/1e-23
    
    M_bh = factor*time_scale  # Mass in gr
    
    M_bh /= M_SOLAR_2_GRAMS
    
    return M_bh




def lognormal_PBH_mass_function(mass, Mc, sigma):
    """
    This PBH mass function is most applicable when the PBHs in 
    question "originate from a smooth, symmetric peak in the 
    inflationary power spectrum under the slow-roll approximation". 
    (Chen and Hall, 2023, p. 5)
    
    Parameters
    ----------    
    mass : PBH mass in [Mo]
    
    Mc : Median mass [Mo]
    
    sigma : Mass distribution width
    """
        
    denominator = np.sqrt(2*np.pi)*sigma*mass
    exponent = - np.log(mass/Mc)**2/(2*sigma**2)
    return np.exp(exponent)/denominator



def power_law_PBH_mass_function(mass, M_min, alpha):
    """
    This PBH mass function emerges as a result of "a broad or 
    flat power spectrum of curvature perturbations during the 
    radiation-dominated era". (Chen and Hall, 2023, p. 5) 
    
    Parameters
    ----------
    mass : PBH mass in [Mo]
    
    M_min : Lower PBH mass limit [Mo]
    
    alpha : Power law index (α > 1)
    """
    
    bool_masses = mass>M_min
    pbh_mass_function = ((alpha-1)/M_min)*(mass/M_min)**-alpha
    
    pbh_mass_function[~bool_masses]=0
    
    return pbh_mass_function



def broken_pl_PBH_mass_function(mass, m_peak, alpha_1, alpha_2):
    """
    This PBH mass function arrises due to the theory that some
    PBHs are "formed by vacuum bubbles nucleating during inflation
    via quantum tunneling". (Chen and Hall, 2023, p. 5)

    Parameters
    ----------
    m : PBH mass in [Mo]
    
    m_peak : Peak mass of mP(m), where P(m) is the mass function [Mo]
    
    alpha_1: First power law index (α_1 > 0)
    
    alpha_2: Second power law index (α_2 > 1)
    """

    coefficient = (m_peak/(alpha_1+1) + m_peak/(alpha_2-1))**-1
    factor_1 = (mass/m_peak)**alpha_1
    factor_2 = (mass/m_peak)**-alpha_2

    mass_function = coefficient*factor_2

    # To create the piecewise break
    bool_threshold = mass < m_peak
    mass_function[bool_threshold] = coefficient*factor_1[bool_threshold]

    return mass_function



def crit_collapse_PBH_mass_function(mass, M_f, alpha):
    """
    This PBH mass function does not contain a lower mass limit. 
    However, M_f acts like an upper mass limit (in that the mass
    function past M_f exhibits an exponential-like decay).
    
    Parameters
    ----------
    mass : PBH mass [Mo]
    
    M_f: Horizon mass at the epoch of collapse [Mo]
    
    alpha: Index related to critical collapse of radiation
    """
    
    fraction = (alpha**2*mass**alpha)/(M_f**(1+alpha)*gamma_func(1/alpha))
    exponent = -(mass/M_f)**alpha
    
    return fraction*np.exp(exponent)



###########################################
### LW radiation and critical intensity ###
###########################################


def J_LW_single_BH(distance, Mbh, E_eV=12.5, units='solar', f_eff=0.2, f_grey=0.24):
    """
    LW specific intensity for a single BH at distance d.
    
    Parameters
    ----------
    distance : distance from the BH [pc]
    
    Mbh : BH mass [Mo]
    
    E_eV : energy of the photons of interest [eV]
    
    Returns
    -------
    J_LW_21_cgs : specific intensity in normalised units (J_21) in cgs
    """
    
    B_LW = Spectrum_freq_mass(photons_from_energy(E_eV)[0], Mbh, units, f_eff, f_grey)
    R_S = Schwarzschild_radius(Mbh)
    
    J_LW_21 = B_LW/4*(R_S/(distance*PARSEC_2_M))**2*1e21  # This is in SI units
    
    # Transforming to cgs
    J_LW_21_cgs = J_LW_21/J_CRIT_2_SI
    
    return J_LW_21_cgs



def J_LW_BH_density_factor(Mbh, rho, E_eV=12.5, units='solar', f_eff=0.2, f_grey=0.24):
    """
    LW specific intensity factor for the cases with a BH density.
    See sections 4.3.2 & 4.3.3 of the paper.
    
    Parameters
    ----------
    rho : density of the halo (at rvir) [kg/m^3]
    
    Mbh : BH mass [Mo]
    
    E_eV : energy of the photons of interest [eV]
    
    Returns
    -------
    J_LW_21_factor_cgs : J_LW factor (not units of specific intensity) in normalised units (J_21) in SI
    """    
    
    B_LW = Spectrum_freq_mass(photons_from_energy(E_eV)[0], Mbh, units, f_eff, f_grey)
    R_S = Schwarzschild_radius(Mbh)
    
    # Calculating the numerator
    J_LW_21_factor = B_LW*R_S**2*1e21
    
    # adding the constribution from density and mass
    J_LW_21_factor *= rho/(Mbh*M_s)
       
    return J_LW_21_factor


def J_LW_critical_density(Mass_halo, Mbh, z, rmin, factor=1, method='isothermal'):
    """
    Toral specific intensity for 
    
    Parameters
    ----------
    Mass_halo : Halo mass [Mo]
       
    Mbh : BH mass [Mo]
    
    z : redshift
    
    rmin : minimum radius of PBHs distribution [in pc]
    
    Returns
    -------
    J_21 : Critical intensity at the centre of the halo, in normalised units
    """    
    
    # Virial radius
    rvir = virial_radius(z, Mass_halo)/PARSEC_2_M  #to make radius in units of parsec
    
    # Virial density - density of the halo (at rvir) [kg/m^3]
    rho_vir = virial_density(z)*factor
    
    # Base J21 factor
    J_base = J_LW_BH_density_factor(Mbh, rho_vir)
    
    if method == 'isothermal':
        J_21 = J_base*(rvir*PARSEC_2_M)*(rvir/rmin-1)/J_CRIT_2_SI  # To transform to cgs
    elif method == 'uniform':
        J_21 = J_base*(rmin*PARSEC_2_M)*(rvir/rmin-1)/J_CRIT_2_SI
        
    return J_21


def LW_Bv_RS(Mbh, E_eV=12.5, units='solar', f_eff=0.2, f_grey=0.24):
    """
    BvRs factor for LW specific intensity for a single BH.
    
    Parameters
    ----------    
    Mbh : BH mass [Mo]
    
    E_eV : energy of the photons of interest [eV]
    
    Returns
    -------
    BvRs : the specific intensity x Schwarzschild radius factor [SI units]
    """
    
    B_LW = Spectrum_freq_mass(photons_from_energy(E_eV)[0], Mbh, units, f_eff, f_grey)
    R_S = Schwarzschild_radius(Mbh)
    
    return B_LW*R_S**2


def LW_background(Mbh, z, E_eV=12.5, units='solar', f_eff=0.2, f_grey=0.24, f_PBH=1, f_z=1.04):
    
    BH_term = c*LW_Bv_RS(Mbh, E_eV, units, f_eff, f_grey)/(Mbh*M_SOLAR_2_GRAMS)
    
    J_crit_norm = 1e21/J_CRIT_2_SI
    
    cosmo_term = 3*HUBBLE_CONST**2/(8*G)*f_PBH*OMEGA_M*(1+z)**3*(f_z-1)/H(z)
    cosmo_term *= KM_2_MPC
    
    return BH_term*cosmo_term*J_crit_norm