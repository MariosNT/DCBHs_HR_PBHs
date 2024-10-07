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



def check_greybody_factor(M_bh, E_eV):
    """
    Function that checks if `greybody` factors need to
    be taken into account when calculating the spectrum
    of a BH
    """
    
    E_J = E_eV/E_JOULE_2_eV
    
    factor_check = h*c**3/(2*G*E_J)
    
    M_kg = M_bh*M_s
    
    bool_greybody = (M_kg >= factor_check)
    
    return bool_greybody



def age_at_redshift(z):
    """
    Function that calculates the age of the Universe in Gyrs
    for a given redshift
    
    Example: z=0 will give the age of the Universe today.
    """
    
    factor = (1/HUBBLE_TIME)*(2/(3*np.sqrt(1-OMEGA_M)))
    argument = (np.sqrt((1-OMEGA_M)*(1+z)**(-3))+np.sqrt((1-OMEGA_M)*(1+z)**(-3)+OMEGA_M))/np.sqrt(OMEGA_M)
    
    return factor*np.log(argument)/1e9



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



def Spectrum_freq_temperature(nu, Temp, f_eff = 1):
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



def Spectrum_freq_mass(nu, Mbh, units='solar', f_eff=1):
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
    return f_eff*2*h*nu**3/(c**2*exp_denominator)




def Greybody_spectrum_freq_mass(nu, Mbh, units='solar'):
    """
    ### WORK IN PROGRESS ###
    
    Blackbody spectrum with frequency for a given temperature
        
    Parameters
    ----------
    nu : frequency [Hz]
    
    Temp : temperature [K]
    
    Returns
    -------
    Blackbody spectrum
    """
    
    A_surface = 16*np.pi*G**2*Mbh**2/c**4
    
    Grey_factor = 4/9*A_surface/np.pi*Mbh**2*nu**4
    
    Temp = Hawking_temperature_from_mass(Mbh, units)
    
    exp_denominator = np.exp(h*nu/(k*Temp))-1
    return Grey_factor*2*h*nu**3/(c**2*exp_denominator)



def total_blackbody_spectrum(v, list_of_masses, units='solar', f_eff=1):
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
        total_spectrum += Spectrum_freq_mass(v, mass, units=units, f_eff=f_eff)
        
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
    """
    
    C_temperature = hbar*c**3/(8*pi*k*G)
    
    evaporation_constant = c_SB*hbar**4*c**6/(256*pi**3*k**4*G**2)*YEARS_2_SEC
    evaporation_constant /= 1e-9
    C_time = 3*evaporation_constant
    
    time_at_redshift = age_at_redshift(z)
    
    initial_mass = (C_time*time_at_redshift+(C_temperature/T_crit)**3)**(1/3)
    
    return initial_mass
    



################################
### PBHs evolution due to HR ###
################################ 


def evaporation_time_from_mass(M_bh):
    """
    Function that returns the evaporation time [in years] of a 
    black hole as a function of its mass [Mo].
    
    Parameters
    ----------
    M_bh : black hole mass [Mo]
    """
    
    evaporation_const = c_SB*hbar**4*c**8/(256*pi**3*k**4*G**2)
    return c**2*(M_s*M_bh)**3/(3*evaporation_const*YEARS_2_SEC)



def effective_evaporation_time_from_mass(M_bh):
    """
    Function that returns the evaporation time [in years] of a 
    black hole as a function of its mass [Mo].
    
    Parameters
    ----------
    M_bh : black hole mass [Mo]
    """
    
    if M_bh*M_s <= 1e18:
        c1 = -0.3015
        c2 = 0.3113
        p = -8e-4
        a_eff = c1 + c2*(M_bh*M_s)**p
    else:
        # a_eff = 2.011*1e-4
        a_eff = 1/(15360*np.pi)
    
    evaporation_const = G**2/(3*a_eff*c**4*hbar)
    return evaporation_const*(M_s*M_bh)**3/YEARS_2_SEC



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
    
    evaporation_constant = c_SB*hbar**4*c**6/(256*pi**3*k**4*G**2)*YEARS_2_SEC
    evaporation_constant /= 1e-9  # Transform to 1/Gyrs
       
    if time_units == 'redshift':
        t_init = age_at_redshift(t0)
        t_final = age_at_redshift(t)
        delta_t = t_final - t_init
    elif time_units == 'Gyrs':
        delta_t = t - t0
    else:
        print("Wrong time units!")
        
    if mass_units == 'solar':
        M_initial = M_bh*M_s
    elif mass_units == 'cgs':
        M_initial = M_bh/1000
    else:
        print("Wrong mass units!")
            
    M_evolved_cubed = M_initial**3 - 3*evaporation_constant*delta_t
    if M_evolved_cubed < 0:
        M_evolved = 0
    else:
        M_evolved = M_evolved_cubed**(1/3)
    
    return M_evolved/M_s



###########################
### PBHs Mass Functions ###
###########################


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