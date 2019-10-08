import numpy as np
from yt import YTQuantity
from yt import YTArray
from astropy import constants as const
from astropy import units as unit
from optparse import OptionParser

########### Change these here or with command-line options ###################

tcool_tff_ratio = 1
halo_profile    = 1  # 1 = isothermal, 2 = isentropic

####### gas parameters ######
# changing these will just rescale the units
T0              = YTQuantity(1e6,   'K')    
rho0            = YTQuantity(1e-27, 'g/cm**3') 
g0              = YTQuantity(5e-10, 'cm/s**2')
gsoft_scale     = 0.1    # a / H



perturb_amplitude     = 0.02   

# Tmin should be ~ T0/20
Tmin = YTQuantity(5e4, 'K')


###### perturbatio default parameters
default_n = 128
default_kmin = 4
default_kmax = 32
default_f_solenoidal = 2./3.
default_alpha = 0
default_seed = 4085281318



###### computing parameters
nodes = 1
tasks_per_node = 32
num_cooling_cycles = 10
num_outputs = 100


######### Don't change these ############
# (unless you really want to)

def calculate_density(rho0, z, a, H, profile):
    scale = np.sqrt(1 + (z/a)**2) - 1 
    if profile == 1: 
        # isothermal
        return rho0 * np.exp(-(a/H) * scale)

    elif profile == 2:
        # isentropic
        base = (1 - (gamma - 1)*a/(gamma*H) * scale)
        return rho0 * np.power(base, 1.0 / (gamma-1))
        
def calculate_temperature(T0, z, a, H, profile):
    if profile == 1:
        # isothermal
        return T0
    elif profile == 2:
        # isentropic
        scale = np.sqrt(1 +(z/a)**2) - 1
        return T0 * (1 - (gamma - 1)*a/(gamma*H) * scale)
    
def calculate_free_fall_time(z, g0):
    return np.sqrt(2*z / g0)


mu          = 1.22
gamma       = 5./3.

# cosntant definitions
G      = YTQuantity(const.G.cgs.value, 'cm**3/g/s**2')
kb     = YTQuantity(const.k_B.cgs.value, 'cm**2*g/s**2/K')
mh     = YTQuantity(const.m_p.cgs.value, 'g')


# solve for the remaining variables 
# scale height and gravitational softening in kpc
H = kb*T0 / mu / mh / g0
a = gsoft_scale * H


# free fall time and density at the scale height, in seconds
tff_H = calculate_free_fall_time(H, g0)
rho_H = calculate_density(rho0, H, a, H, halo_profile)

tcool_over_L0_H = kb*mu*mh*np.power(T0, 5/3) / rho_H / (gamma - 1)
Lambda0 = tcool_over_L0_H / (tcool_tff_ratio * tff_H)

# assuming range of tcool_tff_ratio from 0.1 - 10
LambdaMin = tcool_over_L0_H / (10  * tff_H)
LambdaMax = tcool_over_L0_H / (0.1 * tff_H)


# estimate cell volume 
cell_volume = (6*H / 128)**3

LengthScale = 1    # in units of scale height H
TimeScale   = 1    # in units of free-fall time at the scale height 
MassScale   = 1e3  # the value, in code-units, of the mass in a cell at the scale height

LengthUnits        = LengthScale * H
MassUnits          = MassScale * rho0 * cell_volume
DensityUnits       = MassUnits / np.power(LengthUnits, 3)
TimeUnits          = TimeScale * tff_H
VelocityUnits      = LengthUnits / TimeUnits
AccelerationUnits  = LengthUnits / TimeUnits / TimeUnits
EnergyUnits        = MassUnits * LengthUnits**2 / TimeUnits**2
GravitationalConstant      = 4 * np.pi*G * DensityUnits*TimeUnits**2


velocity_test = 1e7 / VelocityUnits
ethermal_test =  kb*T0 / mu / mh / (gamma - 1) / VelocityUnits**2
ethermal_min = kb*Tmin / mu / mh / (gamma - 1) / VelocityUnits**2
