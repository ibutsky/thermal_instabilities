import numpy as np
from astropy import constants as const
from astropy import units as unit
from optparse import OptionParser

########### Change these here or with command-line options ###################

###### gas default parameters
tcool_tff_ratio = 1
T0              = 1e6    # K 
rho0            = 1e-27  # cgs
g0              = 5e-10  # cm/s**2
gsoft_scale     = 0.1    # a / H

perturb_amplitude     = 0.02   

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


######### Don't change these ############
# (unless you really want to)

mu          = 1.22
gamma       = 5./3.

# cosntant definitions
G      = const.G.cgs.value
kb     = const.k_B.cgs.value
mh     = const.m_p.cgs.value
kpc_cm = const.kpc.cgs.value
yr_s   = unit.yr.to('s')
msun_g = const.M_sun.cgs.value



# solve for the remaining variables 
# scale height and gravitational softening in kpc
H = kb*T0 / mu / mh / g0
a = gsoft_scale * H



# free fall time and density at the scale height, in seconds
tff_H = np.sqrt(2*H / g0)
rho_H = rho0 * np.exp(- (a/H) * (np.sqrt(1 + (H/a)**2) - 1))
tcool_over_L0_H = 3*kb*mu*mh*np.power(T0, 5/3) / rho_H
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
