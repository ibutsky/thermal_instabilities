from yt import YTQuantity
from yt import YTArray

########### These parameters will likely change between runs ########

tcool_tff_ratio = 1    # ratio of t_cool / t_free_fall 
halo_profile    = 3    # 1 = isothermal, 2 = isentropic, 3 = iso-tcool
perturb_type    = 1

# note: beta = p_mag / p_gas 
constant_B_pressure     = 0
magnetic_pressure_ratio = 0
constant_cr_pressure    = 0
cr_pressure_ratio       = 0

cr_streaming = 0
cr_heating   = 0
cr_diffusion = 0

########### These parameters don't really need to change ###########
bfield_direction    = [0, 1, 0]
resolution          = 128  # resolution along z-axis; scaled for x-y axes
grid_rank           = 3

# cooling function parameters
# Tmin should be ~ T0/20    
T_min = YTQuantity(5e4, 'K')
T_max = YTQuantity(1e9, 'K')
T_power_law_index = (-2./3.)
smooth_factor = 0.02

####### gas parameters ######
# changing these will just rescale the units
T0              = YTQuantity(1e6,   'K')    
rho0            = YTQuantity(1e-27, 'g/cm**3') 
g0              = YTQuantity(5e-10, 'cm/s**2')
gsoft_scale     = 0.1    # units of scale height

###### perturbatio default parameters
perturbation_amplitude = 0.02
default_n              = resolution
default_kmin           = 4
default_kmax           = int(resolution / 2)
default_f_solenoidal   = 2./3.
default_alpha          = 0
default_seed           = 4088338814


###### computing parameters
nodes              = 1
tasks_per_node     = 32
num_cooling_cycles = 10  # simulation stop time = 20 * t_cool
num_outputs        = 100
wall_time          = '5:00:00'


