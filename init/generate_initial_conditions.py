from astropy import constants as const
from optparse import OptionParser
import matplotlib.pylab as plt
import numpy as np
import shutil
import os

import perturbation as pert
from constants_and_parameters import *

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


def calculate_cooling_rate(T, Lambda0, logT_min = 5.0):
    cooling_rate = Lambda0 * np.power(T, -2./3.)
    return cooling_rate


def calculate_cooling_time(Lambda0, rho0, T0, z, a, H, profile):
    rho = calculate_density(rho0, z, a, H, profile)
    n = rho / (mu * mh)
    T = calculate_temperature(T0, z, a, H, profile)

    eth = n * kb * T / (gamma - 1)
    d_eth = np.power(n, 2) * calculate_cooling_rate(T, Lambda0)
    cooling_time = eth / d_eth
    return cooling_time
    

def plot_cooling_rate_range():
    T_list = np.logspace(5, 7, 100)
    plt.loglog(T_list, calculate_cooling_rate(T_list, LambdaMin), label = "tcool / tff = 10")
    plt.loglog(T_list, calculate_cooling_rate(T_list, LambdaMax), label = "tcool / tff = 0.1")
    plt.xlim(8e4, 3e6)
    plt.ylim(1e-24, 3e-21)
    plt.xlabel("Temperature (K)")
    plt.ylabel("$\Lambda$(T) (erg cm^3/s)")
    plt.legend()
    plt.savefig("cooling_rates.png")


def generate_enzo_input_file():
    # Generate final enzo parameter file
    outf = open('ThermalInstability.enzo', 'w')
    tff_scale = TimeScale * tff_H
    StopTime = num_cooling_cycles * tcool_tff_ratio
    dtDataDump = StopTime / num_outputs
    outf.write("# Thermal Instability in Turbulent Box\n#\n")
    outf.write("# Initializes a 3D grid of turbulent density perturbations\n")
    outf.write("# in pressure equilibrium\n\n")

    outf.write("# Initialization Parameters\n")
    outf.write("Problem Type \t\t = 417\n")
    outf.write("TopGridRank \t\t = 3\n\n")
    
    res_x = int(resolution * box_x / box_z)
    res_y = int(resolution * box_y / box_z)
    outf.write("TopGridDimensions \t = %i %i %i\n"%(res_x,  res_y, resolution))
    outf.write("DomainLeftEdge \t\t = %i %i %i\n"%(-box_x, -box_y, -box_z))
    outf.write("DomainRightEdge \t = %i %i %i\n\n"%(box_x,  box_y,  box_z))
    
    outf.write("LeftFaceBoundaryCondition   = 3 3 6 # 3 = Periodic, 6 = hydrostatic\n")
    outf.write("RightFaceBoundaryCondition  = 3 3 6\n\n")
    
    outf.write("# Hierarchy Control Parameters\n")
    outf.write("StaticHierarchy           = 1      # 0 = AMR, 1 = No AMR\n")
    outf.write("MaximumRefinementLevel    = 0\n")
    outf.write("UserDefinedRootGridLayout = 1 1 %i"%(int(resolution/2)))
    outf.write(" # note: can't have fewer grid divisions than number of cores\n\n")

    outf.write("# Hydrodynamics Parameters\n")
    outf.write("UseHydro \t\t = 1\n")
    outf.write("HydroMethod \t\t = 4  # 4 = Dedner MHD\n")
    outf.write("RiemannSolver \t\t = 3  # 3 = LLF\n")
    outf.write("FluxCorrection \t\t = 1\n")
    outf.write("DualEnergyFormalism \t = 1\n")
    outf.write("Gamma \t\t\t = 1.66666666666667\n")
    outf.write("CourantSafetyNumber \t = 0.4\n\n")

    outf.write("# I/O Parameters \n")
    outf.write("StopCycle \t\t = 1000000\n")
    outf.write("InitialTime \t\t = 0.0\n")
    outf.write("StopTime \t\t = %f\t # stop time is %0.2f cooling times (at the scale height) \n"%(StopTime, num_cooling_cycles))
    outf.write("dtDataDump \t\t = %f \n"%(dtDataDump))
    outf.write("DataDumpDir \t\t = DD\n")
    outf.write("OutputTemperature \t = 1\n")
    outf.write("OutputCoolingTime \t = 1\n")
    outf.write("WriteExternalAccel \t = 1\n\n")

    outf.write("# Cooling Physics\n")
    outf.write("RadiativeCooling \t = 1\n")
    outf.write("RadiativeCoolingModel  \t = 5\n")
    outf.write("RadiativeHeating \t = 1\n")
    outf.write("RadiativeCoolingFunctionConstant   = %e\n"%Lambda0)
    outf.write("RadiativeCoolingPowerLawIndex      = %0.8f\n"%(T_power_law_index))
    outf.write("RadiativeCoolingMinimumTemperature = %e\n"%(T_min))
    outf.write("RadiativeCoolingMaximumTemperature = %e\n"%(T_max))
    outf.write("MultiSpecies \t\t = 0\n")
    outf.write("MetalCooling \t\t = 0\n")
    outf.write("CIECooling \t\t = 0\n")
    outf.write("UseCoolingTimestep \t = 1\n")
    outf.write("CoolingTimestepSafetyFactor \t = 0.01\n")
    outf.write("CoolingTimestepMinimumTimestep = 0 # years\n\n")

    outf.write("# Global Parameters\n")
    outf.write("tiny_number \t\t = 1.0e-20\n")
    outf.write("MinimumEfficiency \t = 0.4\n")
    outf.write("Mu \t\t\t = 1.22\n\n")
    
    outf.write("# Units Parameters\n")
    outf.write("DensityUnits\t\t = %10e\n"%DensityUnits)
    outf.write("LengthUnits\t\t = %10e\t # 1 code length unit = %f kpc\n"%(LengthUnits, (LengthScale*H).in_units('kpc')))
    outf.write("TimeUnits\t\t = %10e  "%(TimeUnits))
    outf.write("# 1 code time unit = 1 t_ff (s) at the scale height \n"%(tff_scale.in_units('yr')))
    outf.write("GravitationalConstant\t = %10e\t # 4*pi*G*DensityUnits*TimeUnits^2\n\n"%GravitationalConstant)
    
    outf.write("# Gravity Parameters\n")
    outf.write("ExternalGravity\t\t\t = 2\n")
    outf.write("ExternalGravityConstant\t\t = %6e\t # in cm/s**2\n"%g0)
    outf.write("ExternalGravityPosition\t\t = 0 0 0\t # in code units\n")
    outf.write("ExternalGravitySofteningRadius\t = %e\t # in kpc\n\n"%(a.in_units('kpc')))

    outf.write("# Problem-specific Parameters\n")
    outf.write("TIMeanDensity\t\t\t = %e  # g/cc\n"%(rho0))
    outf.write("TIMeanTemperature\t\t = %e  # Kelvin\n"%(T0))
    outf.write("TIDensityPerturbationAmplitude\t = %f\n"%(perturbation_amplitude))
    outf.write("TIHaloProfile \t\t\t = %i\n"%(halo_profile))
    outf.write("TestProblemUseMetallicityField\t = 0 \n")
    outf.close()


def check_units():
    outf = open('unit_sanity_check.dat', 'w')
    outf.write("\n\n************ Unit Sanity Check ***********************\n")
    outf.write("rho0\t  = %e in code units\n"%(rho0/DensityUnits))
    outf.write("g0\t  = %e in code units\n"%(g0 / AccelerationUnits))
    outf.write("e_thermal = %e in code units\n"%ethermal_test)
    outf.write("1 cm/s\t  = %e in code units\n"%(velocity_test/1e7))
    outf.write("100 km/s  = %e in code units\n\n"%velocity_test)

    outf.write("Free fall time at H = %e years\n"%(tff_H.in_units('yr')))
    outf.write("Cooling time at H = %e years\n"%calculate_cooling_time(Lambda0, rho0, T0, H, a, H, halo_profile).in_units('yr'))
    outf.write("Cooling_time at 3H = %e years\n"%calculate_cooling_time(Lambda0, rho0, T0, 3*H, a, H, halo_profile).in_units('yr'))
    outf.close()


def generate_sbatch_file(nodes = 2, tasks_per_node = 32):
    cwd = os.getcwd()
    basename = os.path.basename(cwd)
    outf = open('submit.sbatch', 'w')
    outf.write('#!/bin/bash\n')
    outf.write('#SBATCH --job-name=%s\n'%basename)
    outf.write('#SBATCH --nodes=%i\n'%nodes)
    outf.write('#SBATCH --ntasks-per-node=%i\n'%(tasks_per_node))
    outf.write('#SBATCH --time=1:00:00\n\n')
    
    outf.write('module purge\n')
    outf.write('module load slurm gcc lib/hdf5 openmpi\n\n')

    outf.write('cd %s\n'%cwd)

    outf.write('c=0\nwhile true\ndo\n')
    outf.write('    if [ ! -e estd.out.$c ]; then\n        break\n    fi\n')
    outf.write('    c=$[$c+1]\ndone\nmv estd.out estd.out.$c\n\n')

    procs = nodes * tasks_per_node
    outf.write('mpirun -n %i ~/enzo-dev-cr/src/enzo/enzo.exe -d ThermalInstability.enzo >&estd.out\n'%(procs))
    outf.close()

def generate_clean_file():
    outf = open('clean.sh', 'w')
    outf.write("rm *.out*\n")
    outf.write("rm *~\n")
    outf.write("rm Output*\n")
    outf.write("rm Enzo*\n")
    outf.write("rm Evtime\n")
    outf.write("rm RunFinished\n")
    outf.write("rm -rf DD*\n")
    outf.close()

######### GENERATE INITIAL CONDITIONS AND FILES  ########
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

tcool_over_L0_H = kb*mu*mh*np.power(T0, 1 - T_power_law_index) / rho_H / (gamma - 1)
Lambda0 = tcool_over_L0_H / (tcool_tff_ratio * tff_H)

# assuming range of tcool_tff_ratio from 0.1 - 10                                                                                
LambdaMin = tcool_over_L0_H / (10  * tff_H)
LambdaMax = tcool_over_L0_H / (0.1 * tff_H)

LengthScale  = 1    # in units of scale height H
TimeScale    = 1    # in units of free-fall time at the scale height
DensityScale = 10   # in units of rho0                                                                                           

LengthUnits        = LengthScale * H
DensityUnits       = rho0 / DensityScale
MassUnits          = DensityUnits * np.power(LengthUnits, 3)
TimeUnits          = TimeScale * tff_H
VelocityUnits      = LengthUnits / TimeUnits
AccelerationUnits  = LengthUnits / TimeUnits / TimeUnits
EnergyUnits        = MassUnits * LengthUnits**2 / TimeUnits**2
GravitationalConstant      = 4 * np.pi*G * DensityUnits*TimeUnits**2


velocity_test = 1e7 / VelocityUnits
ethermal_test =  kb*T0 / mu / mh / (gamma - 1) / VelocityUnits**2
ethermal_min = kb*T_min / mu / mh / (gamma - 1) / VelocityUnits**2

# generate perturbation input file: 
pert.generate_perturbation_infile()

plot_cooling_rate_range()

# generate enzo input file
generate_enzo_input_file()
check_units()

# generate submit.sbatch file
generate_sbatch_file(nodes = nodes, tasks_per_node = 32)

generate_clean_file()