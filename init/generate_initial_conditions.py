from astropy import constants as const
from optparse import OptionParser
import matplotlib.pylab as plt
import numpy as np
import shutil
import os

from constants_and_parameters import *

def calculate_density(rho0, z, a, H, profile, inv_beta, eta, T_power_law_index):
    zscale = z/a
    total_pressure_factor = (1.0 + inv_beta + eta)
    halo_scale = a * (np.sqrt(1.0 + zscale*zscale)-1.0) / H / total_pressure_factor

    if profile == 1:
        # isothermal
        rho =  rho0 * np.exp(-halo_scale) 
        return rho

    elif profile == 2:
        # isentropic
        base = (1 - (gamma - 1)/(gamma) * halo_scale)
        return rho0 * np.power(base, 1.0 / (gamma-1))

    elif profile == 3:
        # isocool
        alpha = T_power_law_index
        base = (1 - halo_scale/(2-alpha)) 
        return rho0 * np.power(base, 1-alpha)
        

def calculate_temperature(T0, z, a, H, profile, inv_beta, eta, T_power_law_index):
    zscale = z/a
    total_pressure_factor = (1.0 + inv_beta + eta)
    halo_scale = a * (np.sqrt(1.0 + zscale*zscale)-1.0) / H / total_pressure_factor
    if profile == 1:
        # isothermal
        return T0
    elif profile == 2:
        # isentropic
        scale = np.sqrt(1 +(z/a)**2) - 1
        return T0 * (1 - (gamma - 1)/(gamma)*halo_scale)
    elif profile == 3:
        # isocool
        alpha = T_power_law_index
        return T0 * (1 - halo_scale / (2 - alpha)) 

def calculate_free_fall_time(z, g0):
    return np.sqrt(2*z / g0)


def calculate_cooling_rate(T, Lambda0, T_power_law_index, smooth_factor, Tmin = 5e4, Tmax = 1e8):
    # only the power law
    cooling_rate = Lambda0 * np.power(T, T_power_law_index)
    cooling_rate[T < Tmin] = Lambda0 * np.power(Tmin, T_power_law_index)
    cooling_rate[T > Tmax] = Lambda0 * np.power(Tmax, T_power_law_index)

    logT = np.log10(T)
    logT_min = np.log10(Tmin)
    logT_max = np.log10(Tmax)

    # smooth the edges
    Tminscale =  (logT - logT_min) / (smooth_factor * logT_min) 
    Tmaxscale = -(logT - logT_max) / (smooth_factor * logT_max)
    
    ymin = 1e-20 * Lambda0
    if Lambda0 != 1:
        ymin = ymin.d
    smooth_min = 1
    smooth_max = 1
    smooth_min = (np.tanh(Tminscale - 2.0) + 1.0) / 2. * (1.0 - ymin)  + ymin
    smooth_max = (np.tanh(Tmaxscale - 2.0) + 1.0) / 2. * (1.0 - ymin) + ymin
    return cooling_rate * smooth_min * smooth_max 


def calculate_cooling_time(Lambda0, rho0, T0, z, a, H, profile, T_power_law_index, inv_beta, eta, smooth_factor):
    rho = calculate_density(rho0, z, a, H, profile, inv_beta, eta, T_power_law_index)
    n = rho / (mu * mh)
    T = calculate_temperature(T0, z, a, H, profile, inv_beta, eta, T_power_law_index)

    eth = n * kb * T / (gamma - 1)
    d_eth = np.power(n, 2) * calculate_cooling_rate(T, Lambda0, T_power_law_index, smooth_factor)
    cooling_time = eth / d_eth
    return cooling_time

def plot_cooling_rate(Lambda0, T_power_law_index, smooth_factor, Tmin, Tmax):
    T_list = YTArray(np.logspace(np.log10(Tmin)-1, np.log10(Tmax)+1, 100), 'K')
    
    cool_rate = calculate_cooling_rate(T_list, Lambda0, T_power_law_index, smooth_factor, Tmin = Tmin, Tmax = Tmax)
    plt.loglog(T_list, cool_rate, linewidth = 3)
    plt.loglog(T_list, Lambda0*T_list**(T_power_law_index), color = 'black', linestyle = 'dashed', \
               label = '$ T = \Lambda_0 T^{%0.2f}$'%(T_power_law_index), linewidth = 2)
    
    plt.xlim(Tmin/5, Tmax*5)
    ymax = max(Lambda0*T_list**(T_power_law_index))
    plt.ylim(ymax * 1e-3, ymax)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Cooling Rate (erg cm$^3$ s$^{-1}$)')
    plt.legend()
    plt.savefig("cooling_rates.png", dpi = 300)


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
    outf.write("ProblemType \t\t = 417\n")
    outf.write("TopGridRank \t\t = %i\n\n"%grid_rank)

    if grid_rank == 2:
        res_x = resolution
        res_y = res_x
        box_x = box_z
        box_y = box_z
        outf.write("TopGridDimensions \t = %i %i\n"%(resolution, resolution))
        outf.write("DomainLeftEdge \t\t = %i %i\n"%(-box_x, -box_y))
        outf.write("DomainRightEdge \t = %i %i\n\n"%(box_x,  box_y))
        outf.write("UserDefinedRootGridLayout = 1 64\n")
        outf.write(" # note: can't have fewer grid divisions than number of cores\n\n")

        outf.write("LeftFaceBoundaryCondition   = 3 6 # 3 = Periodic, 6 = hydrostatic\n")
        outf.write("RightFaceBoundaryCondition  = 3 6\n\n")
    else:
        if skinny:
            res_x = int(resolution / 4);
            res_y = res_x
            box_x = box_z / 4.
            box_y = box_z / 4.
        else:
            res_x = resolution
            res_y = res_x
            box_x = box_z
            box_y = box_z

        outf.write("TopGridDimensions \t = %i %i %i\n"%(res_x,  res_y, resolution))
        outf.write("DomainLeftEdge \t\t = %f %f %f\n"%(-box_x, -box_y, -box_z))
        outf.write("DomainRightEdge \t = %f %f %f\n\n"%(box_x,  box_y,  box_z))
        outf.write("UserDefinedRootGridLayout = 1 1 64\n")
        outf.write(" # note: can't have fewer grid divisions than number of cores\n\n")
        outf.write("LeftFaceBoundaryCondition   = 3 3 6 # 3 = Periodic, 6 = hydrostatic\n")
        outf.write("RightFaceBoundaryCondition  = 3 3 6\n\n")
    
    outf.write("# Hierarchy Control Parameters\n")
    outf.write("StaticHierarchy           = 1      # 0 = AMR, 1 = No AMR\n")
    outf.write("MaximumRefinementLevel    = 0\n")

    outf.write("# Hydrodynamics Parameters\n")
    outf.write("UseHydro \t\t = 1\n")
    outf.write("HydroMethod \t\t = 4  # 4 = Dedner MHD\n")
    outf.write("RiemannSolver \t\t = 3  # 3 = LLF\n")
    outf.write("FluxCorrection \t\t = 1\n")
    outf.write("DualEnergyFormalism \t = 1\n")
    outf.write("Gamma \t\t\t = 1.66666666666667\n")
    outf.write("CourantSafetyNumber \t = 0.2\n\n")

    outf.write("# I/O Parameters \n")
    outf.write("StopCycle \t\t = 10000000\n")
    outf.write("InitialTime \t\t = 0.0\n")
    outf.write("StopTime \t\t = %f\t # stop time is %0.2f cooling times (at the scale height) \n"%(StopTime, num_cooling_cycles))
    outf.write("dtDataDump \t\t = %f \n"%(dtDataDump))
    outf.write("DataDumpDir \t\t = DD\n")
    outf.write("DataDumpName \t\t = DD\n")
    outf.write("OutputTemperature \t = 1\n")
    outf.write("OutputCoolingTime \t = 1\n")
    outf.write("WriteExternalAccel \t = 1\n\n")

    outf.write("# Cooling Physics\n")
    outf.write("RadiativeCooling \t = 1\n")
    outf.write("RadiativeCoolingModel  \t = 5\n")
    outf.write("RadiativeHeating \t = 1\n")
    outf.write("RadiativeCoolingFunctionConstant   = %e\n"%Lambda0)
    outf.write("RadiativeCoolingPowerLawIndex      = %0.8f\n"%(T_power_law_index))
    outf.write("RadiativeCoolingSmoothingFactor    = %f\n"%smooth_factor)
    outf.write("RadiativeCoolingMinimumTemperature = %e\n"%(T_min))
    outf.write("RadiativeCoolingMaximumTemperature = %e\n"%(T_max))
    outf.write("MultiSpecies \t\t = 0\n")
    outf.write("MetalCooling \t\t = 0\n")
    outf.write("CIECooling \t\t = 0\n")
    outf.write("UseCoolingTimestep \t = 1\n")
    outf.write("CoolingTimestepSafetyFactor \t = 0.4\n")
    outf.write("CoolingTimestepMinimumTimestep = 20000 # years\n\n")

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
    outf.write("TIPerturbationAmplitude\t = %f\n"%(perturbation_amplitude))
    outf.write("TIPerturbationType\t\t = %i\n"%perturb_type)
    outf.write("TIHaloProfile \t\t\t = %i\n"%(halo_profile))
    outf.write("TIMagneticFieldDirection = %f %f %f\n"%(bfield_direction[0], bfield_direction[1], bfield_direction[2]))
    if constant_B_pressure == 1:
            outf.write("TIMagneticFieldUseConstantBeta    = 0\n")

    outf.write("TIMagneticFieldInverseBeta      = %e\n"%(magnetic_pressure_ratio))
    outf.write("TestProblemUseMetallicityField\t = 0 \n")
    if cr_pressure_ratio > 0 :
        if constant_cr_pressure == 1:
            outf.write("TICosmicRayUseConstantEta    = 0\n")
        outf.write("TICosmicRayPressureRatio   = %e\n"%(cr_pressure_ratio))
        outf.write("CRModel                    = 1\n")
        outf.write("CRCourantSafetyNumber      = 0.2\n")
        outf.write("CRDiffusion                = %i\n"%cr_diffusion)
        outf.write("CRkappa                    = %e\n"%cr_kappa)
        outf.write("CRStreaming                = %i\n"%cr_streaming)
        outf.write("CRStreamStabilityFactor = %e\n"%cr_streaming_stability)
        outf.write("CRHeating\t                = %i\n"%cr_heating)
        outf.write("CRdensFloor                = 1e-20\n")

            
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
    outf.write("Cooling time at H = %e years\n"%calculate_cooling_time(Lambda0, rho0, T0, H, a, H, halo_profile, \
                T_power_law_index, magnetic_pressure_ratio, cr_pressure_ratio, smooth_factor).in_units('yr'))
    outf.write("Cooling_time at 3H = %e years\n"%calculate_cooling_time(Lambda0, rho0, T0, 3*H, a, H, halo_profile, \
                    T_power_law_index, magnetic_pressure_ratio, cr_pressure_ratio, smooth_factor).in_units('yr'))
    outf.close()


def generate_sbatch_file(nodes = 2, tasks_per_node = 32):
    cwd = os.getcwd()
    basename = os.path.basename(cwd)
    outf = open('submit.sbatch', 'w')
    outf.write('#!/bin/bash\n')
    outf.write('#SBATCH --job-name=%s\n'%basename)
    outf.write('#SBATCH --nodes=%i\n'%nodes)
    outf.write('#SBATCH --ntasks-per-node=%i\n'%(tasks_per_node))
    outf.write('#SBATCH --time=%s\n\n'%wall_time)
    
    outf.write('module purge\n')
    outf.write('module load slurm gcc lib/hdf5 openmpi\n\n')

    outf.write('cd %s\n'%cwd)

    outf.write('c=0\nwhile true\ndo\n')
    outf.write('    if [ ! -e estd.out.$c ]; then\n        break\n    fi\n')
    outf.write('    c=$[$c+1]\ndone\nmv estd.out estd.out.$c\n\n')

    procs = nodes * tasks_per_node
    outf.write('mpirun -n %i ~/enzo-dev-cr/src/enzo/enzo.exe -d ThermalInstability.enzo >&estd.out\n'%(procs))
    outf.write('#mpirun -n %i ~/enzo-dev-cr/src/enzo/enzo.exe -d -r DD0100/DD0100 >&estd.out\n'%(procs))
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

mu              = 1.22
gamma           = 5./3.

# cosntant definitions                  
G      = YTQuantity(const.G.cgs.value, 'cm**3/g/s**2')
kb     = YTQuantity(const.k_B.cgs.value, 'cm**2*g/s**2/K')
mh     = YTQuantity(const.m_p.cgs.value, 'g')


# scale height and gravitational softening in kpc
H = kb*T0 / mu / mh / g0
a = gsoft_scale * H


# free fall time and density at the scale height, in seconds 
tff_H = calculate_free_fall_time(H, g0)
rho_H = calculate_density(rho0, H, a, H, halo_profile, \
                          magnetic_pressure_ratio, cr_pressure_ratio, T_power_law_index)
T_H = calculate_temperature(T0, H, a, H, halo_profile, magnetic_pressure_ratio, \
                            cr_pressure_ratio, T_power_law_index)

tcool_over_L0_H = calculate_cooling_time(1, rho0, T0, H, a, H, halo_profile,\
            T_power_law_index, magnetic_pressure_ratio, cr_pressure_ratio, smooth_factor)
Lambda0 = tcool_over_L0_H / (tcool_tff_ratio * tff_H)

# calculate cr diffusion coefficient
tcr = tcr_tff_ratio * tff_H
cr_kappa = H * H / tcr


if halo_profile < 2:
    box_x, box_y, box_z   = [3, 3, 3]  # box dimensions are (2*box_x*box_y*box_z)**3  
else:
    box_x, box_y, box_z   = [2, 2, 2]

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
# pert.generate_perturbation_infile()

plot_cooling_rate(Lambda0, T_power_law_index, smooth_factor, T_min, T_max)

# generate enzo input file
generate_enzo_input_file()
check_units()

# generate submit.sbatch file
generate_sbatch_file(nodes = nodes, tasks_per_node = tasks_per_node)

generate_clean_file()
