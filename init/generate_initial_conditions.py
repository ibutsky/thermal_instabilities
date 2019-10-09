import matplotlib.pylab as plt
import shutil
import os

import perturbation as pert
from constants_and_parameters import *

def calculate_cooling_rate(T, Lambda0, logT_min = 5.0):
    cooling_rate = Lambda0 * np.power(T, -2./3.)
#    cooling_rate[T < np.power(10, logT_min)] = 1e-40
    return cooling_rate

def generate_cooling_input_file(Lambda0, logT_min = 5.0, logT_max = 7.0):
    temperature_list = np.logspace(logT_min, logT_max, 500)
    cooling_rate     = calculate_cooling_rate(temperature_list, Lambda0, logT_min = logT_min)    

    outf = open('cool_rates.in', 'w')
    outf.write("# Cooling table generated for Lambda0 = %e\n"%Lambda0)
    outf.write("# Cooling rate: L(T) = Lambda0 * T**-(2/3) (erg cm^3/s)\n")
    outf.write("#\t log(T)\t\t log(L(T))\talso log(L(T))\n")
    
    for  T, L in zip(np.log10(temperature_list), np.log10(cooling_rate)):
        outf.write("\t%f\t%f\t%f\n"%(T, L, L))
    outf.close()

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
    shutil.copyfile('base_enzo_parameters.enzo', 'ThermalInstability.enzo')
    outf = open('ThermalInstability.enzo', 'a')
    tff_scale = TimeScale * tff_H# / yr_s
    StopTime = num_cooling_cycles * tcool_tff_ratio
    dtDataDump = StopTime / num_outputs
    
    outf.write("# Units Parameters\n")
    outf.write("DensityUnits\t\t = %10e\n"%DensityUnits)
    outf.write("LengthUnits\t\t = %10e\t # 1 code length unit = %f kpc\n"%(LengthUnits, (LengthScale*H).in_units('kpc')))
    outf.write("# 1 code time unit = 1 free-fall time at the scale height (%e years) \n"%(tff_scale.in_units('yr')))
    outf.write("TimeUnits\t\t = %10e \n"%(TimeUnits))
    outf.write("StopTime\t\t = %f\t # stop time is %0.2f cooling times (at the scale height) \n"%(StopTime, num_cooling_cycles))
    outf.write("dtDataDump\t\t = %f \n"%(dtDataDump))
    outf.write("GravitationalConstant\t = %10e\t # 4*pi*G*DensityUnits*TimeUnits^2\n\n"%GravitationalConstant)
    
    outf.write("# Gravity Parameters\n")
    outf.write("ExternalGravity\t\t\t = 2\n")
    outf.write("ExternalGravityConstant\t\t = %6e\t # in cm/s**2\n"%g0)
    outf.write("ExternalGravityPosition\t\t = 0 0 0\t # in code units\n")
    outf.write("ExternalGravitySofteningRadius\t = %e\t # in kpc\n\n"%(a.in_units('kpc')))

    outf.write("# Problem-specific Parameters\n")
    outf.write("TIMeanDensity\t\t\t = %e  # g/cc\n"%(rho0))
    outf.write("TIMeanTemperature\t\t = %e  # Kelvin\n"%(T0))
    outf.write("TIDensityPerturbationAmplitude\t = %f\n"%(perturb_amplitude))
    outf.write("TIHaloProfile\t\t =%i\n"%(halo_profile))
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


######### GENERATE INITIAL CONDITIONS AND FILES  ########

# generate perturbation input file: 
pert.generate_perturbation_infile()

# generate cooling input file
generate_cooling_input_file(Lambda0, np.log10(Tmin))
plot_cooling_rate_range()

# generate enzo input file
generate_enzo_input_file()
check_units()

# generate submit.sbatch file
generate_sbatch_file(nodes = nodes, tasks_per_node = 32)
