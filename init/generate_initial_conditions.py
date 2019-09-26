import matplotlib.pylab as plt
import shutil

import perturbation as pert
from constants_and_parameters import *

def calculate_cooling_rate(T, Lambda0, logT_min = 5.0):
    cooling_rate = Lambda0 * np.power(T, -2./3.)
    cooling_rate[T < np.power(10, logT_min)] = 1e-40
    return cooling_rate

def generate_cooling_input_file(Lambda0, logT_min = 5.0):
    temperature_list = np.logspace(4.99, 7.0, 500)
    cooling_rate     = calculate_cooling_rate(temperature_list, Lambda0, logT_min = logT_min)    

    outf = open('cool_rates.in', 'w')
    outf.write("# Cooling table generated for Lambda0 = %e\n"%Lambda0)
    outf.write("# Cooling rate: L(T) = Lambda0 * T**-(2/3) (erg cm^3/s)\n")
    outf.write("#\t log(T)\t\t log(L(T))\talso log(L(T))\n")
    
    for  T, L in zip(np.log10(temperature_list), np.log10(cooling_rate)):
        outf.write("\t%f\t%f\t%f\n"%(T, L, L))
    outf.close()


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

    outf.write("# Units Parameters\n")
    outf.write("DensityUnits\t\t = %10e\n"%DensityUnits)
    outf.write("LengthUnits\t\t = %10e\t # 1 code length unit = %f kpc\n"%(LengthUnits, LengthScale*H / kpc_cm))
    outf.write("TimeUnits\t\t = %10e\t # 1 code time unit = %e years\n"%(TimeUnits, TimeScale * tff_H / yr_s))
    outf.write("GravitationalConstant\t = %10e\t # 4*pi*G*DensityUnits*TimeUnits^2\n\n"%GravitationalConstant)
    
    outf.write("# Gravity Parameters\n")
    outf.write("SelfGravity\t\t\t = 1\n")
    outf.write("ExternalGravity\t\t\t = 2\n")
    outf.write("ExternalGravityConstant\t\t = %6e\t # in cm/s**2\n"%g0)
    outf.write("ExternalGravityPosition\t\t = 0 0 0\t # in code units\n")
    outf.write("ExternalGravitySofteningRadius\t = %e\t # in kpc\n\n"%(a/kpc_cm))

    outf.write("# Problem-specific Parameters\n")
    outf.write("TIMeanDensity\t\t\t = %e  # g/cc\n"%(rho0))
    outf.write("TIMeanTemperature\t\t = %e  # Kelvin\n"%(T0))
    outf.write("TIDensityPerturbationAmplitude\t = 0.02\n")
    outf.write("TestProblemUseMetallicityField\t = 0 \n")
    outf.close()


def check_units():
    outf = open('unit_sanity_check.dat', 'w')
    outf.write("\n\n************ Unit Sanity Check ***********************\n")
    outf.write("rho0\t  = %e in code units\n"%(rho0/DensityUnits))
    outf.write("g0\t  = %e in code units\n"%(g0 / AccelerationUnits))
    outf.write("e_thermal = %e in code units\n"%ethermal_test)
    outf.write("1 cm/s\t  = %e in code units\n"%(velocity_test/1e7))
    outf.write("100 km/s  = %e in code units\n"%velocity_test)
    outf.close()


######### GENERATE INITIAL CONDITIONS AND FILES  ########

# generate perturbation input file: 
pert.generate_perturbation_infile()

# generate cooling input file
generate_cooling_input_file(Lambda0)
plot_cooling_rate_range()

# generate enzo input file
generate_enzo_input_file()
check_units()
