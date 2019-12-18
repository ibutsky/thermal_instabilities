import yt
from yt import YTArray
from yt import YTQuantity
import sys
import os
import numpy as np
import matplotlib.pylab as plt
import palettable

import plotting_tools as pt

def calculate_mass_flux(sim_folder, tctf, output_list, grid_rank = 3, Tcold = 3.33333e5):
    time_list     = np.array([])
    mass_flux_list = np.array([])
    if not os.path.isdir(sim_folder):
        return time_list, mass_flux_list
    sim_base = os.path.basename(sim_folder)
    print(sim_base)
    
    out_name = '../../data/mass_flux_%s'%(sim_base)
    if os.path.isfile(out_name) and load == True:
        time_list, mass_flux_list = np.loadtxt(out_name, unpack=True)
    
    else:
        for output in output_list:
            ds_path = "%s/DD%04d/DD%04d"%(sim_folder, output, output)

            if os.path.isfile(ds_path):
                ds = yt.load(ds_path)
                print(zstart, zend)

                if (grid_rank == 3):
                    ad = ds.all_data()
                    all_z = ad[('gas', 'z')]
                    zmask1 = (all_z / ds.length_unit.in_units('kpc') > zstart) & (all_z / ds.length_unit.in_units('kpc') < zend)
                    zmask2 = (all_z / ds.length_unit.in_units('kpc') < -zstart) & (all_z / ds.length_unit.in_units('kpc') > -zend)

                    v_cool_flow = ds.length_unit / (tctf * ds.time_unit) # H / tff
                    print(v_cool_flow.in_units('km/s'))
                    cool_flow_flux = (ad[('gas', 'density')] * v_cool_flow).in_units('Msun / kpc**2 / yr')

                    all_flux = (ad[('gas', 'density')] * ad[('gas', 'velocity_z')]).in_units('Msun / kpc**2 / yr')
                    
                    flux_mask1 = (zmask1) &  (ad[('gas', 'velocity_z')] < 0) & (ad[('gas', 'temperature')] <= Tcold)
                    flux_mask2 = (zmask2) &  (ad[('gas', 'velocity_z')] > 0) & (ad[('gas', 'temperature')] <= Tcold)
                    
#                    mass_flux = -np.sum(all_flux[flux_mask1] / cool_flow_flux[flux_mask1])
#                    mass_flux += np.sum(all_flux[flux_mask2] / cool_flow_flux[flux_mask2])

                    mass_flux = -all_flux[flux_mask1] / cool_flow_flux[flux_mask1]
                    mass_flux = np.append(mass_flux, (all_flux[flux_mask2] / cool_flow_flux[flux_mask2]))
                else:

                    ad = ds.all_data()
                    all_y = ad[('gas', 'y')]
                    ymask1 = (all_y / ds.length_unit.in_units('kpc') > zstart) & (all_y / ds.length_unit.in_units('kpc') < zend)
                    ymask2 = (all_y / ds.length_unit.in_units('kpc') < -zstart) & (all_y / ds.length_unit.in_units('kpc') > -zend)

                    mass_flux = 0 

                    dx = ad[('gas', 'dx')]
                    all_flux = (ad[('gas', 'density')] * ad[('gas', 'velocity_y')] * dx * dx).in_units('Msun / yr')

                    cold_flux_mask1 = (ymask1) & (ad[('gas', 'temperature')] <= Tcold) & (ad[('gas', 'velocity_y')] < 0)
                    cold_flux_mask2 = (ymask2) & (ad[('gas', 'temperature')] <= Tcold) & (ad[('gas', 'velocity_y')] > 0)
                    

                    mass_flux = -np.sum(all_flux[cold_flux_mask1])
                    mass_flux += np.sum(all_flux[cold_flux_mask2])


                time_list      = np.append(time_list,    ds.current_time / tctf)
                mass_flux_rms = np.sqrt(np.nanmean(mass_flux**2))
                mass_flux_list = np.append(mass_flux_list, mass_flux_rms)

        if save:
            outf = open(out_name, 'w')
            for i in range(len(time_list)):
                outf.write("%e %e\n"%(time_list[i], mass_flux_list[i]))
            outf.close()

    print(mass_flux_list)
    return time_list, mass_flux_list




def plot_mass_flux(sim, beta = 'inf', tctf_list = None, cr_list = None, diff_list = None, \
                                 beta_list = None, work_dir = '../../simulations/', grid_rank = 3):

    fig, ax = plt.subplots(figsize = (6, 6))
    ax.set_yscale('log')
#    ax.set_ylim(1e-5, 3e-3)
    ax.set_xlim(0, 10)
    
    cpal = palettable.cmocean.sequential.Tempo_7_r.mpl_colors 
    output_list = np.linspace(0, 110, 10)

    for i, tctf in enumerate(tctf_list):
        sim_location = pt.get_sim_location(sim, tctf, beta_list[i], cr_list[i], diff = diff_list[i], work_dir = work_dir)
        if not os.path.isdir(sim_location):
            continue

        label = pt.get_label_name(compare, tctf, beta_list[i], cr_list[i])
        time_list, mass_flux_list = calculate_mass_flux(sim_location, tctf, output_list, grid_rank = grid_rank)
        ax.plot(time_list, mass_flux_list, linewidth = 3, label = label, color = cpal[i])
             

    ax.set_xlabel('t / t$_{cool}$')    
    ax.set_ylabel('Cold Mass Flux (M$_{\odot}$ / yr)')
    ax.legend()
    fig.tight_layout()

    plot_name = pt.get_fig_name('cold_flux', sim, compare, tctf_list[0], beta_list[0], \
                              cr_list[0], diff_list[0])

    plt.savefig(plot_name, dpi = 300)

        

work_dir = '../../simulations'
grid_rank = 3

zstart = 0.9
zend = 1.1
save = False
load = False

crdiff = 0

sim = sys.argv[1]
compare = sys.argv[2]

if sys.argv[3]:
    tctf = float(sys.argv[3])


print(compare)
tctf_list, beta_list, cr_list, diff_list = pt.generate_lists(compare, tctf, crdiff = crdiff)
beta_list = len(tctf_list)*[10.0]
#tctf_list = [1]
#beta_list = ['inf']
#cr_list = [0]
#diff_list = [0]


print(tctf_list, beta_list, cr_list, diff_list)

plot_mass_flux(sim, tctf_list = tctf_list, beta_list = beta_list, grid_rank = grid_rank,\
                                cr_list = cr_list, diff_list = diff_list, work_dir = work_dir)
