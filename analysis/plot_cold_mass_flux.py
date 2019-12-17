import yt
from yt import YTArray
from yt import YTQuantity
import sys
import os
import numpy as np
import matplotlib.pylab as plt
import palettable

import plotting_tools as pt

def calculate_mass_flux(sim_folder, output_list, field = 'density', grid_rank = 3, Tcold = 3.33333e5):
    time_list     = np.array([])
    mass_flux_list = np.array([])
    if not os.path.isdir(sim_folder):
        return time_list, dzfield_rms_list
    sim_base = os.path.basename(sim_folder)
    print(sim_base)
    
    out_name = '../../data/mass_flux_%s_%s'%(sim_base, field)
    if os.path.isfile(out_name) and load == True:
        time_list, mass_flux_list = np.loadtxt(out_name, unpack=True)
    
    else:
        for output in output_list:
            ds_path = "%s/DD%04d/DD%04d"%(sim_folder, output, output)

            if os.path.isfile(ds_path):
                ds = yt.load(ds_path)
                print(zstart, zend)

                if (grid_rank == 3):
                    region1 = ds.r[0, 0, zstart:zend]
                    region2 = ds.r[0, 0, -zend:-zstart]
                    zlist    = np.append(region1[('gas', 'z')].in_units('kpc'),\
                                         region2[('gas', 'z')].in_units('kpc'))
                
                    mass_flux = np.array([])
                    for z in zlist:
                        zslice  = ds.r[:, :, YTQuantity(z, 'kpc')]
                        zfield     = zslice[('gas', field)]
                        zfield_ave = np.mean(zfield)
                        dzfield    = np.append(dzfield, (zfield - zfield_ave) / zfield_ave)
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


                time_list      = np.append(time_list,    ds.current_time)
                mass_flux_list = np.append(mass_flux_list, mass_flux)

        if save:
            outf = open(out_name, 'w')
            for i in range(len(time_list)):
                outf.write("%e %e\n"%(time_list[i], mass_flux_list[i]))
            outf.close()

    print(mass_flux_list)
    return time_list, mass_flux_list




def plot_density_fluctuation_growth(sim, beta = 'inf', tctf_list = None, cr_list = None, diff_list = None, \
                                    field = 'density', beta_list = None, work_dir = '../../simulations/', grid_rank = 3):

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
        time_list, mass_flux_list = calculate_mass_flux(sim_location, output_list, field = field, grid_rank = grid_rank)
        ax.plot(time_list/tctf, mass_flux_list, linewidth = 3, label = label, color = cpal[i])
             

    ax.set_xlabel('t / t$_{cool}$')    
    ax.set_ylabel('Cold Mass Flux (M$_{\odot}$ / yr)')
    ax.legend()
    fig.tight_layout()

    plot_name = pt.get_fig_name('cold_flux', sim, compare, tctf_list[0], beta_list[0], \
                              cr_list[0], diff_list[0])

    plt.savefig(plot_name, dpi = 300)

        

work_dir = '../../simulations/2d_256'
grid_rank = 2

field = 'density'
zstart = 0.5
zend = 1.5
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

plot_density_fluctuation_growth(sim, tctf_list = tctf_list, beta_list = beta_list, grid_rank = grid_rank,\
                                cr_list = cr_list, diff_list = diff_list, field = field, work_dir = work_dir)
