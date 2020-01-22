import yt
from yt import YTArray
from yt import YTQuantity
import sys
import os
import numpy as np
import matplotlib.pylab as plt
import palettable
import seaborn as sns

import plotting_tools as pt


def plot_cold_density(output, sim, compare, tctf, beta, cr, diff = 0, stream = 0, heat = 0,
                              T_cold = 3.3333333e5, zstart = 0.8, zend = 1.2,
                              work_dir = '../../simulations/', grid_rank = 3):

    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list\
                        = pt.generate_lists(compare, tctf, crdiff = diff, cr = cr)
    print(tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list)
    diff_list = [0, 1, 3, 0, 0, 0]
    stream_list = [0, 0, 0, 1, 1, 1]
    heat_list = [0, 0, 0, 0, 1, 1]
    beta_list = [100, 100, 100, 100, 100, 10]
    linestyle_list = ['solid', 'dashed', 'dashed', 'dotted', 'dotted', 'dotted']
    
    fig, ax = plt.subplots(figsize = (4.4, 4))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(.09, 5)
    ax.set_ylim(1, 1e2)
    
    cpal = palettable.cmocean.sequential.Deep_6_r.mpl_colors
    for j, diff in enumerate(diff_list):
        for i, cr in enumerate([0, 0.01, 0.1, 1, 10]):
            tctf_list = []
            cold_dens_list = []
            for tctf in [.1, .3, 1, 3]:
                sim_location = pt.get_sim_location(sim, tctf, beta_list[j], cr, \
                                                   diff = diff_list[j], stream = stream_list[j],
                                                   heat = heat_list[j], work_dir = work_dir)

                sim_base = os.path.basename(sim_location)
                out_name = '../../data/density_contrast_%s_%i'%(sim_base, output)
                if os.path.isfile(out_name) and load == True:
                    mean_rho_cold, mean_rho_hot, median_rho_cold, median_rho_hot = np.loadtxt(out_name, unpack=True)
                    tctf_list.append(tctf)
                    cold_dens_list.append(mean_rho_cold / mean_rho_hot)

                else:
                    ds_loc = '%s/DD%04d/DD%04d'%(sim_location, output, output)
                    if os.path.isfile(ds_loc):
                        ds = yt.load(ds_loc)
                        rho_cold, rho_hot = pt.get_masked_data(ds, 'density')
                        tctf_list.append(tctf)
                        cold_dens_list.append(np.mean(rho_cold) / np.mean(rho_hot))
                        if save:
                            outf = open(out_name, 'w')
                            outf.write("%e %e %e %e\n"%(np.mean(rho_cold), np.mean(rho_hot), np.median(rho_cold), np.median(rho_hot)))
                            outf.close()

            label = pt.get_label_name(compare, tctf, beta_list[j], cr, crdiff = diff_list[j], \
                                  crstream = stream_list[j], crheat = heat_list[j])
                    
            if diff_list[j] > 0:
                label = '$t_{diff} / t_{ff}$ = %i'%diff_list[j]
            if stream_list[j] > 0:
                label = '$\\beta = %i, stream'%beta_list[j]
            if heat_list[j] > 0:
                label += ' + heat'
            if j > 0:
                label = None
            ax.plot(tctf_list, cold_dens_list, color = cpal[i], label = label, linestyle = linestyle_list[j], 
                            linewidth = 2, marker = 'o')

    ax.set_xlabel('P$_c$ / P$_g$')
    ax.set_ylabel('$\\rho_{cold} / \\rho_{hot}$')
    ax.legend(fontsize = 8)
    fig.tight_layout()
    figname = pt.get_fig_name('density_contrast', sim, compare, \
                              tctf, beta, cr, diff_list[0], time = output, \
                              loc = '../../plots/production')
    print(figname)
    plt.savefig(figname, dpi = 300)


work_dir = '../../simulations/production'
load = True
save = True
def make_all_plots(output, compare, beta = 100, cr = 0.1,\
                   tctf = 0.1, crdiff = 0, crstream = 0, crheat = 0):
    for sim in ['isocool', 'isothermal']:
        plot_cold_density(output, sim, compare, tctf, beta, cr, work_dir = work_dir)

                                                                
compare = 'cr'
for output in [30, 50, 80, 100]:
    make_all_plots(output, compare)
