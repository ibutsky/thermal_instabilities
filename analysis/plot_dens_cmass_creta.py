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


def plot_density_fluctuation(output, sim, compare, tctf, beta, cr, diff = 0, stream = 0, heat = 0,
                              T_cold = 3.3333333e5, zstart = 0.8, zend = 1.2, relative = 0, 
                              work_dir = '../../simulations/', grid_rank = 3):

    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list\
                        = pt.generate_lists(compare, tctf, crdiff = diff, cr = cr, beta = beta, cr_only = 1)
    tctf_list = [0.1, 0.3, 1, 3]
    all_cr_list = [0.01, 0.1, 1, 10]
    print(tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list)
    
    fig, ax = plt.subplots(nrows=1, ncols = 3, figsize = (4.4*3, 4), sharex = True, sharey = False)
    for col in range(3):
        ax[col].set_xscale('log')
        ax[col].set_yscale('log')
        ax[col].set_xlim(5e-3, 5e2)
        ax[col].set_xlabel('$ P_c / P_g $')

    if relative == 0:
        ax[0].set_ylim(1e-2, 5)
        ax[1].set_ylim(5e-3, 4)
        ax[2].set_ylim(5e-3, 4)
        ax[0].set_ylabel('Density Fluctuation')
        ax[1].set_ylabel('Cold Mass Fraction')
        ax[2].set_ylabel('Cold Mass Flux')
    else:
        ax[0].set_ylim(1e-2, 10)
        ax[1].set_ylim(1e-2, 100)
        ax[2].set_ylim(1e-2, 10)
        ax[0].set_ylabel('Relative Density Fluctuation')
        ax[1].set_ylabel('Relative Cold Mass Fraction')
        ax[2].set_ylabel('Relative Cold Mass Flux') 

    color_list  = pt.get_color_list('tctf')
    marker = 'o'

    for acr in all_cr_list:
        cr_list = len(beta_list)*[acr]
        for i in range(len(cr_list)):
            for col, plot_type in enumerate(['density_fluctuation', 'cold_fraction', 'cold_flux']):
                for j, tctf in enumerate(tctf_list):
                    data = 0
                    creta = 0
                    err = 0
                    time_list, data_list = pt.get_time_data(plot_type, sim, tctf, beta_list[i], cr_list[i], \
                                           diff = diff_list[i], stream = stream_list[i], heat = heat_list[i],
                                           T_min = T_cold, zstart = zstart, zend = zend, grid_rank = grid_rank,
                                           load = load, save = save, work_dir = work_dir, sim_fam = sim_fam)
                    creta_time_list, creta_list = pt.get_time_data('cold_creta', sim, tctf, beta_list[i], cr_list[i], \
                                           diff = diff_list[i], stream = stream_list[i], heat = heat_list[i],
                                           T_min = T_cold, zstart = zstart, zend = zend, grid_rank = grid_rank,
                                           load = load, save = save, work_dir = work_dir, sim_fam = sim_fam)
                    if len(data_list) > 0 and len(creta_list) > 0:                    
                        data = np.nan_to_num(data_list[output-10:output+10])
                        creta = np.nan_to_num(creta_list[output-10:output+10])                    
                    
                    if acr == 0.01 and i == 0 and col == 0:
                        label = pt.get_label_name('tctf', tctf, beta_list[i], cr_list[i], crdiff = diff_list[i], \
                                                  crstream = stream_list[i], crheat = heat_list[i], counter = i)
                    else:
                        label = None
                    ax[col].scatter(np.mean(creta), np.mean(data), color = color_list[j], label = label, 
                                        marker = marker)
                    ax[col].errorbar(np.mean(creta), np.mean(data), xerr = np.std(creta), yerr = np.std(data),
                                     color = color_list[j])

    ax[0].legend(fontsize = 8)
    fig.tight_layout()
    fig_basename = 'dens_cfrac_cflux_creta'
    figname = pt.get_fig_name(fig_basename, sim, compare, \
                              tctf, beta, cr, crdiff = diff, crstream = stream, \
                              crheat = heat, time = output, sim_fam = sim_fam,\
                              loc = '../../plots')
    print(figname)
    plt.savefig(figname, dpi = 300)

sim_fam = 'production'
work_dir = '../../simulations'
load = True
save = True

def make_all_plots(output, compare, beta = 100, cr = 0.1,\
                   tctf = 0.1, crdiff = 0, crstream = 0, crheat = 0):
    for sim in ['isocool']:
        plot_density_fluctuation(output, sim, compare, tctf, beta, cr, work_dir = work_dir)

                                                                

sim = 'isocool'
tctf = .1

stream = 0
heat = 0
diff = 0
beta = 100

relative = 0
#for compare in ['stream', 'diff']:
#    for cr in [.1, 1, 10]:
#        for output in [40]:

for compare in ['transport']:
    for cr in [0.01, 0.1, 1, 10]:
#    for cr in [1.0]:
        for output in [50]:
            plot_density_fluctuation(output, sim, compare, tctf, beta, cr, diff = diff, stream = stream, heat = heat, \
                                     work_dir = work_dir, relative = relative)

