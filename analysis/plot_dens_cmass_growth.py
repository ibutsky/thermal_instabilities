import yt
from yt import YTArray
from yt import YTQuantity
import sys
import os
import numpy as np
import matplotlib.pylab as plt
import palettable

import plotting_tools as pt
import yt_functions as ytf


def plot_density_fluctuation_growth(sim, compare, tctf, beta, cr, diff = 0, stream = 0, heat = 0,
                                    zstart = 0.8, zend = 1.2, T_cold = 3.33333e5, fs = 12, 
                                    field = 'density', work_dir = '../../simulations/', grid_rank = 3):


    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list \
        = pt.generate_lists(compare, tctf, beta = beta, crdiff = crdiff, cr = cr, crstream = stream, crheat = heat)

    print(tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list)    

    ncols = 3
    fig, ax = plt.subplots(nrows=1, ncols=ncols, figsize = (4*ncols, 3.8), sharex = True, sharey = False)
    for col in range(ncols):
        ax[col].set_yscale('log')
        ax[col].set_xlim(0, 10)
        ax[col].set_xlabel('$t / t_{cool}$', fontsize = fs)

    ax[0].set_ylim(1e-2, 5)
    ax[1].set_ylim(5e-3, 4)
    ax[2].set_ylim(2e-2, 10)
    ax[0].set_ylabel('Density Fluctuation', fontsize = fs)
    ax[1].set_ylabel('Cold Mass Fraction', fontsize = fs)
    ax[2].set_ylabel('Cold Mass Flux', fontsize = fs)
    
    gamma = 5./3.
    time_list = np.arange(0, 12, 1)
    wcool = 1.0 / (gamma * 1.0)
    pi = (5./3.) * wcool
    linecolor = 'black'

    if pt.dark_mode:
        linecolor = 'white'
    if compare == 'tctf' or compare == 'cr':
        pi = 1.0
        ax[0].plot(time_list, 0.02*np.exp(pi*time_list), color = linecolor,\
                   linestyle = 'dashed', label = 'Linear Theory', linewidth = 3)
#    if compare == 'cr':
#        pi = 2./3.
#        ax[0].plot(time_list, 0.02*np.exp(pi*time_list), color = linecolor,\
#                   linestyle = 'dotted', label = 'Linear Theory, $\\eta \\gg 1$', linewidth = 3)

 #       pi = 1./3.
 #       ax[0].plot(time_list, 0.02*np.exp(pi*time_list), color = linecolor,\
 #                  linestyle = 'dotted', label = 'Linear Theory, $\\eta \\gg 1$', linewidth = 3)
#
  #      pi = 1./12.
 #       ax[0].plot(time_list, 0.02*np.exp(pi*time_list), color = linecolor,\
 #                  linestyle = 'dotted', label = 'Linear Theory, $\\eta \\gg 1$', linewidth = 3)

    cpal = pt.get_color_list(compare)

    for col, plot_type in enumerate(['density_fluctuation', 'cold_fraction', 'cold_flux']):
        for i, tctf in enumerate(tctf_list):
            time_list, data_list = pt.get_time_data(plot_type, sim, tctf, beta_list[i], cr_list[i], use_mpi = use_mpi, \
                                           diff = diff_list[i], stream = stream_list[i], heat = heat_list[i], 
                                           field = field, zstart = zstart, zend = zend, grid_rank = grid_rank, 
                                           load = load, save = save, work_dir = work_dir, sim_fam = sim_fam)

            label = pt.get_label_name(compare, tctf, beta_list[i], cr_list[i], crdiff = diff_list[i], \
                                      crstream = stream_list[i], crheat = heat_list[i], counter = i)
#            linestyle = pt.get_linestyle(compare, tctf, beta_list[i], cr_list[i], crdiff = diff_list[i], \
#                                              crstream = stream_list[i], crheat = heat_list[i], counter = i)
            linestyle = 'solid'
            ax[col].plot(time_list/tctf, data_list, linewidth = 3, linestyle = linestyle, label = label, color = cpal[i])
            ax[col].tick_params(labelsize = fs)

            if resolution_compare:
                time_list, data_list = pt.get_time_data(plot_type, sim, tctf, beta_list[i], cr_list[i], use_mpi = use_mpi,\
                                           diff = diff_list[i], stream = stream_list[i], heat = heat_list[i],
                                           field = field, zstart = zstart, zend = zend, grid_rank = grid_rank,
                                           load = load, save = save, work_dir = work_dir, sim_fam = 'production/high_res')
                ax[col].plot(time_list/tctf, data_list, linewidth = 2, linestyle = 'dotted', label = label, color = cpal[i])
            if mhd_compare:
                linestyle_list = ['dashed', 'dotted']
                alpha_list = [.8, .8]
                for j,compare_beta in enumerate([3, 'inf']):
                    time_list, data_list = pt.get_time_data(plot_type, sim, tctf, compare_beta, cr_list[i], use_mpi =use_mpi,\
                                                        diff = diff_list[i], stream = stream_list[i], heat = heat_list[i],
                                           field = field, zstart = zstart, zend = zend, grid_rank = grid_rank,
                                                            load = load, save = save, work_dir = work_dir, sim_fam = sim_fam)
                    label = None
                    ax[col].plot(time_list/tctf, data_list, linewidth = 2, linestyle = linestyle_list[j], 
                                 alpha = alpha_list[j], label = label, color = cpal[i])
            if compare == 'stream' and stream_list[i] >0:
                linestyle = 'dotted'
                time_list, data_list = pt.get_time_data(plot_type, sim, tctf, beta_list[i], cr_list[i], use_mpi =use_mpi,\
                                                        diff = diff_list[i], stream = stream_list[i], heat = 0,
                                           field = field, zstart = zstart, zend = zend, grid_rank = grid_rank,
                                                           load = load, save = save, work_dir = work_dir, sim_fam = sim_fam)
                label = None
                ax[col].plot(time_list/tctf, data_list, linewidth = 2, linestyle = linestyle,
                                 label = label, color = cpal[i])
                    
    

    if compare == 'transport' or compare == 'transport_relative':
        ax[0].legend(fontsize = 7, loc = 3, ncol = 2)
    else:
        ax[1].legend(fontsize = 8, loc = 2)
    fig.tight_layout()
    figname = pt.get_fig_name('dens_cfrac_cflux_growth', sim, compare, \
                              tctf, beta, cr, diff,  sim_fam = sim_fam)
    plt.savefig(figname, dpi = 300)


def make_all_plots(compare, beta = 100, cr = 0, field = 'density', stream = 0, heat = 0):
    all_tctf = [.1, 0.3, 1, 3, 10]
    all_cr = [0.01, .1, 1, 3, 10]
    for sim in ['isocool']:
        if compare == 'diff' or compare == 'stream' or compare == 'transport' or compare == 'transport_relative':
            for tctf in all_tctf:
                for cr in all_cr:
                    plot_density_fluctuation_growth(sim, compare, tctf, beta, cr, work_dir = work_dir, field = field)

        elif compare == 'cr':
            for tctf in all_tctf:
                plot_density_fluctuation_growth(sim, compare, tctf, beta, cr, stream = stream, heat = heat, work_dir = work_dir, field = field)
                
        elif compare == 'tctf':
            tctf = 0.1
            plot_density_fluctuation_growth(sim, compare, tctf, beta, cr, work_dir = work_dir, field = field)
            
        elif compare == 'beta':
            cr = 0
            for tctf in all_tctf:
                plot_density_fluctuation_growth(sim, compare, tctf, beta, cr, work_dir = work_dir, field = field)


sim_fam = 'production/constant_crp'
work_dir = '../../simulations'
save = True
load = True
use_mpi = True
resolution_compare = 0
mhd_compare = 0

crdiff = 0
compare = sys.argv[1]

for beta in ['inf']:#, 100, 10, 3]:
    make_all_plots(compare, beta = beta)#, stream = 0, heat = 0)

