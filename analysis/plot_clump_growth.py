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
        = pt.generate_lists(compare, tctf, beta = beta, crdiff = crdiff, cr = cr)

    print(tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list)
    
    ncols = 2
    fig, ax = plt.subplots(nrows=1, ncols=ncols, figsize = (4*ncols, 3.8), sharex = True, sharey = False)
    for col in range(ncols):
        ax[col].set_yscale('log')
        ax[col].set_xlim(0, 10)
        ax[col].set_xlabel('$t / t_{cool}$', fontsize = fs)

    ax[0].set_ylim(1, 100)
    ax[1].set_ylim(1, 300)
    ax[0].set_ylabel('Clump Size', fontsize = fs)
    ax[1].set_ylabel('Number of Clumps', fontsize = fs)

    cpal = pt.get_color_list(compare)

    for i, tctf in enumerate(tctf_list):
        time_list, clump_data = pt.get_time_data('clump', sim, tctf, beta_list[i], cr_list[i], use_mpi = True,\
                                           diff = diff_list[i], stream = stream_list[i], heat = heat_list[i], 
                                           field = field, zstart = zstart, zend = zend, grid_rank = grid_rank, 
                                           load = load, save = save, work_dir = work_dir, sim_fam = sim_fam)
#        print(clump_data[0])
#        clump_data = list(clump_data)
        print(len(clump_data), len(clump_data[0]))
#        n_clumps, clump_size, clump_std = list(zip(*clump_data))
        n_clumps, clump_size, clump_std = clump_data
        print(n_clumps, clump_size)

        label = pt.get_label_name(compare, tctf, beta_list[i], cr_list[i], crdiff = diff_list[i], \
                                      crstream = stream_list[i], crheat = heat_list[i], counter = i)

        linestyle = 'solid'
        ax[0].plot(time_list/tctf, clump_size, linewidth = 3, linestyle = linestyle, label = label, color = cpal[i])
        ax[0].tick_params(labelsize = fs)
        ax[1].plot(time_list/tctf, n_clumps, linewidth = 3, linestyle = linestyle, color = cpal[i])
    
     
    ax[0].legend()
    fig.tight_layout()
    figname = pt.get_fig_name('clump_growth', sim, compare, \
                              tctf, beta, cr, diff,  sim_fam = sim_fam)
    plt.savefig(figname, dpi = 300)


def make_all_plots(compare, beta = 100, cr = 0, field = 'density'):
    all_tctf = [.1, 0.3, 1, 3, 10]
    all_cr = [0.01, .1, 1, 10]
    for sim in ['isocool']:
        if compare == 'diff' or compare == 'stream' or compare == 'transport':
            for tctf in all_tctf:
                for cr in all_cr:
                    plot_density_fluctuation_growth(sim, compare, tctf, beta, cr, work_dir = work_dir, field = field)

        elif compare == 'cr':
            for tctf in all_tctf:
                plot_density_fluctuation_growth(sim, compare, tctf, beta, cr, work_dir = work_dir, field = field)
                
        elif compare == 'tctf':
            tctf = 0.1
            plot_density_fluctuation_growth(sim, compare, tctf, beta, cr, work_dir = work_dir, field = field)
            
        elif compare == 'beta':
            cr = 0
            for tctf in all_tctf:
                plot_density_fluctuation_growth(sim, compare, tctf, beta, cr, work_dir = work_dir, field = field)


sim_fam = 'production'
work_dir = '../../simulations'
save = True
load = True
resolution_compare = 0
mhd_compare = 0

crdiff = 0
compare = sys.argv[1]


make_all_plots(compare, cr = 1)
