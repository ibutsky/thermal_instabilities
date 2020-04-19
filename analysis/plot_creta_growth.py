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
                                    zstart = 0.8, zend = 1.2, T_cold = 3.33333e5, 
                                    field = 'density', work_dir = '../../simulations/', grid_rank = 3):


    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list \
                        = pt.generate_lists(compare, tctf, beta = beta, crdiff = crdiff, cr = cr)

    mask = cr_list > 0
    tctf_list = tctf_list[mask]
    beta_list = beta_list[mask]
    cr_list = cr_list[mask]
    diff_list = diff_list[mask]
    stream_list = stream_list[mask]
    heat_list = heat_list[mask]

    ncols = 1
    fig, ax = plt.subplots(nrows=1, ncols=ncols, figsize = (4.4*ncols, 4), sharex = True, sharey = False)

    ax.set_yscale('log')
    ax.set_xlim(0, 10)
    ax.set_xlabel('$t / t_{cool}$')

    ax.set_ylim(1e-2, 1e4)
    ax.set_ylabel('$ \\langle P_{c} / P_{g} \\rangle $')
    
    cpal = pt.get_color_list(compare)

    for i, tctf in enumerate(tctf_list):
        time_list, data_list = pt.get_time_data('cold_creta', sim, tctf, beta_list[i], cr_list[i], \
                                           diff = diff_list[i], stream = stream_list[i], heat = heat_list[i], 
                                           field = field, zstart = zstart, zend = zend, grid_rank = grid_rank, 
                                           load = load, save = save, work_dir = work_dir, sim_fam = sim_fam)

        label = pt.get_label_name(compare, tctf, beta_list[i], cr_list[i], crdiff = diff_list[i], \
                                      crstream = stream_list[i], crheat = heat_list[i], counter = i)
        linestyle = pt.get_linestyle(compare, tctf, beta_list[i], cr_list[i], crdiff = diff_list[i], \
                                              crstream = stream_list[i], crheat = heat_list[i], counter = i)
            
        ax.plot(time_list/tctf, data_list, linewidth = 3, linestyle = linestyle, label = label, color = cpal[i])

     
    ax.legend()
    fig.tight_layout()
    figname = pt.get_fig_name('cold_creta_growth', sim, compare, \
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

crdiff = 0
compare = sys.argv[1]



make_all_plots(compare)
