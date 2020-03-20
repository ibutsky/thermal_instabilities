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

import multiprocessing as mp

def plot_cold_fraction_growth(sim, compare, tctf, beta, cr, diff = 0, stream = 0, heat = 0,
                              T_min = 3.3333333e5, zstart = 0.8, zend = 1.2,
                              work_dir = '../../simulations', grid_rank = 3):

    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list\
                        = pt.generate_lists(compare, tctf, crdiff = crdiff, beta = beta, cr = cr)

    print(tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list)

    fig, ax = plt.subplots(figsize = (4.4, 4))
    ax.set_yscale('log')
    ax.set_ylim(5e-3, 5)
    ax.set_xlim(0, 10)
    
    if compare == 'tctf':
        cpal = palettable.cmocean.sequential.Tempo_5_r.mpl_colors
    else:
        cpal = palettable.scientific.sequential.Batlow_8.mpl_colors

    for i, tctf in enumerate(tctf_list):
        time_list, cold_fraction_list = pt.get_time_data('cold_fraction', sim, tctf, beta_list[i], cr_list[i], \
                                           diff = diff_list[i], stream = stream_list[i], heat = heat_list[i],
                                           T_min = T_min, zstart = 0.8, zend = zend, grid_rank = grid_rank,
                                           load = load, save = save, work_dir = work_dir, sim_fam = sim_fam)


        label = pt.get_label_name(compare, tctf, beta_list[i], cr_list[i], crdiff = diff_list[i], \
                           crstream = stream_list[i], crheat = heat_list[i])
        ax.plot(time_list / tctf, cold_fraction_list, linewidth = 2, label = label, color = cpal[i])

     
    ax.set_xlabel('t/t$_{cool}$')
    ax.set_ylabel('Cold Mass Fraction')

    ax.legend()
    fig.tight_layout()
    figname = pt.get_fig_name('cold_fraction_growth', sim, compare, tctf, beta, cr, diff, 
                              sim_fam = sim_fam)
    plt.savefig(figname, dpi = 300)

def make_all_plots(compare, beta = 100, cr = 0, tctf = 0.3):
    all_tctf = [.1, 0.3, 1, 3]
    all_cr = [0.01, .1, 1, 10]
    for sim in ['isocool']:
        if compare == 'diff' or compare == 'stream' or compare == 'transport':
            for tctf in all_tctf:
                for cr in all_cr:
                    plot_cold_fraction_growth(sim, compare, tctf, beta, cr, work_dir = work_dir)

        elif compare == 'cr':
            for tctf in all_tctf:
                    plot_cold_fraction_growth(sim, compare, tctf, beta, cr, work_dir = work_dir)
        else:
            plot_cold_fraction_growth(sim, compare, tctf, beta, cr, work_dir = work_dir)
        

sim_fam = 'production'
work_dir = '../../simulations'
load = True
save = True

crdiff = 0

compare = sys.argv[1]

make_all_plots(compare)
