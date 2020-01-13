import yt
from yt import YTArray
from yt import YTQuantity

import numpy as np
import sys
import os
import glob

import matplotlib.pylab as plt
from matplotlib.colors import SymLogNorm, LogNorm

import palettable
import plotting_tools as pt
import yt_functions as ytf


def plot_power_spectrum(sim, compare, tctf, beta, cr, diff = 0, stream = 0, heat = 0,
                        output = 50, field = 'drho', work_dir = '../../simulations/production'):

    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list\
                        = pt.generate_lists(compare, tctf, crdiff = diff, cr = cr)
    print(tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list)

    color_list  = palettable.scientific.sequential.Batlow_6.mpl_colors
    fig, ax = plt.subplots(figsize = (6, 6))
    for i, tctf in enumerate(tctf_list):
        color = color_list[i]
        # load the simulation
        sim_loc = pt.get_sim_location(sim, tctf, beta_list[i], cr_list[i], \
                                           diff = diff_list[i], stream = stream_list[i],
                                           heat = heat_list[i], work_dir = work_dir)
        ds_path = '%s/DD%04d/DD%04d'%(sim_loc, output, output)
        if not os.path.isfile(ds_path):
            print('nope')
            continue
        ds = yt.load(ds_path)
        
        label = pt.get_label_name(compare, tctf, beta_list[i], cr_list[i], crdiff = diff_list[i], \
                       crstream = stream_list[i], crheat = heat_list[i])
        k, P_k = pt.make_power_spectrum(ds, field)
        ax.loglog(k, P_k, label = label, linewidth = 3, color = color)
     
    ax.set_xlabel('k')
    ax.set_ylabel('P(k)')
    ax.set_xlim(1, 200)
    ax.set_ylim(1e-6, 1e2)
    ax.legend(loc = 1)
    fig.tight_layout()
    figname = pt.get_fig_name('power_spectrum', sim, compare, tctf_list[0], beta_list[0], \
                              cr_list[0], diff_list[0], loc = '../../plots/production')
    plt.savefig(figname, dpi = 300)


def make_all_plots(compare, beta = 100, cr = 0.1, output_list = [50]):
    all_tctf = [.1, 0.3, 1, 3]
    all_cr = [0.01, .1, 1, 10]
    for sim in ['isothermal', 'isocool']:
        for output in output_list:
            if compare == 'diff' or compare == 'stream':
                for tctf in all_tctf:
                    for cr in all_cr:
                        plot_power_spectrum(sim, compare, tctf, beta, cr, output = output, work_dir = work_dir)

            elif compare == 'cr':
                for tctf in all_tctf:
                    plot_power_spectrum(sim, compare, tctf, beta, cr, output = output, work_dir = work_dir)
    
work_dir = '../../simulations/production'

compare = sys.argv[1]

make_all_plots(compare, output_list  = [30, 50, 80])
