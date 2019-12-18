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


def plot_power_spectrum(sim, tctf_list, beta_list, cr_list, diff_list = [0], stream_list = [0], heat_list = [0], \
                        output = 50, field = 'drho'):
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
    ax.set_xlim(1, 1e2)
    ax.set_ylim(1e-6, 1e2)
    ax.legend(loc = 1)
    fig.tight_layout()
    figname = pt.get_fig_name('power_spectrum', sim, compare, tctf_list[0], beta_list[0], \
                              cr_list[0], diff_list[0])
    plt.savefig(figname, dpi = 300)

    
work_dir = '../../simulations'
sim = sys.argv[1]
compare = sys.argv[2]
tctf = float(sys.argv[3])

output = 80

crdiff = 3.0

tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list = \
                        pt.generate_lists(compare, tctf, crdiff = crdiff)

print(tctf_list, beta_list, cr_list, diff_list)
plot_power_spectrum(sim, tctf_list, beta_list, cr_list, diff_list = diff_list, output = output,\
                    stream_list = stream_list, heat_list = heat_list)
