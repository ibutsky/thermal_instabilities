import yt
from yt import YTArray
from yt import YTQuantity
import sys
import os
import numpy as np
import matplotlib.pylab as plt
import palettable

import plotting_tools as pt


def calculate_cold_fraction_list(sim, beta, cr, diff = 0, T_min = 3.333333e5, z_min = 0.1,\
                                 tctf_list = [0.1, 0.3, 1.0, 3.0, 10.0], work_dir = '../../simulations'):
    cold_fraction_list = np.array([])
    for i, tctf in enumerate(tctf_list):
        sim_location = pt.get_sim_location(sim, tctf, beta, cr, diff = diff, work_dir = work_dir)
        if os.path.isdir(sim_location):
            cold_fraction_list = np.append(cold_fraction_list, \
                                           pt.calculate_averaged_cold_fraction(output_list, \
                                           sim_location, T_min = T_min, z_min = z_min))
    return cold_fraction_list

def plot_cold_fraction(model, beta_list = ['inf'], tctf_list = None, output_list = [50], compare = False,\
                              cr_list = None, diff_list = None,  work_dir = '../../simulations', T_min = 3.3333333e5, z_min = 0.1):
    if tctf_list == None:
        tctf_list = np.array([0.1, 0.3, 1.0, 3.0, 10.0])
    else:
        tctf_list = np.array(tctf_list)

    fig, ax = plt.subplots(figsize = (6, 6))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.08, 15)
    ax.set_ylim(1e-2, 1)
    
    
    cpal = palettable.cmocean.sequential.Deep_7_r.mpl_colors
    for i, beta in enumerate(beta_list):
        cold_fraction_list = calculate_cold_fraction_list(model, beta, cr_list[i], diff_list[i], tctf_list = tctf_list, \
                                                          work_dir = work_dir, T_min = T_min, z_min = z_min)
        label = pt.get_label_name(compare, tctf_list[0], beta, cr_list[i])
        mask = cold_fraction_list > 0
        ax.plot(tctf_list[mask], cold_fraction_list[mask], color = cpal[i], label = label, linewidth = 3, marker = 'o')

    
    ax.set_xlabel('t$_{cool}$ / t$_{ff}$')
    ax.set_ylabel('Cold Fraction')
    ax.legend()
    fig.tight_layout()

    time = np.mean(output_list)
    figname = pt.get_fig_name('cold_fraction', sim, compare, tctf_list[0], beta_list[0], \
                              cr_list[0], diff_list[0], time = time)
    plt.savefig(figname, dpi = 300)

sim = sys.argv[1]
compare = sys.argv[2]
time = int(sys.argv[3])
tctf = 0.1
crdiff = 0

tctf_list, beta_list, cr_list, diff_list = pt.generate_lists(compare, tctf, crdiff = crdiff)
                                                                
T_min = 1e6 / 3.0
z_min = 0.1
output_list = [time]
plot_cold_fraction(sim, output_list = output_list, beta_list = beta_list, \
                   cr_list = cr_list, diff_list = diff_list, \
                   T_min = T_min, z_min = z_min, compare = compare)
