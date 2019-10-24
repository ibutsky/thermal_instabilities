import yt
from yt import YTArray
from yt import YTQuantity
import sys
import os
import numpy as np
import matplotlib.pylab as plt
import palettable

import plotting_tools as pt


def calculate_cold_fraction_list(model, beta, tctf_list, work_dir = '.', Tmin = 3.333333e5):
    cold_fraction_list = np.array([])
    for tctf in tctf_list:
        if beta == 'inf':
            sim_location = '%s/%s_tctf_%.1f'%(work_dir, model, tctf)
        else:
            sim_location = '%s/%s_tctf_%.1f_beta_%.1f'%(work_dir, model, tctf, beta)
        if os.path.isdir(sim_location):
            cold_fraction_list = np.append(cold_fraction_list, 
                        pt.calculate_averaged_cold_fraction(output_list, sim_location, Tmin = Tmin))
    return cold_fraction_list

def plot_cold_fraction(model, beta_list = ['inf'], tctf_list = None, output_list = [50], compare = False,\
                                    work_dir = '../../simulations', Tmin = 3.3333333e5):
    if tctf_list == None:
        tctf_list = np.array([0.1, 0.3, 1.0, 10.0])
    else:
        tctf_list = np.array(tctf_list)

    fig, ax = plt.subplots(figsize = (6, 6))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.08, 15)
    ax.set_ylim(1e-2, 1)
    
    
    cpal = palettable.cmocean.sequential.Deep_4_r.mpl_colors
    for i, beta in enumerate(beta_list):
        if beta == 'inf':
            label = 'Hydro'
        else:
            label = '$\\beta = %i$'%(beta)
        cold_fraction_list = calculate_cold_fraction_list(model, beta, tctf_list, work_dir = work_dir, Tmin = Tmin)
        mask = cold_fraction_list > 0
        ax.plot(tctf_list[mask], cold_fraction_list[mask], color = cpal[i], label = label, linewidth = 3, marker = 'o')

        if beta == 'inf' and compare == True:
            # add in a couple measurements by hand
            sim_location = '%s/%s_64'%(work_dir, model)
            cold_fraction = pt.calculate_averaged_cold_fraction(output_list, sim_location, Tmin = Tmin)
            ax.scatter(1, cold_fraction, marker = 'v', color = cpal[i], label = 'res = 64$^3$')

            sim_location = '%s/%s_256'%(work_dir, model)
            cold_fraction = pt.calculate_averaged_cold_fraction(output_list, sim_location, Tmin = Tmin)
            ax.scatter(1, cold_fraction, marker = '^', color = cpal[i], label= 'res = 256$^3$')


    ax.set_xlabel('t$_{cool}$ / t$_{ff}$')
    ax.set_ylabel('Cold Fraction')
    ax.legend()
    fig.tight_layout()
    plt.savefig('../../plots/cold_fraction_tcool_%.1f_%s.png'%(np.mean(output_list)/10, model), dpi = 300)


tctf_list = [0.1, 0.3, 1.0, 3.0, 10]
beta_list = [4, 100, 'inf']
                                                                
Tmin = 1e6 / 3.0
model = sys.argv[1]
output_list = [90, 95, 100]
#output_list = [40, 45, 50, 55, 60]
plot_cold_fraction(model, output_list = output_list, \
                   beta_list = beta_list, tctf_list = tctf_list, Tmin = Tmin)
