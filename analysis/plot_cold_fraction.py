import yt
from yt import YTArray
from yt import YTQuantity
import sys
import os
import numpy as np
import matplotlib.pylab as plt
import palettable

import plotting_tools as pt


def calculate_cold_fraction_list(model, beta, tctf_list, work_dir = '.', T_min = 3.333333e5, z_min = 0.1):
    cold_fraction_list = np.array([])
    for tctf in tctf_list:
        if beta == 'inf':
            sim_location = '%s/%s_tctf_%.1f'%(work_dir, model, tctf)
        else:
            sim_location = '%s/%s_tctf_%.1f_beta_%.1f'%(work_dir, model, tctf, beta)
        if os.path.isdir(sim_location):
            cold_fraction_list = np.append(cold_fraction_list, \
                                           pt.calculate_averaged_cold_fraction(output_list, \
                                           sim_location, T_min = T_min, z_min = z_min))
    return cold_fraction_list

def plot_cold_fraction(model, beta_list = ['inf'], tctf_list = None, output_list = [50], compare = False,\
                                    work_dir = '../../simulations', T_min = 3.3333333e5, z_min = 0.1):
    if tctf_list == None:
        tctf_list = np.array([0.1, 0.3, 1.0, 10.0])
    else:
        tctf_list = np.array(tctf_list)

    fig, ax = plt.subplots(figsize = (6, 6))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.08, 15)
    ax.set_ylim(1e-2, 1)
    
    
    cpal = palettable.cmocean.sequential.Deep_7_r.mpl_colors
    for i, beta in enumerate(beta_list):
        if beta == 'inf':
            label = 'Hydro'
        else:
            label = '$\\beta = %i$'%(beta)
        cold_fraction_list = calculate_cold_fraction_list(model, beta, tctf_list, \
                                                          work_dir = work_dir, T_min = T_min, z_min = z_min)
        mask = cold_fraction_list > 0
        ax.plot(tctf_list[mask], cold_fraction_list[mask], color = cpal[i], label = label, linewidth = 3, marker = 'o')

        if beta == 'inf' and compare == True:
            # add in a couple measurements by hand
            sim_location = '%s/%s_64'%(work_dir, model)
            cold_fraction = pt.calculate_averaged_cold_fraction(output_list, sim_location, T_min = T_min)
            ax.scatter(1, cold_fraction, marker = 'v', color = cpal[i], label = 'res = 64$^3$')

            sim_location = '%s/%s_256'%(work_dir, model)
            cold_fraction = pt.calculate_averaged_cold_fraction(output_list, sim_location, T_min = T_min)
            ax.scatter(1, cold_fraction, marker = '^', color = cpal[i], label= 'res = 256$^3$')


    ax.set_xlabel('t$_{cool}$ / t$_{ff}$')
    ax.set_ylabel('Cold Fraction')
    ax.legend()
    fig.tight_layout()
    if model.__contains__('no_center'):
        model = 'no_center_heating_isothermal'
    figname = '../../plots/cold_fraction_tcool_%.1f_%s'%(np.mean(output_list)/10, model)
    print(compare)
    if compare == 'beta':
        figname += '_beta_compare'
    print(figname)
    plt.savefig(figname+'.png', dpi = 300)

compare = 'beta'
if compare == 'beta':
    tctf_list = [0.1, 0.3, 1.0, 3.0, 10]
    beta_list = [3, 10, 30, 100, 300, 'inf']
    cr_list = 6*[0]
                                                                
T_min = 1e6 / 3.0
z_min = 0.1
#model = 'no_center_heating/isothermal'
model = 'isothermal'
output_list = [90,92, 95,97, 100]
output_list = [80, 85, 90, 95, 100]
#output_list = [40, 45, 50, 55, 60]
#output_list = [20, 25, 30, 35, 40]
output_list = [25, 30, 35]
plot_cold_fraction(model, output_list = output_list, \
                   beta_list = beta_list, tctf_list = tctf_list, T_min = T_min, z_min = z_min)
