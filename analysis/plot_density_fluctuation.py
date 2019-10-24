import yt
from yt import YTArray
from yt import YTQuantity
import sys
import os
import numpy as np
import matplotlib.pylab as plt
import palettable

import plotting_tools as pt

def calculate_drho_rms_list(model, beta, tctf_list, work_dir = '.'):
    drho_rms_list = []
    for tctf in tctf_list:
        if beta == 'inf':
            sim_location = '%s/%s_tctf_%.1f'%(work_dir, model, tctf)
        else:
            sim_location = '%s/%s_tctf_%.1f_beta_%.1f'%(work_dir, model, tctf, beta)
        if os.path.isdir(sim_location):    
            drho_rms_list.append(pt.calculate_averaged_drho_rms(output_list, sim_location))
    return drho_rms_list

def plot_density_fluctuation(model, beta_list = ['inf'], tctf_list = None, output_list = [50], \
                                    work_dir = '../../simulations', compare = False):
    if tctf_list == None:
        tctf_list = [0.1, 0.3, 1.0, 3.0, 10.0]

    fig, ax = plt.subplots(figsize = (6, 6))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.08, 15)
    ax.set_ylim(1e-2, 10)
    
#    cpal = palettable.wesanderson.Darjeeling4_5.mpl_colors
    cpal = palettable.cmocean.sequential.Deep_4_r.mpl_colors

    for i, beta in enumerate(beta_list):
        if beta == 'inf':
            label = 'Hydro'
        else:
            label = '$\\beta = %i$'%(beta)
        drho_rms_list = calculate_drho_rms_list(model, beta, tctf_list, work_dir = work_dir)
        ax.plot(tctf_list, drho_rms_list, color = cpal[i], label = label, linewidth = 3, marker = 'o')
        
        if beta == 'inf' and compare == True:
            sim_location = '%s/%s_64'%(work_dir, model)
            drho_rms = pt.calculate_averaged_drho_rms(output_list, sim_location)
            ax.scatter(1, drho_rms, marker = 'v', color = 'black', label = 'res = 64$^3$')

            sim_location = '%s/%s_256'%(work_dir, model)
            drho_rms = pt.calculate_averaged_drho_rms(output_list, sim_location)
            ax.scatter(1, drho_rms, marker = '^', color = 'black', label = 'res = 256$^3$')

    ax.set_xlabel('t$_{cool}$ / t$_{ff}$')
    ax.set_ylabel('$ \\langle \\delta \\rho / \\rho\\rangle_{\\mathrm{rms}}$ ')
    ax.legend()
    fig.tight_layout()
    plt.savefig('../../plots/density_fluctuation_tcool_%.1f_%s.png'%(np.mean(output_list)/10, model), dpi = 300)


tctf_list = [0.1, 0.3, 1.0, 10]
    
model = sys.argv[1]
output_list = [90, 95, 100]
output_list = [40, 45, 50, 55, 60]
beta_list = [4, 100, 'inf']
plot_density_fluctuation(model, output_list = output_list, beta_list = beta_list)
