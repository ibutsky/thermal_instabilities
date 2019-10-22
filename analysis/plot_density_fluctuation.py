import yt
from yt import YTArray
from yt import YTQuantity
import sys
import numpy as np
import matplotlib.pylab as plt
import palettable

import plotting_tools as pt


def plot_density_fluctuation(model, beta_list = None, tctf_list = None, output_list = [50], \
                                    work_dir = '../../simulations'):
    if tctf_list == None:
        tctf_list = [0.1, 0.3, 1.0, 10.0]

    fig, ax = plt.subplots(figsize = (6, 6))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.08, 15)
    ax.set_ylim(1e-2, 10)
    
    cpal = palettable.wesanderson.Darjeeling4_5.mpl_colors

    drho_rms_list = []
    for tctf in tctf_list:        
        sim_location = '%s/%s_tctf_%.1f'%(work_dir, model, tctf)
        drho_rms_list.append(pt.calculate_averaged_drho_rms(output_list, sim_location))

    ax.plot(tctf_list, drho_rms_list, color = 'black', label = 'Hydro', marker = 'o')

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
    plt.savefig('../../plots/density_fluctuation_%s.png'%model, dpi = 300)


tctf_list = [0.1, 0.3, 1.0, 10]
    
model = sys.argv[1]
output_list = [90, 95, 100]
plot_density_fluctuation(model, output_list = output_list)
