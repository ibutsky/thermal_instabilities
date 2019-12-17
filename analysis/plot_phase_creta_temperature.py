import yt

import numpy as np
import glob
import os
import sys

import matplotlib.pylab as plt
from matplotlib.colors import LogNorm, SymLogNorm
import palettable
import astropy.constants as const

import yt_functions as ytf
import plotting_tools as pt

half_range = 1
rho0 = 1e-27
T0 = 1e6
mu = 1.22
mh = const.m_p.cgs.value
kb = const.k_B.cgs.value
p0 = (rho0 / mu / mh) * kb*T0

workdir = '../../simulations'
sim = 'isothermal'
tctf = float(sys.argv[1])
cr_list = [0.01, 0.1, 1.0, 10.0]
beta_list = [10, 10, 10, 10]
diff = 3.0

title_list = []
for cr, beta in zip(cr_list, beta_list):
    if beta == 'inf':
        title_list.append('Hydro')
    else:
        if cr == 0:
            title_list.append('$\\beta = $ %.1f'%beta)
        elif cr == 0.01:
            title_list.append('P$_c$/P$_g$ = %.2f'%cr)
        else:
            title_list.append('P$_c$/P$_g$ = %.1f'%cr)

print(title_list)
    
def plot_phase(output, folder = '.'):

    cmap_list = [palettable.cmocean.sequential.Tempo_20.mpl_colormap]
    ncols = int(len(cr_list) / 2)
    nrows = 2
    fig, ax = plt.subplots(ncols = ncols, nrows = nrows, figsize=(4*ncols, 4.4*nrows))

    for i, cr in enumerate(cr_list):
        sim_loc = pt.get_sim_location(sim, tctf, beta_list[i], cr, diff = diff)
        ds = ytf.load('%s/DD%04d/DD%04d'%(sim_loc, output, output))

        ad = ds.all_data()
        ph = yt.PhasePlot(ad, ('gas', 'temperature'), ('gas', 'cr_eta'), ('gas', 'cell_mass'),\
                          weight_field = None, fractional = True)

        ph.set_ylim(cr/ 1e3, cr * 1e3)
        prof = ph.profile

        xbins = prof.x
        ybins = prof.y
        data  = prof[('gas', 'cell_mass')].T


        vmin = 1e-5
        vmax = 1e-2
        cmap = cmap_list[0]
        cmap = palettable.cubehelix.jim_special_16_r.mpl_colormap
        
        row = int(i / ncols)
        col = i - row*ncols
        print(i, row, col)
        ax[row][col].set_xscale('log')
        ax[row][col].set_yscale('log')
#        ax[row][col].set_xlim(1.5e-28, 8e-26)
        ax[row][col].set_xlim(4e4, 1e7)

        pcm = ax[row][col].pcolormesh(xbins, ybins, data, norm = LogNorm(), cmap = cmap, \
                               vmax = vmax, vmin = vmin)

        ax[row][col].set_title(title_list[i], fontsize = 16)
#        ax[row][col].set_ylabel('Density (g cm$^{-3}$)')
        ax[row][col].set_ylabel('CR Pressure / Gas Pressure')
        ax[row][col].set_xlabel('Temperature (K)')

    
    figname = '../../plots/phase_temperature_cr_eta_tctf_%.1f_%.1f.png'%(tctf, output)
    if diff > 0: 
        figname = '../../plots/phase_temperature_cr_eta_mass_tctf_%.1f_cr_diff_%.1f_%.1f.png'%(tctf, diff, output)
    fig.tight_layout()
    plt.savefig(figname, dpi = 300)



                    
output = 80
plot_phase(output = output)
