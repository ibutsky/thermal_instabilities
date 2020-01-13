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


    
def plot_phase(output, sim, tctf, beta_list, cr_list, diff_list, stream_list, heat_list, title_list):

    cmap_list = [palettable.cmocean.sequential.Tempo_20.mpl_colormap]
    ncols = int(len(cr_list) / 2)
    nrows = 2
    fig, ax = plt.subplots(ncols = ncols, nrows = nrows, figsize=(4*ncols, 4.4*nrows))

    for i, cr in enumerate(cr_list):
        sim_loc = pt.get_sim_location(sim, tctf, beta_list[i], cr, \
                                      diff = diff_list[i], stream = stream_list[i], heat = heat_list[i], work_dir = workdir)
        ds = ytf.load('%s/DD%04d/DD%04d'%(sim_loc, output, output))

        ad = ds.all_data()
        ad_cut = ad.cut_region(["(obj[('gas', 'z')].in_units('kpc') > 4.3) | (obj[('gas', 'z')].in_units('kpc') < -4.3)"])
        ph = yt.PhasePlot(ad_cut, ('gas', 'density'), ('gas', 'temperature'), ('gas', 'cell_mass'),\
                          weight_field = None, fractional = True)
        ph.save()
        prof = ph.profile
        xbins = prof.x
        ybins = prof.y
        data  = prof[('gas', 'cell_mass')].T

        print(xbins, ybins, data)
        vmin = 1e-5
        vmax = 1e-2
        cmap = cmap_list[0]
        cmap = palettable.cubehelix.jim_special_16_r.mpl_colormap
        
        row = int(i / ncols)
        col = i - row*ncols
        print(i, row, col)
        ax[row][col].set_xscale('log')
        ax[row][col].set_yscale('log')
        ax[row][col].set_xlim(1.5e-28, 8e-26)
        ax[row][col].set_ylim(4e4, 1e7)

        pcm = ax[row][col].pcolormesh(xbins, ybins, data, norm = LogNorm(), cmap = cmap, \
                               vmax = vmax, vmin = vmin)

        ax[row][col].set_title(title_list[i], fontsize = 16)
        ax[row][col].set_xlabel('Density (g cm$^{-3}$)')
        ax[row][col].set_ylabel('Temperature (K)')

    
    figname = '../../plots/production/phase_density_temperature_mass_%s_tctf_%.1f_cr_%.2f_%.1f.png'%(sim, tctf, cr, output)

    fig.tight_layout()
    print(figname)
    plt.savefig(figname, dpi = 300)



                    
workdir = '../../simulations/production'

output_list = [30, 50, 80]
for sim in ['isothermal', 'isocool']:
    for tctf in [.1, .3, 1, 3]:
        for cr in [.01, .1, 1, 10]:
            cr_list = [0, cr, cr, cr, cr, cr]
            beta_list = len(cr_list) * [100]
            diff_list = [0, 0, 3.0, 1.0, 0, 0]
            stream_list = [0, 0, 0, 0, 1, 1]
            heat_list = [0, 0, 0, 0, 0, 1]

            title_list = ['No CR', 'P$_c$/P$_g$ = %.2f'%cr, 'P$_c$/P$_g$ = %.2f, tdiff = 3tff'%cr, \
              'P$_c$/P$_g$ = %.2f, tdiff = tff'%cr, 'P$_c$/P$_g$ = %.2f + stream'%cr, \
                          'P$_c$/P$_g$ = %.2f + stream + heat'%cr]

            for output in output_list:
                plot_phase(output, sim, tctf, beta_list, cr_list, diff_list, stream_list, heat_list, title_list)
