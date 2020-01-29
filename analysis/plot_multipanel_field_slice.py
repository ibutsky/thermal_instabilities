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


def plot_multipanel_slices(field, output, ds_loc_list, label_list, folder = '.'):

    cmap_list = [palettable.cmocean.sequential.Tempo_20.mpl_colormap]
    fig, ax = plt.subplots(ncols = len(ds_loc_list), nrows = 1, figsize=(1.2*len(ds_loc_list), 3.8))

    for i, ds_loc in enumerate(ds_loc_list):
        print(ds_loc)
        if not os.path.isfile('%s/DD%04d/DD%04d'%(ds_loc, output, output)):
            continue
        ds = ytf.load('%s/DD%04d/DD%04d'%(ds_loc, output, output))

        s = yt.SlicePlot(ds, 'x', ('gas', field))
        frb = s.frb

        xbins = frb['y'].in_units('kpc')
        ybins = frb['z'].in_units('kpc')
        if field == 'density':
            data_norm = rho0
        else:
            data_norm = 1

        data   = frb[field] / data_norm
        
        if field == 'density':
            vmin = 3e-2
            vmax = 1
        elif field == 'temperature':
            vmin = 5e4
            vmax = 5e6
        elif field == 'cr_eta':
            cr_eta = pt.get_cr_eta(ds)
            vmin = cr_eta / 100
            vmax = cr_eta * 100

        cmap = pt.get_cmap(field)
        pcm = ax[i].pcolormesh(xbins, ybins, data, norm = LogNorm(), cmap = cmap, \
                               vmax = vmax, vmin = vmin)


        ax[i].set_aspect('equal')
        ax[i].set_xticks([])
        ax[i].set_yticks([])
        ax[i].set_title(label_list[i], fontsize = 8)

    
    figname = pt.get_fig_name('%s_multipanel_slice'%field, sim, compare, \
                           tctf, cr=cr, time = output, loc = '../../plots/production')
    plt.subplots_adjust(wspace=0.02, top = 0.9)
    plt.savefig(figname, dpi = 300)


def get_iteration_lists(compare, tctf = 0.1, beta = 100, cr = 0.1):
    all_tctf = [.1, 0.3, 1, 3]
    all_cr = [0, 0.01, .1, 1, 10]
    all_sim = ['isothermal', 'isocool']
    if compare == 'cr':
        all_cr = [cr]  #won't actually matter
    elif compare == 'tctf':
        all_tctf = [tctf]

    return all_sim, all_tctf, all_cr

work_dir = '../../simulations/production'
compare = sys.argv[1]
field = 'density'
all_sim, all_tctf, all_cr = get_iteration_lists(compare)      
output_list = [30, 40, 50]
field_list = ['density', 'temperature']
for sim in ['isocool']:
    for tctf in all_tctf:
        for cr in all_cr:
            for field in field_list:
                for output in output_list:
                    sim_list, label_list  = pt.get_sim_list(sim, compare, tctf, cr = cr, \
                                                        work_dir = work_dir)

                    plot_multipanel_slices(field, output, sim_list, label_list)



