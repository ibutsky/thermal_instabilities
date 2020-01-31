import yt

import numpy as np
import glob
import os
import sys

import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import AxesGrid
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


def plot_multipanel_slices(field, output, ds_loc_list, label_list, beta = 100, cr = 0, folder = '.'):

    fig, ax = plt.subplots(ncols = len(ds_loc_list), nrows = 1, figsize=(1.5*len(ds_loc_list), 3.8), constrained_layout = True)    

    for i, ds_loc in enumerate(ds_loc_list):
        print(ds_loc)
        if not os.path.isfile('%s/DD%04d/DD%04d'%(ds_loc, output, output)):
            continue
        ds = ytf.load('%s/DD%04d/DD%04d'%(ds_loc, output, output))

        s = yt.SlicePlot(ds, 'x', ('gas', field), center = (0, 0, 1), width = (1, 2))
        frb = s.frb

        xbins = frb['y'].in_units('kpc')
        ybins = frb['z'].in_units('kpc')
        if field == 'density':
            data_norm = rho0
        else:
            data_norm = 1

        data   = frb[field] / data_norm
        
        if field == 'density':
            vmin = 1e-1
            vmax = 1
            label = '$\\rho / \\rho_0 $'
        elif field == 'temperature':
            vmin = 5e4
            vmax = 5e6
            label = 'T / T$_0$'
        elif field == 'cr_eta':
            cr_eta = pt.get_cr_eta(ds)
            vmin = cr_eta / 100
            vmax = cr_eta * 100
            label = 'P$_c$ / P$_g$'

        cmap = pt.get_cmap(field)
#        ax[i].set_axis_off()
        pcm = ax[i].pcolormesh(xbins, ybins, data, norm = LogNorm(), cmap = cmap, \
                               vmax = vmax, vmin = vmin)
        
        
        ax[i].set_aspect('equal')
        ax[i].set_xticks([])
        ax[i].set_yticks([])
        ax[i].set_title(label_list[i], fontsize = 10)

    #fig.subplots_adjust(left = 0.2, wspace = 0.1)
    fig.tight_layout()
    fig.subplots_adjust(bottom = 0.2)

    pos_l = ax[0].get_position()
    pos_r = ax[-1].get_position()

#    cbax = fig.add_axes([.15, pos.y0, .03, .95*pos.x1*2])
#    cbar = fig.colorbar(pcm,  cax=cbax)
#    cbax.yaxis.set_ticks_position('left')
#    cbax.yaxis.set_label_position('left')
     
    dx = .8
    
    cbax = fig.add_axes([.5-dx/2, 0.17, dx,  0.03])
    cbar = fig.colorbar(pcm,  cax=cbax, orientation = 'horizontal')
#    cbax.yaxis.set_ticks_position('left')
#    cbax.yaxis.set_label_position('left')

    cbar.set_label(label)

    # figname not working right !!!!!!!!!!!!
    figname = pt.get_fig_name('%s_multipanel_slice'%field, sim, compare, tctf, beta = beta, use_tctf = 1, \
                           cr=cr, time = output, loc = '../../plots/production')
#    plt.subplots_adjust(wspace=0.02, top = 0.9)
#    fig.tight_layout()
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
output_list = [40]
field_list = ['density', 'temperature']
#field_list = ['temperature']
all_tctf = [0.1, 0.3, 1, 3]
all_cr = [0]
all_beta = [100]
for sim in ['isocool']:
    for tctf in all_tctf:
        for cr in all_cr:
            for field in field_list:
                for output in output_list:
                    for beta in all_beta:
                        sim_list, label_list  = pt.get_sim_list(sim, compare, tctf, beta = beta,  cr = cr, \
                                                            work_dir = work_dir)
                        print(sim_list)
                        plot_multipanel_slices(field, output, sim_list, label_list, beta = beta, cr = cr)



