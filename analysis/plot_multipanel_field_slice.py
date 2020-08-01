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

#import seaborn as sns
pt.sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

half_range = 1
rho0 = 1e-27
T0 = 1e6
mu = 1.22
mh = const.m_p.cgs.value
kb = const.k_B.cgs.value
p0 = (rho0 / mu / mh) * kb*T0

def _inv_T(field, data):
    invT = 1.0 / np.log10(data[('gas', 'temperature')].d)
    return invT

def plot_multipanel_slices(field, output, sim, compare, tctf, beta = 100, cr = 0,\
                           crdiff = 0, crstream = 0, crheat = 0, fixed_time = 0, 
                           weight_field = 'density', projection = False, work_dir = '.'):
    ds_loc_list, label_list  = pt.get_sim_list(sim, compare, tctf, beta = beta,  cr = cr, \
                crdiff = diff, crstream = stream, crheat = heat, work_dir = work_dir, sim_fam = sim_fam)
    print(ds_loc_list)
    fig, ax = plt.subplots(ncols = len(ds_loc_list), nrows = 1, figsize=(1.5*len(ds_loc_list), 3.8), constrained_layout = True)    
    for i, ds_loc in enumerate(ds_loc_list):
        print(ds_loc)
        if fixed_time:
            output_list = [100, 33, 10, 3, 1]
            output = output_list[i]
        if not os.path.isfile('%s/DD%04d/DD%04d'%(ds_loc, output, output)):
            continue
        ds = ytf.load('%s/DD%04d/DD%04d'%(ds_loc, output, output))
        ds.add_field(('gas', 'invT'), function = _inv_T, units = '')
        if projection:
            s = yt.ProjectionPlot(ds, 'x', ('gas', field), center = (0, 0, 1), width = (1, 1.8),
                                  weight_field = weight_field)
        else:
            s = yt.SlicePlot(ds, 'x', ('gas', field), center = (0, 0, 1), width = (1, 1.8))
        s.set_buff_size(1024)
        frb = s.frb

        xbins = frb['y'].in_units('kpc')
        ybins = frb['z'].in_units('kpc')
        if field == 'density':
            data_norm = rho0
        else:
            data_norm = T0

        data   = frb[field] / data_norm
        
        if field == 'density':
            vmin = 1e-1
            vmax = 3   
            if projection:
                vmax = 1.5
            label = '$\\rho / \\rho_0 $'
        elif field == 'temperature':
            vmin = 5e4 / T0
            vmax = 5e6 / T0
            label = 'T / T$_0$'
        elif field == 'cr_eta':
            cr_eta = pt.get_cr_eta(ds)
            vmin = cr_eta / 100
            vmax = cr_eta * 100
            label = 'P$_c$ / P$_g$'

        cmap = pt.get_cmap(field)
        pcm = ax[i].pcolormesh(xbins, ybins, data, norm = LogNorm(), cmap = cmap, \
                               vmax = vmax, vmin = vmin, zorder = 1)
        
        
        ax[i].set_aspect('equal')
#        ax[i].tick_params(direction='in', top=True, right=True, zorder = 10)
        ax[i].set_xticklabels([])
        if i == 0:
            ax[i].set_ylabel('z (kpc)')
        else:
            ax[i].set_yticklabels([])
        H_kpc = 43.85 # scale height in kpc
        ax[i].axhline(0.8*H_kpc, linestyle = 'dashed', color = 'black', linewidth = 0.7)
        ax[i].axhline(1.2*H_kpc, linestyle = 'dashed', color = 'black', linewidth = 0.7)
        if field == 'temperature':
            ax[i].set_xlabel(label_list[i], fontsize = 10)
        else:
            ax[i].set_title(label_list[i], fontsize = 10)

    fig.tight_layout()
    fig.subplots_adjust(bottom = 0.2)
    pos_l = ax[0].get_position().get_points()
    pos_r = ax[-1].get_position().get_points()

    dx = 0.02
    dy = pos_r[1][1] - pos_r[0][1]
    cbar_y = pos_r[0][1]
    cbar_x = pos_r[1][0] + .02
    print(cbar_x, cbar_y)
    cbax = fig.add_axes([cbar_x, cbar_y, dx, dy])
    cbar = fig.colorbar(pcm,  cax=cbax, orientation = 'vertical')
    cbar.set_label(label)

    if projection:
        fig_base = '%s_multipanel_projection'%field
    else:
        fig_base = '%s_multipanel_slice'%field
    figname = pt.get_fig_name(fig_base, sim, compare, tctf, beta = beta, use_tctf = 1, \
                           cr=cr, crdiff = crdiff, crstream = crstream, crheat = crheat, \
                              time = output, sim_fam = sim_fam)
    plt.savefig(figname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)


sim_fam = 'production/high_res'
work_dir = '../../simulations'

compare = 'cr'


sim = 'isocool'
tctf = 1.0
beta = 100
cr = 0

diff = 0
stream = 0
heat = 0

fixed_time = False
projection = False
for output in [50, 60]:
    for field in ['density', 'temperature']:
        for tctf in [0.3, 1]:
            for compare in ['cr', 'tctf']:#, 'transport_multipanel']:
                if compare.__contains__('transport'):
                    cr = 1.0
                if compare.__contains__('tctf'):
                    projection = True
                else:
                    projection = False
                plot_multipanel_slices(field, output, sim, compare, tctf, beta = beta, cr = cr, \
                                       crdiff = diff, crstream = stream, crheat = heat, projection = projection,
                                       work_dir = work_dir, fixed_time = fixed_time, weight_field = ('index', 'ones'))



# tctf = 0.1, .3, 1, 3, 10
# output: 100, 33, 10, 3 1
