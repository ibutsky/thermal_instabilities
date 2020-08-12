import yt

import numpy as np
import glob
import os
import sys

import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
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


def plot_multipanel_slices(output, sim_name, sim_loc = '../simulations', sim_fam = 'production', plot_loc = '../../plots',
                           field_list = ['density', 'temperature', 'cr_eta'], weight_field = ('index', 'ones'), projection = False):
    fig, ax = plt.subplots(ncols = len(field_list), nrows = 1, figsize=(1.5*len(field_list), 3.8), constrained_layout = True)    
    ds_loc = '%s/%s/%s'%(sim_loc, sim_fam, sim_name)
    if not os.path.isfile('%s/DD%04d/DD%04d'%(ds_loc, output, output)):
        return
    ds = ytf.load('%s/DD%04d/DD%04d'%(ds_loc, output, output))

    yt_field_list = []
    for field in field_list:
        yt_field_list.append(('gas', field))
    if projection:
        s = yt.ProjectionPlot(ds, 'x', yt_field_list, center = (0, 0, 1), width = (1, 1.8),
                          weight_field = weight_field)
    else:
        s = yt.SlicePlot(ds, 'x', yt_field_list, center = (0, 0, 1), width = (1, 1.8))
    s.set_buff_size(1024)
    frb = s.frb

    xbins = frb['y'].in_units('kpc')
    ybins = frb['z'].in_units('kpc')

    for i, field in enumerate(field_list):
        if field == 'density':
            data_norm = rho0
            vmin = 1e-1
            vmax = 3   
            if projection:
                vmax = 1.5
            label = '$\\rho / \\rho_0 $'
        elif field == 'temperature':
            data_norm = T0
            vmin = 5e4 / T0
            vmax = 5e6 / T0
            label = 'T / T$_0$'
        elif field == 'cr_eta':
            cr_eta = pt.get_cr_eta(ds)
            data_norm = 1
            vmin = cr_eta / 10
            vmax = cr_eta * 10
            label = 'P$_c$ / P$_g$'
        elif field == 'cr_pressure':
            data_norm = p0
            label = 'P$_c$ / P$_{c,0}$'
            vmin = 0.2
            vmax = 2
        data   = frb[field] / data_norm
        cmap = pt.get_cmap(field)
        pcm = ax[i].pcolormesh(xbins, ybins, data, cmap = cmap, norm = LogNorm(),\
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

        cbax = inset_axes(ax[i], width = '90%', height = '4%', loc = 'lower center', bbox_to_anchor = (0.0, -.075, 1, 1),
                          bbox_transform=ax[i].transAxes, borderpad = 0)
        cbar = fig.colorbar(pcm, cax = cbax, orientation = 'horizontal')
        cbar.set_ticks([vmin, vmax])
        cbar.set_ticklabels([vmin, vmax])
        cbar.set_label(label, labelpad = -2)

    fig_base = sim_name
    for field in field_list:
        fig_base += '_%s'%field
        
    if projection:
        fig_base += '_projection'
    else:
        fig_base += '_slice'
        
    figname = '%s/%s.png'%(plot_loc, fig_base)

    plt.savefig(figname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.1)


sim_fam = ''
sim_loc = '/Users/irynabutsky/simulations'
plot_loc = 'plots'

sim_name = 'isocool_tctf_0.3_beta_100.0_cr_1.0'


for output in [40]:
    plot_multipanel_slices(output, sim_name, sim_loc = sim_loc, sim_fam = sim_fam, plot_loc = plot_loc)



# tctf = 0.1, .3, 1, 3, 10
# output: 100, 33, 10, 3 1
