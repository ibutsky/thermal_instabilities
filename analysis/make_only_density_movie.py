import yt

import numpy as np
import glob
import os
import sys

import matplotlib.pylab as plt
from matplotlib.colors import LogNorm, SymLogNorm
import palettable
import astropy.constants as const

import multiprocessing as mp

import yt_functions as ytf
import plotting_tools as pt

#plt.style.use('dark_background')

half_range = 1
rho0 = 1e-27
T0 = 1e6
mu = 1.22
mh = const.m_p.cgs.value
kb = const.k_B.cgs.value
p0 = (rho0 / mu / mh) * kb*T0

def make_movie_plots(output, folder = '.'):

    ds_loc_list, label_list  = pt.get_sim_list('isocool', compare, tctf, beta = beta,  cr = cr, \
                crdiff = diff, crstream = stream, crheat = heat, work_dir = work_dir, sim_fam = sim_fam)

    fig, ax = plt.subplots(ncols = len(ds_loc_list), nrows = 1, figsize=(1.5*len(ds_loc_list), 3.8), constrained_layout = True)
    cmap = pt.get_cmap(field)

    fig, ax = plt.subplots(ncols = len(ds_loc_list), nrows = 1, figsize=(1.5*len(ds_loc_list), 3.8), constrained_layout = True)
    for i, ds_loc in enumerate(ds_loc_list):
        print(ds_loc)
        if os.path.isfile('%s/DD%04d/DD%04d'%(ds_loc, output, output)):
            
            ds = ytf.load('%s/DD%04d/DD%04d'%(ds_loc, output, output))
            if projection:
                s = yt.ProjectionPlot(ds, 'x', ('gas', field), center = (0, 0, 1), width = (1, 1.8),
                                  weight_field = ('index', 'ones'))
            else:
                s = yt.SlicePlot(ds, 'x', ('gas', field), center = (0, 0, 1), width = (1, 1.8))
            s.set_buff_size(512)
            frb = s.frb
        
            xbins = frb['y'].in_units('kpc')
            ybins = frb['z'].in_units('kpc')
            if field == 'density':
                data_norm = rho0
            elif field == 'temperature':
                data_norm = T0
            elif field == 'cr_pressure':
                data_norm = p0 * cr

            data   = frb[field] / data_norm

            if field == 'density':
                vmin = 1e-1
                vmax = 3
                if projection:
                    vmax = 1.5
            elif field == 'temperature':
                vmin = 5e4 / T0
                vmax = 5e6 / T0
                if projection:
                    vmin = 1e5 / T0
            elif field == 'cr_eta':
                cr_eta = pt.get_cr_eta(ds)
                vmin = cr_eta / 100
                vmax = cr_eta * 100
            elif field == 'cr_pressure':
                vmin = 0.2
                vmax = 2.0
                

            cmap = pt.get_cmap(field)
            pcm = ax[i].pcolormesh(xbins, ybins, data, norm = LogNorm(), cmap = cmap, \
                               vmax = vmax, vmin = vmin, zorder = 1)
        else:
            pcm = ax[i].scatter(0, 0, color = 'white')

        ax[i].set_aspect('equal')
        ax[i].set_xticks([])
        ax[i].set_yticks([])

        ax[i].set_xticklabels([])
        ax[i].set_yticklabels([])
        H_kpc = 43.85 # scale height in kpc                                                                                              


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
    if field == 'cr_pressure':
        cbar = fig.colorbar(pcm,  cax=cbax, orientation = 'vertical', ticks =[0.2, 2])
        cbar.ax.set_yticklabels(['0.2', '2'])
    else:
        cbar = fig.colorbar(pcm,  cax=cbax, orientation = 'vertical')

    if field == 'density':
        label =  '$\\rho / \\rho_0 $'
    elif field == 'temperature':
        label = 'T / T$_0$'
    elif field == 'cr_pressure':
        label = 'P$_c$ / P$_{c,0}$'
    cbar.set_label(label)

    figname = '%s/%04d.png'%(plot_folder, output)
    plt.savefig(figname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.2 )


def make_movie():                    

    output_list = np.arange(0, 101, 1)
#    output_list = np.arange(30, 63) 
    pool = mp.Pool(mp.cpu_count())
    print("Number of processors: ", mp.cpu_count())

    pool.map(make_movie_plots, [output for output in output_list])
    pool.close()

    cwd = os.getcwd()
    os.chdir(plot_folder)
    os.system('ffmpeg -framerate 9 -pattern_type glob -i "*.png" -c:v mpeg4 -pix_fmt yuv420p -q:v 3 %s.mov'%field)
    if projection:
        movie_name_base =  'multipanel_projection_%s_only'%field
    else:
        movie_name_base =  'multipanel_slice_%s_only'%field
    movie_name = pt.get_fig_name(movie_name_base, profile, compare, tctf, beta = beta, use_tctf = 1, \
                           cr=cr, crdiff = diff, crstream = stream, crheat = heat, \
                              sim_fam = '',  loc = '../.')
    print(movie_name[6:-4])
    os.rename('%s.mov'%field, '../%s.mov'%movie_name[6:-4])
    png_files = glob.glob('*.png')
    print(cwd)
   # print(png_files)
    
    for pic in png_files:
        os.remove(pic)
    os.chdir(cwd)


work_dir = '../../simulations'
sim_fam = 'production/high_res'
plot_folder = '../../movies/%s/temp'%sim_fam

projection = False

#field = sys.argv[1]

profile  = 'isocool'
tctf = 1.0
beta = 100
cr = 0

diff = 0
stream = 0
heat = 0

compare = 'beta'
beta = 100
compare = 'transport_multipanel'

cr = 1

for tctf in [0.1, 0.3]:#, 1.0, 3.0]:
    for projection in [False]:
        for cr in [0.01, 0.1, 10]:#for field in ['cr_pressure']:#,'temperature']:
            field = 'cr_pressure'
            make_movie()

