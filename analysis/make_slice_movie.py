import yt

import numpy as np
import glob
import os
import sys

import matplotlib.pylab as plt
from matplotlib.colors import LogNorm, SymLogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import palettable
import astropy.constants as const

import multiprocessing as mp

import yt_functions as ytf
import plotting_tools as pt

half_range = 1
rho0 = 1e-27
T0 = 1e6
mu = 1.22
mh = const.m_p.cgs.value
kb = const.k_B.cgs.value
p0 = (rho0 / mu / mh) * kb*T0

sim_family = 'production/constant_crp'
#sim_family = 'skinny'
workdir = '../../simulations/%s'%sim_family
plot_folder = '../../movies/temp2'
movie_loc = '../../movies/%s'%sim_family


def get_zlims(field, cr_frac = 1):
    if field == 'density':
        zlims = (rho0*3e-2, 3*rho0)
    elif field == 'pressure':
        zlims = (p0*1e-2, p0*1e2)
    elif field == 'temperature':
        zlims = (5e4, 5e6)
    elif field == 'cr_eta':
        zlims = (cr_frac*1e-2, cr_frac*1e2)
    elif field == 'cr_pressure':
        zlims = (p0*cr_frac / 50, p0*cr_frac*50)
    elif field == 'velocity_z':
        zlims = (-200, 200)
    elif field == 'magnetic_field_strength':
        zlims = (3e-8, 5e-7)
    return zlims

def plot_density_slices(ds, folder = '.'):
    cr_eta = pt.get_cr_eta(ds)
#    cr_eta = 10
    field_list = [('gas', 'density'), ('gas', 'pressure'), ('gas', 'temperature'), \
                  ('gas', 'velocity_z')]

    if ds.directory.__contains__('beta') or cr_eta > 0:
        field_list.append(('gas', 'magnetic_field_strength'))
        if cr_eta > 0:
            field_list.append(('gas', 'cr_eta'))
            field_list.append(('gas', 'cr_pressure'))

    s = yt.SlicePlot(ds, 'x', field_list)
    frb = s.frb

    ncols = len(field_list)
    fig, ax = plt.subplots(ncols = ncols, nrows = 1, figsize=(2*ncols, 8))

    xbins = frb['y'].in_units('kpc')
    ybins = frb['z'].in_units('kpc')
    rho   = frb['density'] / rho0

    for i, field in enumerate(field_list):
        data = frb[field[1]]
        zlims = get_zlims(field[1], cr_frac = cr_eta)
        if field[1] == 'velocity_z':
            pcm = ax[i].pcolormesh(xbins, ybins, data.in_units('km/s'), norm = SymLogNorm(1), cmap = pt.get_cmap(field[1]), vmin = zlims[0], vmax = zlims[1])
        else:
            pcm = ax[i].pcolormesh(xbins, ybins, data, norm = LogNorm(), cmap = pt.get_cmap(field[1]), vmin = zlims[0], vmax = zlims[1])

        ax[i].set_aspect('equal')
        ax[i].get_xaxis().set_visible(False)
        ax[i].get_yaxis().set_visible(False)
        ax[i].set_title(pt.get_title(field[1]), fontsize = 12)

        cbax = inset_axes(ax[i], width = "90%", height = "3%", loc = 9)
        cbar = fig.colorbar(pcm, cax=cbax, orientation = 'horizontal')

    plt.subplots_adjust(wspace=0.02, top = 0.9)
    return fig, ax

def make_movie_plots(output):
    basename = os.path.basename(output)
    figname = '%s/%s.png'%(plot_folder, basename[2:])
    if not os.path.isfile(figname):
        ds = ytf.load('%s/%s'%(output, basename))
        fig, ax = plot_density_slices(ds)
        plt.savefig(figname, dpi = 300)

                    

sim_list = glob.glob('%s/*tctf*'%workdir)
print(sim_list)
#sim_list = glob.glob('%s/*tdiff_3.0'%workdir)
#sim_list = ['%s/isocool_tctf_0.1_beta_100.0_cr_1.0_tdiff_1.0'%workdir]
#sname = sys.argv[1]
#sim_list = ['%s/%s'%(workdir, sname)]

for sim_loc in sim_list:
    sim_base = os.path.basename(sim_loc)
    movie_name = '%s/slice_movie_%s.mov'%(movie_loc, sim_base)
    make_movie = False
    if not os.path.isfile(movie_name):
        if os.path.isdir('%s/DD0100'%sim_loc):
            make_movie = True
    if make_movie:    
        output_list = glob.glob('%s/DD*'%(sim_loc))
        pool = mp.Pool(mp.cpu_count())
        print("Number of processors: ", mp.cpu_count())

        pool.map(make_movie_plots, [output for output in output_list])
        pool.close()

        cwd = os.getcwd()
        os.chdir(plot_folder)
        os.system('ffmpeg -r 10 -f image2 -s 1920x1080 -i %04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p slice_movie.mov')
        os.rename('slice_movie.mov', '../%s/slice_movie_%s.mov'%(sim_family, sim_base))
        png_files = glob.glob('*.png')
        for pic in png_files:
            os.remove(pic)
        os.chdir(cwd)

# ffmpeg -framerate 12 -pattern_type glob -i *.png -c:v mpeg4 -pix_fmt yuv420p -q:v 0 -b 512k movie.mov
