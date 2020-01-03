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

sim = sys.argv[1]
half_range = 1
rho0 = 1e-27
T0 = 1e6
mu = 1.22
mh = const.m_p.cgs.value
kb = const.k_B.cgs.value
p0 = (rho0 / mu / mh) * kb*T0

workdir = '../../simulations/'
plot_folder = '../../movies/temp2'


output_list = glob.glob('%s/%s/DD*'%(workdir, sim))

def plot_density_slices(ds, folder = '.'):

    s = yt.SlicePlot(ds, 'x', [('gas', 'density'), ('gas', 'pressure')])
    frb_s = s.frb
    p = yt.ProjectionPlot(ds, 'x', [('gas', 'density'), ('gas', 'velocity_z')], weight_field = 'ones')
    p.set_unit(('gas', 'velocity_z'), 'km/s')
    frb_p = p.frb

    ad = ds.all_data()
    ph = yt.PhasePlot(ad, ('gas', 'pressure'), ('gas', 'entropy'), ('gas', 'cell_mass'), weight_field = None, x_bins = 128, y_bins = 128)
    ph2 = yt.PhasePlot(ad, ('gas', 'z'), ('gas', 'tcool_tff_ratio'), ('gas', 'cell_mass'), weight_field = None, x_bins = 128, y_bins = 128)
    ph2.set_log(('gas', 'z'), False)

    ph.set_unit(('gas', 'cell_mass'), 'Msun')
    ph2.set_unit(('gas', 'cell_mass'), 'Msun')

    cmap_list = [palettable.cmocean.sequential.Tempo_20.mpl_colormap, \
                 palettable.cmocean.diverging.Curl_9.mpl_colormap]

    fig, ax = plt.subplots(ncols = 4, nrows = 2, figsize=(30,14))

    for i, frb in enumerate([frb_s, frb_p]):
        print(i)
        xbins = frb['y'].in_units('kpc')
        ybins = frb['z'].in_units('kpc')
        rho   = frb['density']

        data = rho/rho0

        pcm = ax[i][0].pcolormesh(xbins, ybins, data, norm = LogNorm(), cmap = cmap_list[0], vmin = 3e-2, vmax = 1)
        cbar = fig.colorbar(pcm, ax = ax[i][0], pad=0)
        if i == 0:
            cbar.set_label('Normalized Density (Slice)')
        elif i == 1:
            cbar.set_label('Normalized Density (Projection)')
        ax[i][0].set_xlabel('y (kpc)')
        ax[i][0].set_ylabel('z (kpc)')
    
        
        # calculate density fluctuation
        data = []
        for rho_slice in rho:
            ave_rho = np.mean(rho_slice)
            data.append((rho_slice - ave_rho) / rho_slice)

        pcm = ax[i][1].pcolormesh(xbins, ybins, data, norm=SymLogNorm(0.01), \
                                  cmap = cmap_list[1], vmin = -half_range, vmax = half_range)
        cbar = fig.colorbar(pcm, ax = ax[i][1], pad=0)
        if i == 0:
            cbar.set_label('Density Fluctuation (Slice)')
        elif i == 1:
            cbar.set_label('Density Fluctuation (Projection)')
        ax[i][1].set_xlabel('y (kpc)')
        ax[i][1].set_ylabel('z (kpc)')

    # plot pressure
#    data = []
#    for p_slice in frb_p[('gas',' pressure')]:
#        ave_p = np.mean(p_slice)
#        data.append((p_slice - ave_p) / p_slice)
    data = []
    pres = frb_s[('gas', 'pressure')] / p0
    for p_slice in pres:
        ave_p = np.mean(p_slice)
        data.append((p_slice - ave_p) / p_slice)
    pcm = ax[0][2].pcolormesh(xbins, ybins, data, norm=SymLogNorm(0.01), \
                                  cmap = 'magma',\
                                   vmin = -half_range, vmax = half_range)
    cbar = fig.colorbar(pcm, ax = ax[0][2], pad=0)
    cbar.set_label('Gas Pressure Fluctuation (Slice)')
    ax[0][2].set_xlabel('y (kpc)')
    ax[0][2].set_ylabel('z (kpc)')

    # plot velocity
    pcm = ax[1][2].pcolormesh(xbins, ybins, frb_p[('gas', 'velocity_z')], norm=SymLogNorm(1), \
                                  cmap = palettable.scientific.diverging.Vik_20.mpl_colormap, \
                                   vmin = -100, vmax = 100)
    cbar = fig.colorbar(pcm, ax = ax[1][2], pad=0)
    cbar.set_label('Z-Velocity (Projection)')
    ax[1][2].set_xlabel('y (kpc)')
    ax[1][2].set_ylabel('z (kpc)')
        
    # plot phase plots
    prof = ph.profile
    xbins = prof.x
    ybins = prof.y
    data  = prof[('gas', 'cell_mass')].T
    ax[0][3].set_xscale('log')
    ax[0][3].set_yscale('log')
#    ax[0][3].set_xlim(5e-29, 5e-26)
#    ax[0][3].set_ylim(1e4, 1e7)
    ax[0][3].set_xlim(1e-15, 3e-13)
    ax[0][3].set_ylim(3e-2, 1e3)
    pcm = ax[0][3].pcolormesh(xbins, ybins, data, norm=LogNorm(), \
                                  cmap = palettable.scientific.sequential.Bilbao_16.mpl_colormap,\
                                   vmin = 1e5, vmax = 1e9)
    cbar = fig.colorbar(pcm, ax = ax[0][3], pad=0)
    cbar.set_label('Cell Mass (M$_{\odot}$)')
 #   ax[0][3].set_xlabel('Density (g/cm$^3$)')
 #   ax[0][3].set_ylabel('Temperature (K)')
    ax[0][3].set_xlabel('Gas Pressure ($\\frac{\\mathrm{dyn}}{\\mathrm{cm}^2}$)')
    ax[0][3].set_ylabel('Gas Entropy ($\\mathrm{cm}^2 \cdot \\mathrm{keV}$)')

    
    prof = ph2.profile
    xbins = prof.x.in_units('kpc')
    ybins = prof.y
    data  = prof[('gas', 'cell_mass')].T

    ax[1][3].set_yscale('log')
    ax[1][3].set_xlim(0, max(xbins))
    ax[1][3].set_ylim(1e-3, 1e3)

    pcm = ax[1][3].pcolormesh(xbins, ybins, data, norm=LogNorm(), \
                                  cmap = palettable.cubehelix.jim_special_16_r.mpl_colormap,
                                   vmin = 1e5, vmax = 1e9)
    cbar = fig.colorbar(pcm, ax = ax[1][3], pad=0)
    cbar.set_label('Cell Mass (M$_{\odot}$)')
    ax[1][3].set_xlabel('z (kpc)')
    ax[1][3].set_ylabel('Cooling Time / Free-Fall Time')

    fig.tight_layout()
    return fig, ax

def make_movie_plots(output):
    basename = os.path.basename(output)
    figname = '%s/%s.png'%(plot_folder, basename[2:])
    if not os.path.isfile(figname):
        ds = ytf.load('%s/%s'%(output, basename))
        fig, ax = plot_density_slices(ds)
        plt.savefig(figname, dpi = 300)

                    

pool = mp.Pool(mp.cpu_count())
print("Number of processors: ", mp.cpu_count())

pool.map(make_movie_plots, [output for output in output_list])
pool.close()

cwd = os.getcwd()
os.chdir(plot_folder)
os.system('ffmpeg -r 10 -f image2 -s 1920x1080 -i %04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p density.mov')
os.rename('density.mov', '../%s_density.mov'%(sim))
png_files = glob.glob('*.png')
for pic in png_files:
    os.remove(pic)
os.chdir(cwd)

# ffmpeg -framerate 12 -pattern_type glob -i *.png -c:v mpeg4 -pix_fmt yuv420p -q:v 0 -b 512k movie.mov
