import yt

import numpy as np
import glob
import os
import sys

import matplotlib.pylab as plt
from matplotlib.colors import LogNorm, SymLogNorm
import palettable

import multiprocessing as mp


sim = sys.argv[1]
half_range = 1
rho0 = 1e-27

workdir = '../../simulations'
plot_folder = '../../movies/%s'%sim

if not os.path.isdir(plot_folder):
    os.mkdir(plot_folder)

output_list = glob.glob('%s/%s/DD*'%(workdir, sim))

def plot_density_slices(ds, folder = '.'):

    s = yt.SlicePlot(ds, 'x', ('gas', 'density'))
    frb_s = s.frb
    p = yt.ProjectionPlot(ds, 'x', ('gas', 'density'), weight_field = 'ones')
    frb_p = p.frb

    cmap_list = [palettable.cmocean.sequential.Tempo_20.mpl_colormap, \
                 palettable.cmocean.diverging.Curl_9.mpl_colormap]

    fig, ax = plt.subplots(ncols = 2, nrows = 2, figsize=(16,14))

    for i, frb in enumerate([frb_s, frb_p]):
        print(i)
        xbins = frb['y'].in_units('kpc')
        ybins = frb['z'].in_units('kpc')
        rho   = frb['density']

        data = rho/rho0

        pcm = ax[i][0].pcolormesh(xbins, ybins, data, norm = LogNorm(), cmap = cmap_list[0], vmin = 1e-1, vmax = 1)
        cbar = fig.colorbar(pcm, ax = ax[i][0], pad=0)
        cbar.set_label('Normalized Density')
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
        cbar.set_label('Density Fluctuation')
        ax[i][1].set_xlabel('y (kpc)')
        ax[i][1].set_ylabel('z (kpc)')

    fig.tight_layout()
    return fig, ax

def make_movie_plots(output):
    basename = os.path.basename(output)
    figname = '%s/%s.png'%(plot_folder, basename[2:])
    if not os.path.isfile(figname):
        ds = yt.load('%s/%s'%(output, basename))
        fig, ax = plot_density_slices(ds)
        plt.savefig(figname, dpi = 300)

                    

pool = mp.Pool(mp.cpu_count())
print("Number of processors: ", mp.cpu_count())

pool.map(make_movie_plots, [output for output in output_list])
pool.close()

cwd = os.getcwd()
os.chdir(plot_folder)
os.system('ffmpeg -r 10 -f image2 -s 1920x1080 -i %04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p density.mov')
os.rename('density.mov', '%s_density.mov'%(sim))
png_files = glob.glob('*.png')
for pic in png_files:
    os.remove(pic)
os.chdir(cwd)

# ffmpeg -framerate 12 -pattern_type glob -i *.png -c:v mpeg4 -pix_fmt yuv420p -q:v 0 -b 512k movie.mov
