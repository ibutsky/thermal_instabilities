import yt

import numpy as np
import glob
import os

import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import palettable

import multiprocessing as mp


workdir = '/mnt/ceph/users/ibutsky/simulations/test/bomb_debug'
sim = 'low_rg_courant'
plot_folder = '/mnt/ceph/users/ibutsky/plots/movies/%s'%sim
if not os.path.isdir(plot_folder):
    os.mkdir(plot_folder)

output_list = glob.glob('%s/%s/DD*'%(workdir, sim))

def calculate_density_fluctuation(ad, rho0 = 1e-27):    
    x = ad[('gas', 'x')].in_units('kpc')
    y = ad[('gas', 'y')].in_units('kpc')
    z = ad[('gas', 'z')].in_units('kpc')
    rho = ad[('gas', 'density')].in_units('g/cm**3')
    
    drho_over_rho = np.zeros(len(rho))
    rho_norm = rho / rho0
    resolution = int( abs(z[0] - z[-1]) / abs(z[0] - z[1])) + 2
    for i in range(resolution):
        mask = z == z[i]
        ave_rho = np.mean(rho[mask])
        drho_over_rho[mask] = (rho[mask] - ave_rho) / ave_rho

            
    mask = x == x[64]
    xbins = y[mask].reshape(resolution,resolution)
    ybins = z[mask].reshape(resolution,resolution)
    drho = drho_over_rho[mask].reshape(resolution,resolution)
    rho_norm = rho_norm[mask].reshape(resolution, resolution)
    return xbins, ybins, rho_norm, drho

def plot_density_slices(ds, half_range = 0.1, folder = '.', savefig = False):
    ad = ds.all_data()
    xbins, ybins, rho_norm, drho = calculate_density_fluctuation(ad)
    fig, ax = plt.subplots(ncols = 2, nrows = 1, figsize=(16,7))

    cmap = palettable.cmocean.sequential.Tempo_20.mpl_colormap
    pcm = ax[0].pcolormesh(xbins, ybins, rho_norm, norm = LogNorm(),\
                           cmap = cmap, vmin = 1e-2, vmax = 10)
    cbar = fig.colorbar(pcm, ax = ax[0], pad=0)
    cbar.set_label('Normalized Density')
    ax[0].set_xlabel('y (kpc)')
    ax[0].set_ylabel('z (kpc)')
    
    cmap = palettable.cmocean.diverging.Curl_9.mpl_colormap
    pcm = ax[1].pcolormesh(xbins, ybins, drho, \
                           cmap = cmap, vmin = -half_range, vmax = half_range)
    cbar = fig.colorbar(pcm, ax = ax[1], pad=0)
    cbar.set_label('Density Fluctuation')
    ax[1].set_xlabel('y (kpc)')
    ax[1].set_ylabel('z (kpc)')

    fig.tight_layout()
    return fig, ax

def make_movie_plots(output):
    basename = os.path.basename(output)
    figname = '%s/%s.png'%(plot_folder, basename[2:])
    if not os.path.isfile(figname):
        ds = yt.load('%s/%s'%(output, basename))
        fig, ax = plot_density_slices(ds, half_range = 1)
        plt.savefig(figname, dpi = 300)

                    

pool = mp.Pool(mp.cpu_count())
print("Number of processors: ", mp.cpu_count())

pool.map(make_movie_plots, [output for output in output_list])
pool.close()

cwd = os.getcwd()
os.chdir(plot_folder)
os.system('ffmpeg -r 5 -f image2 -s 1920x1080 -i %04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p density.mov')
os.rename('density.mov', '%s_density.mov'%(sim))
os.chdir(cwd)
