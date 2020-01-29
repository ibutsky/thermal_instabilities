import multiprocessing as mp
import numpy as np
import h5py as h5
import glob
import os

import matplotlib.pylab as plt
from matplotlib.colors import SymLogNorm, LogNorm

import palettable
import plotting_tools as pt

def plot_pdf_time(data_loc):
    sim_base = os.path.basename(data_loc)
    
    data = h5.File(data_loc, 'r')
    time = data['time']

    for field in ['density', 'temperature']:
        field_bins = data['%s_bins'%field]
        field_data = data['%s_data'%field]

        cmap = pt.get_cmap(field)
        fig, ax = plt.subplots(figsize=(5, 4))
        im = ax.pcolormesh(time, field_bins, field_data, norm = LogNorm(), 
                           cmap = cmap, vmin = 1e-2, vmax = 10)
        cbar = fig.colorbar(im, ax = ax, pad = 0.01)
        cbar.set_label('Frequency')
        ax.set_xlabel('t / t$_{cool}$')
        if field == 'temperature':
            ax.set_ylabel('Log Temperature (K)')
        elif field == 'density':
            ax.set_ylabel('Log Density (g cm$^{-3}$)')

        fig.tight_layout()
        plt.savefig('%s/pdf_time_%s_%s.png'%(plot_loc, field, sim_base), dpi = 300)
        plt.clf()

workdir = '../../data/production'
data_loc_list = glob.glob('%s/pdf*.h5'%(workdir))
plot_loc = '../../plots/production'

pool = mp.Pool(mp.cpu_count())
print("Number of processors: ", mp.cpu_count())
pool.map(plot_pdf_time, [data_loc for data_loc in data_loc_list])
pool.close()


