import yt
import numpy as np
import h5py as h5

import multiprocessing as mp
import glob
import os

import plotting_tools as pt
import yt_functions as ytf


def get_histogram(ds, field, field_range = None, weighted = True, use_log = True, nbins = 100):
    ad = ds.all_data()
    z_abs_code = np.abs(ad[('gas', 'z')] / ds.length_unit.in_units('kpc'))
    z_max = 1.2
    z_min = 0.8
    zmask = (z_abs_code >= z_min) & (z_abs_code <= z_max)
        
    final = ad[('gas', field)][zmask].d
    if weighted:
        weights = np.array(ad[('gas', 'cell_mass')][zmask].in_units('Msun'))
    else:
        weights = np.ones(len(final))
    
    if use_log:
        final = np.log10(final)
        
    if field_range == None:
        if field == 'temperature':
            field_range = (4.25, 7)
        elif field == 'density':
            field_range = (-29, -25)
    
    hist, bin_edges = np.histogram(final, weights = weights, range = field_range, density = True, bins = nbins)
    return hist, bin_edges



def generate_pdf_data(sim):
    time_list = []
    rho = np.ndarray(shape = (100, 100))
    T = np.ndarray(shape = (100, 100))
    for i in range(100):
        time_list.append(float(i) / 10.)
        if os.path.isdir('%s/DD%04d'%(sim, i)):
            ds = ytf.load('%s/DD%04d/DD%04d'%(sim, i, i))
            rho_hist, rho_bin_edges = get_histogram(ds, 'density', weighted = True)
            T_hist, T_bin_edges = get_histogram(ds, 'temperature', weighted = True)

            for j in range(100):
                rho[j][i] = rho_hist[j]
                T[j][i] = T_hist[j]

    temperature_bins = 0.5*(T_bin_edges[1:] + T_bin_edges[:-1])
    density_bins     = 0.5*(rho_bin_edges[1:] + rho_bin_edges[:-1])

    sim_base = os.path.basename(sim)
    pdf_file = h5.File('../../data/production/pdf_data_%s.h5'%sim_base, 'a')
    pdf_file.create_dataset('time', data = np.array(time_list))
    pdf_file.create_dataset('temperature_bins', data = temperature_bins)
    pdf_file.create_dataset('density_bins', data = density_bins)
    pdf_file.create_dataset('temperature_data', data = T)
    pdf_file.create_dataset('density_data', data = rho)
    pdf_file.close()



workdir = '../../simulations/production'
sim_list = glob.glob('%s/*tctf*'%(workdir))

#for sim in sim_list:
#    generate_pdf_data(sim)
pool = mp.Pool(mp.cpu_count())
print("Number of processors: ", mp.cpu_count())
pool.map(generate_pdf_data, [sim for sim in sim_list])
pool.close()
