import yt
from yt import YTQuantity
from yt import YTArray

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



def calculate_p(rho, T):
    mu = 1.22
    mh = const.m_p.cgs.value
    kb = const.k_B.cgs.value
    return (rho / mu / mh) * kb*T

def calculate_drho_rms(ds, zmin = 0.9, zmax = 1.1):
    region1 = ds.r[0, 0, zmin:zmax]
    region2 = ds.r[0, 0, -zmax:-zmin]
    zlist    = np.append(region1[('gas', 'z')].in_units('kpc'),\
                         region2[('gas', 'z')].in_units('kpc'))

    drho = np.array([])
    for z in zlist:
        zslice  = ds.r[:, :, YTQuantity(z, 'kpc')]
        rho     = zslice[('gas', 'density')]
        rho_ave = np.mean(rho)
        drho    = np.append(drho, (rho - rho_ave) / rho_ave)

    drho_rms = np.sqrt(np.mean(drho**2))
    return drho_rms

def calculate_averaged_drho_rms(output_list, sim_location = '.', zmin = 0.9, zmax = 1.1):
    drho_rms_list = []
    for output in output_list:
        ds_path = '%s/DD%04d/DD%04d'%(sim_location, output, output)
        if os.path.isfile(ds_path):
            ds = yt.load(ds_path)
            drho_rms_list.append(calculate_drho_rms(ds, zmin = zmin, zmax = zmax))
        else: 
            print("warning: no such file '%s'"%ds_path)
            drho_rms_list.append(0)
    return np.mean(drho_rms_list)

def calculate_cold_fraction(ds, T_min = 3.333333e5, z_min = 0.1, z_max = 1.2):
    ad = ds.all_data()
    z_min = 0.1
#    z_max = 1.2

    region1 = ds.r[:, :, z_min:]
    region2 = ds.r[:, :, :-z_min]

    total_mass = np.sum(region1[('gas', 'cell_mass')].in_units('Msun')) + \
                 np.sum(region2[('gas', 'cell_mass')].in_units('Msun'))

    cold1 = region1[('gas', 'temperature')] <= T_min
    cold2 = region2[('gas', 'temperature')] <= T_min
    

    cold_mass = np.sum(region1[('gas', 'cell_mass')][cold1].in_units('Msun'))+ \
                np.sum(region2[('gas', 'cell_mass')][cold2].in_units('Msun'))

    return cold_mass / total_mass

def calculate_averaged_cold_fraction(output_list, sim_location = '.', T_min = 3.333333e5, z_min =0.1):
    cold_fraction_list = []
    for output in output_list:
        ds_path = '%s/DD%04d/DD%04d'%(sim_location, output, output)
        if os.path.isfile(ds_path):
            ds = yt.load(ds_path)
            cold_fraction_list.append(calculate_cold_fraction(ds, T_min = T_min, z_min = z_min))
        else:
            print("warning: no such file '%s'"%(ds_path))
            cold_fraction_list.append(0)
    return np.mean(cold_fraction_list)

def fft_comp(ds, field, level = 0 ):
    cube = ds.covering_grid(level, left_edge=ds.domain_left_edge,
                            dims=ds.domain_dimensions)

    rho = cube[('gas', 'density')].d
    if field == 'rho':
        fft_field = rho
    elif field == 'drho':
        drho = np.ndarray(shape = rho.shape)
        for i in range(len(rho)):
            rho_slice = rho[:, :, i]
            rho_ave = np.mean(rho_slice)
            drho[:, :, i]  = (rho_slice - rho_ave) / rho_ave
        fft_field = drho

    nx, ny, nz = rho.shape

    # do the FFTs -- note that since our data is real, there will be 
    # too much information here.  fftn puts the positive freq terms in
    # the first half of the axes -- that's what we keep.  Our      
    # normalization has an '8' to account for this clipping to one 
    # octant.

    ru = np.fft.fftn(fft_field)[0:nx//2+1,0:ny//2+1,0:nz//2+1]
    ru = 8.0*ru/(nx*ny*nz)

    return np.abs(ru)**2


def make_power_spectrum(ds, field = 'drho'):
    dims = ds.domain_dimensions
    nx, ny, nz = dims

    Kk = np.zeros( (nx//2+1, ny//2+1, nz//2+1))
    Kk = fft_comp(ds, field)

    # wavenumbers in units of box length 
    L = np.array([1.0, 1.0, 1.0])

    kx = np.fft.rfftfreq(nx)*nx/L[0]
    ky = np.fft.rfftfreq(ny)*ny/L[1]
    kz = np.fft.rfftfreq(nz)*nz/L[2]

    # physical limits to the wavenumbers
    kmin = np.min(1.0/L)
    kmax = np.min(0.5*np.array(dims)/L)


    kbins = np.arange(kmin, kmax, kmin)
    N = len(kbins)

    # bin the Fourier KE into radial kbins
    kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")
    k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

    whichbin = np.digitize(k.flat, kbins)
    ncount = np.bincount(whichbin)

    P_spectrum = np.zeros(len(ncount)-1)

    for n in range(1,len(ncount)):
        P_spectrum[n-1] = np.sum(Kk.flat[whichbin==n])

    k = kbins[1:N]
    P_spectrum = P_spectrum[1:N]

    return k, P_spectrum
