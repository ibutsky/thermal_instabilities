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

def plot_density_slices(ds, folder = '.', rho0 = 1e-27, T0 = 1e6, half_range = 1):
    p0 = calculate_p(rho0, T0)

    s = yt.SlicePlot(ds, 'x', [('gas', 'density'), ('gas', 'pressure')])
    frb_s = s.frb
    p = yt.ProjectionPlot(ds, 'x', [('gas', 'density'), ('gas', 'velocity_z')], weight_field = 'ones')
    p.set_unit(('gas', 'velocity_z'), 'km/s')
    frb_p = p.frb

    ad = ds.all_data()
    ph = yt.PhasePlot(ad, ('gas', 'density'), ('gas', 'temperature'), ('gas', 'cell_mass'), weight_field = None)
    ph2 = yt.PhasePlot(ad, ('gas', 'z'), ('gas', 'tcool_tff_ratio'), ('gas', 'cell_mass'), weight_field = None)
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

    data = []
    pres = frb_p[('gas', 'pressure')] / p0
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
    ax[0][3].set_xlim(1e-29, 5e-26)
    ax[0][3].set_ylim(1e4, 1e9)
    pcm = ax[0][3].pcolormesh(xbins, ybins, data, norm=LogNorm(), \
                                  cmap = palettable.scientific.sequential.Bilbao_16.mpl_colormap,\
                                   vmin = 1e5, vmax = 1e9)
    cbar = fig.colorbar(pcm, ax = ax[0][3], pad=0)
    cbar.set_label('Cell Mass (M$_{\odot}$)')
    ax[0][3].set_xlabel('Density (g/cm$^3$)')
    ax[0][3].set_ylabel('Temperature (K)')

    
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


def add_plot(ax, ds, field, field_type = 'slice', view = 'x', cmap = 'viridis', \
                 norm_factor = 1.0, cbar = True, vmin = None, vmax = None, weight_field = 'ones'):
        
    if not field_type.__contains__('phase'):
        if field_type.__contains('slice'):
            s = yt.SlicePlot(ds, view, field)
        elif field_type.__contains('proj'):
            s = yt.ProjectionPlot(ds, view, field, weight_field = weight_field)
        frb = s.frb
        xbins = frb['y'].in_units('kpc')
        ybins = frb['z'].in_units('kpc')
        field_frb = frb[field]
        
        if field_type.__contains__('fluct'):
            data = []
            for field_slice in field_frb:
                ave_field = np.mean(field_slice)
                data.append((field_slice - ave_field) / field_slice)
        else:
            data  = field_frb / norm_factor
                
        if vmin == None:
            vmin = min(data)
        if vmax == None: 
            vmax = max(data)

        pcm = ax.pcolormesh(xbins, ybins, data, norm = LogNorm(), cmap = cmap, vmin = vmin, vmax = vmax)
        if cbar:
            cbar = fig.colorbar(pcm, ax = ax, pad=0)
            cbar.set_label('Normalized Density (Slice)')
        elif i == 1:
            cbar.set_label('Normalized Density (Projection)')
        ax[i][0].set_xlabel('y (kpc)')
        ax[i][0].set_ylabel('z (kpc)')
        

        

def make_movie_plots(output):
    basename = os.path.basename(output)
    figname = '%s/%s.png'%(plot_folder, basename[2:])
    if not os.path.isfile(figname):
        ds = ytf.load('%s/%s'%(output, basename))
        fig, ax = plot_density_slices(ds)
        plt.savefig(figname, dpi = 300)


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


def generate_lists(compare, tctf, crdiff = 0):
    if compare == 'tctf':
        tctf_list = [0.1, 0.3, 1.0, 3.0, 10]
        cr_list = 5*[0]
        beta_list = 5*['inf']
        diff_list = 6*[0]
    elif compare == 'beta':
        tctf_list = 6*[tctf]
        cr_list = 6*[0]
        beta_list = [3, 10, 30, 100, 300, 'inf']
        diff_list = 6*[0]
    elif compare == 'cr':
        tctf_list = 6*[tctf]
        cr_list = [0, 0.01, 0.1, 1.0, 10.0, 100]
        beta_list = 6*[10.0]
        diff_list = 6*[crdiff]

    return tctf_list, beta_list, cr_list, diff_list

def get_sim_location(sim, tctf, beta, cr, diff = 0, work_dir = '../../simulations/'):
    if beta == 'inf':
        sim_location = '%s/%s_tctf_%.1f'%(work_dir, sim, tctf)
    else:
        sim_location =  '%s/%s_tctf_%.1f_beta_%.1f'%(work_dir, sim, tctf, beta)
    if cr > 0:
        if cr < 0.1:
            sim_location += '_cr_%.2f'%(cr)
        else:
            sim_location += '_cr_%0.1f'%(cr)
        if diff > 0:
            sim_location += '_diff_%0.1f'%(diff)
    print(sim_location)
    return sim_location

def get_label_name(compare, tctf, beta, cr):
    label = '$t_{cool}/t_{ff}$ = %.1f'%tctf
    print(compare)
    if compare == 'cr':
        label = 'P$_c$ / P$_g$ = %.1f'%cr
        if cr < 0.1:
            label = 'P$_c$ / P$_g$ = %.2f'%cr
    elif compare == 'beta':
        if beta == 'inf':
            label = 'Hydro'
        else:
            label = '$\\beta = $%.1f'%beta
    return label

def get_fig_name(base, sim, compare, tctf, beta, cr, diff, time = -1, use_tctf = 0, loc = '../../plots'):
    plot_name = '%s/%s_%s'%(loc, base, sim)
    if compare == 'tctf':
        if beta != 'inf':
            plot_name += '_beta_%.1f'%(beta)
    if time < 0 or use_tctf:
        plot_name += '_tctf_%.1f'%tctf
    if compare == 'beta':
        plot_name += '_beta_compare'
    elif compare == 'cr':
        plot_name += '_beta_%.1f_cr_compare'%(beta)
        if diff > 0:
            plot_name += '_diff_%.1f'%diff
    if time > 0:
        plot_name += '_%.1f'%time

    plot_name += '.png'
    print(plot_name)
    return plot_name
