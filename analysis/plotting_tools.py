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

dark_mode = 0

if dark_mode:
    plt.style.use('dark_background')


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

def calculate_cold_fraction(ds, T_min = 3.333333e5, z_min = 0.1, z_max = None, grid_rank = 3):
    ad = ds.all_data()
    
    if (grid_rank == 3):
        z_abs_code = np.abs(ad[('gas', 'z')] / ds.length_unit.in_units('kpc'))

        if z_max == None:
            z_max = ds.domain_right_edge[2].d
        zmask = (z_abs_code >= z_min) & (z_abs_code <= z_max)
        total_mass = np.sum(ad[('gas', 'cell_mass')][zmask].in_units('Msun'))

        cold_mask = (zmask) & (ad[('gas', 'temperature')] <= T_min)
        cold_mass = np.sum(ad[('gas', 'cell_mass')][cold_mask].in_units('Msun'))
    
    elif (grid_rank == 2):
        y_abs_code = np.abs(ad[('gas', 'y')] / ds.length_unit.in_units('kpc'))
        if z_max == None:
            z_max = ds.domain_right_edge[1].d
        ymaxk = (y_abs_code >= z_min) & (y_abs_code <= z_max)
        total_mass = np.sum(ad[('gas', 'cell_mass')][ymask].in_units('Msun'))

        cold_mask = (ymask) & (ad[('gas', 'temperature')] <= T_min)
        cold_mass = np.sum(ad[('gas', 'cell_mass')][cold_mask].in_units('Msun'))

    print(cold_mass, total_mass)
    return cold_mass / total_mass

def calculate_averaged_cold_fraction(output_list, sim_location = '.', T_min = 3.333333e5, z_min =0.1, grid_rank = 3):
    cold_fraction_list = []
    for output in output_list:
        ds_path = '%s/DD%04d/DD%04d'%(sim_location, output, output)
        if os.path.isfile(ds_path):
            ds = yt.load(ds_path)
            cold_fraction_list.append(calculate_cold_fraction(ds, T_min = T_min, z_min = z_min, grid_rank = grid_rank))
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
  #  cube = ds.arbitrary_grid(ds.domain_left_edge, ds.domain_right_edge, ds.domain_dimensions)

    rho = cube[('gas', 'density')].d
    nx, ny, nz = rho.shape
    if field == 'rho':
        fft_field = rho
    elif field == 'drho':
        drho = np.ndarray(shape = rho.shape)
        for i in range(len(rho)):
            if nz == 1:
                rho_slice = rho[:, i, 0]
                rho_ave = np.mean(rho_slice)
                drho[:, i, 0]  = (rho_slice - rho_ave) / rho_ave
            else:
                rho_slice = rho[:, :, i]
                rho_ave = np.mean(rho_slice)
                drho[:, :, i]  = (rho_slice - rho_ave) / rho_ave
        fft_field = drho

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
    L = (ds.domain_right_edge - ds.domain_left_edge).d
    L = np.divide(L, L[-1])
    print(L)
    

    kx = np.fft.rfftfreq(nx)*nx/L[0]
    ky = np.fft.rfftfreq(ny)*ny/L[1]
    kz = np.fft.rfftfreq(nz)*nz/L[2]

    # physical limits to the wavenumbers
    kmin = np.min(1.0/L)
 #   kmax = np.min(0.5*np.array(dims)/L)
    kmax = 0.5*dims[-1]/L[-1]
    print(kmin, kmax)
    if nz == 1:
        kmax = np.max(0.5*np.array(dims)/L)


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


def generate_lists(compare, tctf, crdiff = 0, crstream = 0, crheat=0, cr = 1.0, beta = 100.0):
    if compare == 'tctf':
        tctf_list = [0.1, 0.3, 1.0, 3.0, 10]
        num = len(tctf_list)
        cr_list = num*[cr]
        beta_list = num*[beta]
        diff_list = num*[crdiff]
        stream_list = num*[crstream]
        heat_list = num*[crheat]
    elif compare == 'beta':
        beta_list = [3, 10, 30, 100, 300, 'inf']
        num = len(beta_list)
        tctf_list = num*[tctf]
        cr_list = num*[0]
        diff_list = num*[crdiff]
        stream_list = num*[crstream]
        heat_list = num*[crheat]
    elif compare == 'cr':
        cr_list = [0, 0.01, 0.1, 1.0, 10.0, 100]
        num = len(cr_list)
        tctf_list = num*[tctf]
        beta_list = num*[beta]
        diff_list = num*[crdiff]
        stream_list = num*[crstream]
        heat_list = num*[crheat]
    elif compare == 'diff':
        diff_list = [0, 0, 3, 1]
        num = len(diff_list)
        tctf_list = num*[tctf]
        beta_list = num*[beta]
        cr_list = num*[cr]
        cr_list[0]=0
        stream_list = num*[0]
        heat_list = num*[0]
    elif compare == 'stream':
        beta_list = [100, 100, 100, 10, 10]
        num = len(beta_list)
        cr_list = num*[cr]
        cr_list[0] = 0
        tctf_list = num*[tctf]
        stream_list = num*[1]
        stream_list[0] = 0
        diff_list = num*[0]
        heat_list = num*[0]
        heat_list[-1] = 1
        heat_list[-3] = 1
    elif compare == 'transport':
        tctf_list = [0.1, 0.1, 0.1, 0.1, 0.1]
        cr_list  = [0, 1, 1, 1, 1]
        beta_list = 5*[beta]
        diff_list = [0, 0, 3, 0, 0]
        stream_list = [0, 0, 0, 1, 1]
        heat_list = [0, 0, 0, 0, 1]
    tctf_list = np.array(tctf_list)
    beta_list = np.array(beta_list)
    cr_list = np.array(cr_list)
    diff_list = np.array(diff_list)
    stream_list = np.array(stream_list)
    heat_list = np.array(heat_list)
    return tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list


def generate_sim_list(sim, compare, tctf, crdiff = 0, crstream = 0, crheat=0, \
                      cr = 1.0, beta = 100.0, work_dir = '../../simulations'):
    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list = \
                    generate_lists(compare, tctf, crdiff = crdiff, crstream = crstream,\
                                   crheat=crheat, cr = cr, beta = beta)
    sim_list = []
    for i in range(len(tctf)):
        sim_list.append(get_sim_location(sim, tctf_list[i], beta_list[i], cr_list[i], \
                                         stream_list[i], heat_list[i], work_dir = work_dir))

    return sim_list
                                         

def get_sim_location(sim, tctf, beta, cr, diff = 0, \
            stream = 0, heat = 0, work_dir = '../../simulations'):
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
            sim_location += '_tdiff_%0.1f'%(diff)
        if stream > 0:
            sim_location += '_stream'
        if heat > 0:
            sim_location += '_heat'
    print(sim_location)
    return sim_location

def get_label_name(compare, tctf, beta, cr, crdiff = 0, \
                   crstream = 0, crheat = 0):
    label = '$t_{cool}/t_{ff}$ = %.1f'%tctf
    print(compare)
    if compare == 'cr':
        if cr == 0:
            label = 'No CR'
        else:
            label = 'P$_c$ / P$_g$ = %.1f'%cr
            if cr < 0.1:
                label = 'P$_c$ / P$_g$ = %.2f'%cr
    elif compare == 'beta':
        if beta == 'inf':
            label = 'Hydro'
        else:
            label = '$\\beta = $%.1f'%beta
    elif compare == 'transport':
        if cr > 0:
            label = 'P$_c$ / P$_g$ = %.1f'%cr
        else:
            label = 'No CR'
        if crdiff > 0:
            label = 'diffusion'
        if crstream > 0:
            label = 'streaming'
        if crheat > 0:
            label = 'stream + heat'
    elif compare == 'diff':
        if cr == 0:
            label = 'No CRs'
        elif crdiff == 0:
            label = 't$_{diff}$/t$_{ff}$ = $\infty$'
        else:
            label = 't$_{diff}$/t$_{ff}$ = %.1f'%crdiff
    elif compare == 'stream':
        if cr == 0:
            label = 'No CRs'
        elif crstream:
            label = 'stream'
            if crheat:
                label += ' + heat'
            label += ', $\\beta = %i$'%beta
            
    return label

def get_fig_name(base, sim, compare, tctf, beta, cr=0, crdiff=0, crheat=0, time = -1, use_tctf = 0, loc = '../../plots/'):
    plot_name = '%s/%s_%s'%(loc, base, sim)
    if compare == 'tctf':
        if beta != 'inf':
            plot_name += '_beta_%.1f'%(beta)
        if cr >0:
            plot_name += '_cr_%.1f'%(cr)
        if crdiff > 0:
            plot_name += '_diff_%.1f'%crdiff
    if time < 0 or use_tctf:
        plot_name += '_tctf_%.1f'%tctf
    if compare == 'beta':
        plot_name += '_beta_compare'
    elif compare == 'cr':
        plot_name += '_beta_%.1f_cr_compare'%(beta)
        if crdiff > 0:
            plot_name += '_diff_%.1f'%diff
    elif compare == 'diff':
        plot_name += '_beta_%.1f_cr_%.1f_diff_compare'%(beta, cr)
    elif compare == 'stream':
        plot_name += '_beta_%.1f_cr_%.1f_stream_compare'%(beta, cr)
        if crheat > 0:
            plot_name += '_heat'
    elif compare == 'transport':
        plot_name += '_transport_compare'
    if time > 0:
        plot_name += '_%i'%time

    plot_name += '.png'
    print(plot_name)
    return plot_name


def get_cr_eta(ds):
    base = ds.directory
    if base.__contains__('cr_0.01'):
        cr_eta = 0.01
    elif base.__contains__('cr_0.1'):
        cr_eta = 0.1
    elif base.__contains__('cr_1.0'):
        cr_eta = 1.0
    elif base.__contains__('cr_10.0'):
        cr_eta = 10.0
    else:
        cr_eta = 0
    return cr_eta

def get_title(field, cr_frac = 1):
    if field == 'density':
        title = '$\\rho$'
    elif field == 'pressure':
        title = 'P$_g$'
    elif field == 'temperature':
        title = 'T'
    elif field == 'cr_eta':
        title = 'P$_c$/P$_g$'
    elif field == 'cr_pressure':
        title = 'P$_c$'
    elif field == 'velocity_z':
        title = 'V$_z$'
    elif field == 'magnetic_field_strength':
        title = '|B|'
    return title

def get_zlims(field, cr_frac = 1):
    if field == 'density':
        zlims = (rho0*3e-2, rho0)
    elif field == 'pressure':
        zlims = (p0*1e-2, p0*1e2)
    elif field == 'temperature':
        zlims = (T0 / 20, T0*10)
    elif field == 'cr_eta':
        zlims = (cr_frac*1e-2, cr_frac*1e2)
    elif field == 'cr_pressure':
        zlims = (p0*cr_frac / 50, p0*cr_frac*50)
    elif field == 'velocity_z':
        zlims = (-200, 200)
    elif field == 'magnetic_field_strength':
        zlims = (3e-8, 5e-7)
    return zlims

def get_cmap(field):
    if field =='density':
        cmap = palettable.cmocean.sequential.Tempo_20.mpl_colormap
    elif field == 'pressure':
        cmap = 'magma'
    elif field == 'temperature':
        cmap = palettable.scientific.sequential.LaJolla_20_r.mpl_colormap
    elif field == 'cr_eta':
        cmap = palettable.scientific.sequential.Tokyo_20.mpl_colormap
    elif field == 'cr_pressure':
        cmap = palettable.scientific.sequential.Turku_20.mpl_colormap
    elif field == 'velocity_z':
        cmap = palettable.scientific.diverging.Vik_20.mpl_colormap
    elif field == 'magnetic_field_strength':
        cmap = palettable.scientific.sequential.LaPaz_20.mpl_colormap
    return cmap


def get_masked_data(ds, field, z_min = 0.8, z_max = 1.2, T_cold = 3.33333e5):
    ad = ds.all_data()
    z_abs_code = np.abs(ad[('gas', 'z')] / ds.length_unit.in_units('kpc'))
    temp = np.array(ad[('gas', 'temperature')])
    zmask_cold = (z_abs_code >= z_min) & (z_abs_code <= z_max) & (temp <= T_cold)
    zmask_hot = (z_abs_code >= z_min) & (z_abs_code <= z_max) & (temp > T_cold)

    cold = np.array(ad[('gas', field)][zmask_cold].d)
    hot = np.array(ad[('gas', field)][zmask_hot].d)

    return cold, hot


