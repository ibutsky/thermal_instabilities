import yt
from yt import YTQuantity
from yt import YTArray
from yt.data_objects.level_sets.api import *

import numpy as np
import h5py as h5
import glob
import os
import sys

import matplotlib.pylab as plt
from matplotlib.colors import LogNorm, SymLogNorm
import palettable
import seaborn as sns
sns.set_style("ticks",{'axes.grid': True, 'grid.linestyle': '--'})

import astropy.constants as const

import multiprocessing as mp
import yt_functions as ytf

dark_mode = 0

if dark_mode:
    plt.style.use('dark_background')

def _neg_log_T(field, data):
    log_rho = -np.log10(data[('gas', 'temperature')])
    return log_rho

def calculate_pressure(rho, T):
    mu = 1.22
    mh = const.m_p.cgs.value
    kb = const.k_B.cgs.value
    return (rho / mu / mh) * kb*T

def calculate_entropy(rho, T):
    gammam1 = 2./3.
    mh    = YTQuantity(const.m_p.cgs.value, 'g')
    mu    = 1.22
    mu = 1
    kb    = YTQuantity(const.k_B.cgs.value, 'erg / K')    
    T     = YTQuantity(T, 'K')
    rho   = YTQuantity(rho, 'g/cm**3')

    kT = kb*T
    n = rho / (mu * mh)
    e = kT / np.power(n, gammam1)
    return  e.in_units('cm**2 * keV')

def calculate_density(e, p):
    invgamma = 3./5.
    mh    = YTQuantity(const.m_p.cgs.value, 'g')
    mu    = 1.22
    temp = (p / e).in_units('cm**-5')
    n     = np.power(temp, invgamma)

    return n * mu * mh

def calculate_temperature(e, p):
    gamma     = 5./3.
    gammam1   = 2./3.
    invgamma  = 3./5.
    e = e.in_units('erg * cm**2')
    p = p.in_units('dyn / cm**2')

    temp1 = np.power(p, gammam1 / gamma)
    temp2 = np.power(e, invgamma) 
    kb    = YTQuantity(const.k_B.cgs.value, 'erg / K')
    T = temp1*temp2 / kb
#    print (temp1 * temp2)
    return T

def calculate_cold_clump_properties(ds, field = 'neg_log_T'):
    ds.add_field(('gas', 'neg_log_T'), function = _neg_log_T, units = '')
    ad = ds.all_data()
    ad_slice = ds.r[:, :, 0.8:1.2]
    master_clump = Clump(ad, ('gas', field))

    master_clump.add_validator("min_cells", 3)
    master_clump.add_info_item("total_cells")

    c_min = -5
    c_max = ad["gas", field].max()
    step = -0.1
    find_clumps(master_clump, c_min, c_max, step)
    
    leaf_clumps = master_clump.leaves

    radius_list = []
    n_clumps = 0
    H = 43.8535772
    z_min = 0.8 * H
    z_max = 1.2 * H
    for i, clump in enumerate(leaf_clumps):
        z = np.abs(clump[('gas', 'z')].in_units('kpc'))
        z_ave = np.mean(z)
        if (z_ave >= z_min) and (z_ave <= z_max):
            n_clumps += 1
            ncells = len(z)
            radius = np.power(ncells, 1./3.)
            print(i, ncells, radius)
            radius_list.append(radius)

    return n_clumps, np.mean(radius_list), np.std(radius_list)
    
def get_2d_hist_data(xfield, yfield, sim, weighted = True,
                     field = 'density', zstart = 0.8, zend = 1.2, grid_rank = 3,
                     T_min = 3.33333e5, save = True, load = True, data_loc = '../../data',
                     work_dir = '../../simulations', sim_fam = 'production'):

    sim_location = '%s/%s/%s'%(work_dir, sim_fam, sim)
    print(sim_location)
    out_name = '%s/%s/histogram_%s.h5'%(data_loc, sim_fam, sim)
    logx_list = np.array([])
    logy_list = np.array([])
    mass_list = np.array([])
    if os.path.isfile(out_name) and load == True:
        data = h5.File(out_name, 'r')
        if xfield in data.keys() and yfield in data.keys() and 'mass' in data.keys():
            logx_list = data[xfield].value
            logy_list = data[yfield].value
            mass_list = data['mass'].value
        else:
            load = False
        data.close()
    if not os.path.isfile(out_name) or load == False:
        if not os.path.isdir(sim_location):
            return logx_list, logy_list, mass_list
        args_list = []
        output_list = np.arange(40, 61, 1)
        #testing:
#        output_list = [40]
        for output in output_list:
            ds_loc = '%s/DD%04d/DD%04d'%(sim_location, output, output)
            args_list.append((ds_loc, xfield, yfield, weighted, T_min, zstart, zend))
        pool = mp.Pool(mp.cpu_count())
        results = pool.map(calculate_histogram_data_wrapper, args_list)
        pool.close()

        for result in results:
            logx, logy, mass = result
            logx_list = np.append(logx_list, logx)
            logy_list = np.append(logy_list, logy)
            mass_list = np.append(mass_list, mass)
        
        print(logx_list)
        if save and len(logx_list) > 0:
            outf = h5.File(out_name, 'a')
            for key, data in zip([xfield, yfield, 'mass'], [logx_list, logy_list, mass_list]):
                if key not in outf.keys():
                    outf.create_dataset(key, data = data)
            outf.close()

    return logx_list, logy_list, mass_list


def get_time_data(data_type, sim='isocool', tctf=0.1, beta=100, cr=0, diff = 0, stream = 0, heat = 0, use_mpi = True,  
                                  field = 'density', zstart = 0.8, zend = 1.2, grid_rank = 3, 
                                  T_min = 3.33333e5, save = True, load = True, data_loc = '../../data', 
                                  work_dir = '../../simulations', sim_fam = 'production', 
                                   sim_location = None):
    if sim_location is None:
        sim_location = get_sim_location(sim, tctf, beta, cr, diff = diff, stream = stream, heat = heat,
                                    work_dir = work_dir, sim_fam = sim_fam)

    time_list = []
    data_list = []
    if data_type == 'rms_fluctuation' or data_type == 'density_fluctuation':
        out_name = '%s/%s/%s_fluctuation_growth_%s.dat'%(data_loc, sim_fam, field, os.path.basename(sim_location))
    elif data_type == 'cold_fraction':
        out_name = '%s/%s/cold_fraction_growth_%s.dat'%(data_loc, sim_fam, os.path.basename(sim_location))
    elif data_type == 'cold_flux':
        out_name = '%s/%s/cold_mass_flux_growth_%s.dat'%(data_loc, sim_fam, os.path.basename(sim_location))
    elif data_type == 'cold_creta':
        out_name = '%s/%s/cold_creta_growth_%s.dat'%(data_loc, sim_fam, os.path.basename(sim_location))
    elif data_type == 'clump':
        out_name = '%s/%s/clump_growth_%s.dat'%(data_loc, sim_fam, os.path.basename(sim_location))
    else:
        print("ERROR: Data type %s not recognized"%data_type)
    if os.path.isfile(out_name) and load == True:
        if data_type == 'clump':
            time_list, n_clumps, clump_size, clump_std  = np.loadtxt(out_name, skiprows = 1, unpack = True)
            data_list = [n_clumps, clump_size, clump_std]
        else:
            time_list, data_list = np.loadtxt(out_name, skiprows = 1, unpack=True)
            if len(time_list) < 100:
                os.remove(out_name)
                print("WARNING: removing %s"%out_name)
    if not os.path.isfile(out_name) or load == False:
        if not os.path.isdir(sim_location):
            if data_type == 'clump':
                data_list = np.array([[], [], []])
            return time_list, data_list

        args_list = []
        output_list = glob.glob('%s/DD*'%sim_location)
        for output_loc in output_list:
            ds_loc = '%s/%s'%(output_loc, os.path.basename(output_loc))
            args_list.append((ds_loc, data_type, field, T_min, zstart, zend, grid_rank))
        
        if use_mpi:
            pool = mp.Pool(mp.cpu_count())
            results = pool.map(calculate_time_data_wrapper, args_list)
            pool.close()       
#            time_list, data_list = zip(*sorted(results))
            if len(results) > 0:
                time_list, data_list = zip(*sorted(results))
            else:
                time_list = []
                data_list = []
        else:
            time_list = []
            data_list = []
            for args in args_list:
                time, data = calculate_time_data_wrapper(args)
                time_list.append(time)
                data_list.append(data)
            results = zip(time_list, data_list)
            time_list, data_list = zip(*sorted(results))

        if save and len(data_list) > 0:
            outf = open(out_name, 'w')
            for time, data in zip(time_list, data_list):
                if data_type == 'clump':
                    outf.write("%e %e %e %e\n"%(time, data[0], data[1], data[2]))
                else:
                    outf.write("%e %e\n"%(time, data))
            outf.close()
        if data_type == 'clump':
            data_list = list(zip(*data_list))
    return time_list, data_list


def calculate_time_data_wrapper(args):
    output_loc, data_type, field, T_min, zstart, zend, grid_rank = args
    ds = ytf.load(output_loc)
    if data_type == 'rms_fluctuation' or data_type == 'density_fluctuation':
        data = calculate_rms_fluctuation(ds, field = field, zstart = zstart, zend = zend, grid_rank = grid_rank).d
    elif data_type == 'cold_fraction':
        data = calculate_cold_fraction(ds, T_min = T_min, z_min = zstart, z_max = zend, grid_rank = grid_rank).d
    elif data_type == 'cold_flux':
        data = calculate_mass_flux(ds, T_min = T_min, z_min = zstart, z_max = zend, grid_rank = grid_rank)
    elif data_type == 'cold_creta':
        data = calculate_cold_creta(ds, T_min = T_min, z_min = zstart, z_max = zend, grid_rank = grid_rank)
    elif data_type == 'clump':
        data = calculate_cold_clump_properties(ds)
        
    time = ds.current_time.d
    return time, data

def calculate_histogram_data_wrapper(args):
    output_loc, xfield, yfield, weighted, T_min, zstart, zend = args
    ds = ytf.load(output_loc)
    logx, logy, mass = get_log_phase_data(ds, xfield = xfield, yfield = yfield, z_min = zstart, z_max = zend)
    return (logx, logy, mass)


def get_log_phase_data(ds, xfield = 'density', yfield = 'temperature', z_min = 0.8, z_max = 1.2):
    ad = ds.all_data()
    z_abs_code = np.abs(ad[('gas', 'z')] / ds.length_unit.in_units('kpc'))
    zmask = (z_abs_code >= z_min) & (z_abs_code <= z_max)
    xfield = ad[('gas', xfield)][zmask]# / 1e27                                                                       
    yfield = ad[('gas', yfield)][zmask]# / 1e-27                                                                      
    mass = ad[('gas', 'cell_mass')].in_units('Msun')[zmask]
    
    logx = np.log10(xfield)
    logy = np.log10(yfield)

    return logx, logy, mass

def calculate_rms_fluctuation(ds, field = 'density', zstart = .8, zend = 1.2, grid_rank = 3):
    if (grid_rank == 3):
        region1 = ds.r[0, 0, zstart:zend]
        region2 = ds.r[0, 0, -zend:-zstart]
        zlist    = np.append(region1[('gas', 'z')].in_units('kpc'),\
                             region2[('gas', 'z')].in_units('kpc'))

        dzfield = np.array([])
        for z in zlist:
            zslice  = ds.r[:, :, YTQuantity(z, 'kpc')]
            zfield     = zslice[('gas', field)]
            zfield_ave = np.mean(zfield)
            dzfield    = np.append(dzfield, (zfield - zfield_ave) / zfield_ave)
    else:
        all_y = ds.ortho_ray('y', (0, 0))[('gas', 'y')].in_units('kpc')
        ymask = np.abs(all_y / ds.length_unit.in_units('kpc') > zstart) & np.abs(all_y / ds.length_unit.in_units('kpc') < zend)
        ylist = all_y[ymask]
        dzfield = np.array([])
        for y in ylist:
            yray = ds.ortho_ray('x', (y, 0))
            zfield = yray[('gas', field)]
            zfield_ave = np.mean(zfield)
            dzfield    = np.append(dzfield, (zfield - zfield_ave) / zfield_ave)

    dzfield_rms = np.sqrt(np.mean(dzfield**2))
    return dzfield_rms

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

def calculate_mass_flux(ds, z_min = 0.8, z_max = 1.2, grid_rank = 3, T_min = 3.33333e5):
    if (grid_rank == 3):
        ad = ds.all_data()
        z_abs_code = np.abs(ad[('gas', 'z')] / ds.length_unit.in_units('kpc'))

        if z_max == None:
            z_max = ds.domain_right_edge[2].d

        # convert negative velocities to positive in bottom half                           
        vz = ad[('gas', 'velocity_z')]
        vz[ad[('gas', 'z')] < 0] *= -1

        zmask = (z_abs_code >= z_min) & (z_abs_code <= z_max)
        cold_influx_mask = zmask & (vz < 0) & (ad[('gas', 'temperature')] <= T_min)

        rho0 = YTQuantity(1e-27, 'g/cm**3')
        v_cool_flow = ds.length_unit / ds.time_unit # H / tff                                                                               
        cool_flow_flux = (rho0 * v_cool_flow).in_units('Msun / kpc**2 / yr')
#        cool_flow_flux = (ad[('gas', 'density')] * v_cool_flow).in_units('Msun / kpc**2 / yr')
        all_flux = (ad[('gas', 'density')] * ad[('gas', 'velocity_z')]).in_units('Msun / kpc**2 / yr')
        all_flux_norm =(all_flux / cool_flow_flux).d

        cold_influx = all_flux_norm[cold_influx_mask]
        return np.sqrt(np.mean(cold_influx**2))

def calculate_cold_creta(ds, T_min = 3.333333e5, z_min = 0.8, z_max = 1.2, grid_rank = 3):
    ad = ds.all_data()

    if (grid_rank == 3):
        z_abs_code = np.abs(ad[('gas', 'z')] / ds.length_unit.in_units('kpc'))

        if z_max == None:
            z_max = ds.domain_right_edge[2].d
        zmask = (z_abs_code >= z_min) & (z_abs_code <= z_max)
        cold_mask = (zmask) & (ad[('gas', 'temperature')] <= T_min)
    creta = ad[('gas', 'cr_eta')][cold_mask]
    return np.mean(creta)#, np.std(creta)

def plot_density_slices(ds, folder = '.', rho0 = 1e-27, T0 = 1e6, half_range = 1):
    p0 = calculate_pressure(rho0, T0)

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


def generate_lists(compare, tctf, crdiff = 0, crstream = 0, crheat=0, cr = 1.0, beta = 100.0, cr_only = 0):
    if compare == 'tctf':
        tctf_list = [0.1, 0.3, 1.0, 3.0, 10]
        num = len(tctf_list)
        cr_list = num*[cr]
        beta_list = num*[beta]
        diff_list = num*[crdiff]
        stream_list = num*[crstream]
        heat_list = num*[crheat]
    elif compare == 'beta':
        beta_list = [3, 10, 100, 'inf']
        num = len(beta_list)
        tctf_list = num*[tctf]
        cr_list = num*[cr]
        diff_list = num*[crdiff]
        stream_list = num*[crstream]
        heat_list = num*[crheat]
    elif compare == 'cr':
        cr_list = [0, 0.01, 0.1, 1.0, 10.0]
        num = len(cr_list)
        tctf_list = num*[tctf]
        beta_list = num*[beta]
        diff_list = num*[crdiff]
        stream_list = num*[crstream]
        heat_list = num*[crheat]
    elif compare == 'diff':
        diff_list = [0, 0, 10, 3, 1]
        num = len(diff_list)
        tctf_list = num*[tctf]
        beta_list = num*[beta]
        cr_list = num*[cr]
        cr_list[0]=0
        stream_list = num*[0]
        heat_list = num*[0]
    elif compare == 'stream':
        beta_list = [100, 100, 10, 10, 3, 3]
        num = len(beta_list)
        tctf_list = num*[tctf]
        diff_list = num*[0]
        cr_list     = num*[cr]#[cr, cr, cr, cr, cr, cr, cr, cr, cr]
        stream_list = num*[1]#[0, 1, 1, 0, 1, 1, 0, 1, 1]
        heat_list   = [0, 1, 0, 1, 0, 1]

    elif compare == 'transport':
        diff_list   = [0, 0, 10, 3, 1, 0, 0, 0]
        stream_list = [0, 0, 0, 0, 0, 1, 1, 1]
        heat_list   = [0, 0, 0, 0, 0, 1, 1, 1]
        num = len(diff_list)
        tctf_list = num*[tctf]
        cr_list  = num*[cr]
        cr_list[0] = 0
        beta_list = num*[100]
        beta_list[-2] = 10
        beta_list[-1] = 3

    elif compare == 'transport_relative':
        diff_list   = [0, 0, 0, 0,10, 3, 1, 0, 0, 0]
        stream_list = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1]
        heat_list   = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1]
        num = len(diff_list)
        tctf_list = num*[tctf]
        cr_list  = num*[cr]
        cr_list[0] = 0
#        cr_list[1] = 0
#        cr_list[2] = 0
#        cr_list[3] = 0
        beta_list = num*[100]
        beta_list[-2] = 10
        beta_list[-1] = 3
        beta_list[2] = 10
        beta_list[3] = 3

    elif compare == 'transport_pdf':
        diff_list   = [0, 0, 0, 0,10, 3, 1, 0, 0, 0, 0, 0, 0]
        stream_list = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]
        heat_list   = [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0]
        num = len(diff_list)
        tctf_list = num*[tctf]
        cr_list  = num*[cr]
        cr_list[0] = 0
        beta_list = num*[100]
        beta_list[-4] = 3
        beta_list[-3] = 3
        beta_list[-2] = 3
        beta_list[-1] = 3
        beta_list[2] = 10
        beta_list[3] = 3

    elif compare == 'transport_multipanel':
        diff_list   = [0, 0, 3, 0]
        stream_list = [0, 0, 0, 1]
        heat_list   = [0, 0, 0, 1]
        num = len(diff_list)
        tctf_list = num*[tctf]
        cr_list  = num*[cr]
        cr_list[0] = 0
        beta_list = num*[100]
    else:
        print("Unrecognized compare keyword: %s"%compare)
    tctf_list = np.array(tctf_list)
    beta_list = np.array(beta_list)
    cr_list = np.array(cr_list)
    diff_list = np.array(diff_list)
    stream_list = np.array(stream_list)
    heat_list = np.array(heat_list)
    if cr_only:
        mask = cr_list > 0
        tctf_list = tctf_list[mask]
        beta_list = beta_list[mask]
        cr_list   = cr_list[mask]
        diff_list = diff_list[mask]
        stream_list = stream_list[mask]
        heat_list = heat_list[mask]
    return tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list

def generate_sim_list(compare, sim = 'isocool', tctf = 0.3, beta = 100, 
                      cr = 0, diff = 0, stream = 0, heat =  0):
    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list \
        = generate_lists(compare, tctf, beta = beta, cr = cr, 
                         crdiff = diff, crstream = stream, crheat = heat)
    sim_list = []
    for i in range(len(tctf_list)):
        sim_list.append(get_sim_name(sim, tctf_list[i], beta_list[i], cr_list[i], 
                        diff = diff_list[i], stream = stream_list[i], heat = heat_list[i]))
    return sim_list

def generate_label_list(compare, sim = 'isocool', tctf = 0.3, beta = 100,
                      cr = 0, diff = 0, stream = 0, heat =  0):
    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list \
        = generate_lists(compare, tctf, beta = beta, cr = cr,
                         crdiff = diff, crstream = stream, crheat = heat)
    label_list = []
    for i in range(len(tctf_list)):
        label_list.append(get_label_name(compare, tctf_list[i], beta_list[i], cr_list[i],
                        crdiff = diff_list[i], crstream = stream_list[i], crheat = heat_list[i]))
    return label_list

def get_sim_list(sim, compare, tctf=1.0, crdiff = 0, crstream = 0, crheat=0, \
                      cr = 1.0, beta = 100.0, work_dir = '../../simulations', sim_fam = 'production'):
    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list = \
                    generate_lists(compare, tctf, crdiff = crdiff, crstream = crstream,\
                                   crheat=crheat, cr = cr, beta = beta)
    print(tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list)
    sim_list = []
    label_list = []
    
    for i in range(len(tctf_list)):
        sim_list.append(get_sim_location(sim, tctf_list[i], beta_list[i], cr_list[i], \
                                         diff = diff_list[i], stream = stream_list[i], \
                                         heat = heat_list[i], work_dir = work_dir, sim_fam = sim_fam))
        label_list.append(get_label_name(compare, tctf_list[i], beta_list[i], cr_list[i], \
                                         crdiff = diff_list[i], crstream = stream_list[i], \
                                         crheat = heat_list[i]))
    return sim_list, label_list
                                         
def get_sim_name(sim, tctf, beta, cr, diff = 0, stream = 0, heat = 0):
    if beta == 'inf':
        sim_name = '%s_tctf_%.1f'%(sim, tctf)
    else:
        beta = float(beta)
        sim_name =  '%s_tctf_%.1f_beta_%.1f'%(sim, tctf, beta)
    if cr > 0:
        if cr < 0.1:
            sim_name += '_cr_%.2f'%(cr)
        else:
            sim_name += '_cr_%0.1f'%(cr)
        if diff > 0:
            sim_name += '_tdiff_%0.1f'%(diff)
        if stream > 0:
            sim_name += '_stream'
        if heat > 0:
            sim_name += '_heat'
    return sim_name

def get_sim_location(sim, tctf, beta, cr, diff = 0, \
                     stream = 0, heat = 0, work_dir = '../../simulations', sim_fam = 'production'):
    sim_name = get_sim_name(sim, tctf, beta, cr, diff = diff, stream = stream, heat = heat)
    sim_location = '%s/%s/%s'%(work_dir, sim_fam, sim_name)
    print(sim_location)
    return sim_location

def get_label_name(compare, tctf, beta, cr, crdiff = 0, \
                   crstream = 0, crheat = 0, use_crheat = 0, counter = 0):
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
            beta = float(beta)
            label = '$\\beta = $%.1f'%beta
    elif compare.__contains__('transport'):
        if cr == 0:
            label = 'No CR'
            if compare == 'transport_relative':
                label += ', $\\beta = %i$'%beta
        if cr > 0:
            label = 'P$_c$ / P$_g$ = %.1f'%cr
            if compare == 'transport_relative' or compare == 'transport_pdf':
                label += ', $\\beta = %i$'%beta
            if crdiff > 0:
                label = 'Diffusion, $t_{diff} / t_{ff}$ = %i'%crdiff
            elif crstream > 0:
                label = 'Streaming, $\\beta = %i$'%beta
                if crheat == 0:
                    label += ' + no heat'
    elif compare == 'diff':
        if cr == 0:
            label = 'No CRs'
            if counter > 0:
                label = None
        elif crdiff == 0:
            label = 'CR Advection'
            if counter > 1:
                label = None
        else:
            label = 't$_{diff}$/t$_{ff}$ = %.1f'%crdiff
    elif compare == 'stream':
        if cr == 0:
            label = 'No CRs'
            if counter > 0:
                label = None
        elif crstream:
            label = 'stream'
            if use_crheat:
                label += ' + heat'
            label += ', $\\beta = %i$'%beta
        else:
            label = 'CR Advection'
            if counter > 1:
                label = None
            
    return label


def get_linestyle(compare, tctf, beta, cr, crdiff = 0, \
                   crstream = 0, crheat = 0, use_crheat = 0, counter = 0):
    linestyle = 'solid'
    if compare == 'diff' or compare == 'stream':
        if cr == 0:
            linestyle = 'dotted'
        if crheat == 0:
            linestyle = 'dashed'
    elif compare == 'transport':
        if crdiff:
            linestyle = 'dashed'
        elif crstream:
            linestyle = 'dashdot'
    return linestyle


def get_fig_name(base, sim, compare, tctf, beta=100.0, cr=0, crdiff=0, crstream = 0, crheat=0, 
                 time = -1, use_tctf = 0, sim_fam = 'production', loc = '../../plots'):
    plot_name = '%s/%s/%s_%s'%(loc, sim_fam, base, sim)
    if compare == 'tctf':
        if beta != 'inf':
            plot_name += '_beta_%.1f'%(beta)
        if cr >0:
            plot_name += '_cr_%.1f'%(cr)
        if crdiff > 0:
            plot_name += '_diff_%.1f'%crdiff
        plot_name += '_tctf_compare'
    if (time < 0 or use_tctf) and compare != 'tctf':
        plot_name += '_tctf_%.1f'%tctf
    if compare == 'beta':
        plot_name += '_beta_compare'
    elif compare == 'cr':
        plot_name += '_beta_%.1f_cr_compare'%(beta)
        if crdiff > 0:
            plot_name += '_diff_%.1f'%crdiff
        if crstream > 0:
            plot_name += '_stream'
        if crheat > 0:
            plot_name += '_heat'
    elif compare == 'diff':
        plot_name += '_beta_%.1f_cr_%.1f_diff_compare'%(beta, cr)
    elif compare == 'stream':
        plot_name += '_beta_%.1f_cr_%.1f_stream_compare'%(beta, cr)
        if crheat > 0:
            plot_name += '_heat'
    elif compare.__contains__('transport'):
        plot_name += '_cr_%.2f_%s_compare'%(cr, compare)
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

def get_color_list(compare):
    if compare == 'tctf':
        color_list = palettable.cmocean.sequential.Tempo_6_r.mpl_colors
    elif compare == 'beta':
        color_list = palettable.cmocean.sequential.Tempo_5_r.mpl_colors
    elif compare == 'cr':
        color_list = palettable.scientific.sequential.Batlow_6.mpl_colors
    elif compare == 'diff':
        color_list = palettable.scientific.sequential.Bamako_6.mpl_colors
    elif compare == 'stream':
        cpal = palettable.scientific.sequential.Bamako_4.mpl_colors
        color_list = [cpal[0], cpal[0], cpal[0], cpal[1], cpal[1], cpal[1], cpal[2], cpal[2], cpal[2]]
        color_list = [cpal[0], cpal[0], cpal[1], cpal[1], cpal[2], cpal[2]]
    elif compare == 'transport':
        mhd_color = palettable.wesanderson.Darjeeling2_5.mpl_colors[-1]
#        adv_color = palettable.wesanderson.Darjeeling1_4.mpl_colors[1]
#        diff_color = palettable.wesanderson.Chevalier_4.mpl_colors[0]
#        stream_color = palettable.wesanderson.Darjeeling2_5.mpl_colors[1]
        adv_color = palettable.scientific.sequential.Batlow_6.mpl_colors[3]
        ap = palettable.scientific.diverging.Tofino_12.mpl_colors
        adv_color = ap[8]
        #adv_color = 'black'#mhd_color
        tp = palettable.scientific.diverging.Berlin_10.mpl_colors
        diff0_color = tp[6]
        diff1_color = tp[7]
        diff2_color = tp[9]
        
        stream0_color = tp[3]
        stream1_color = tp[2]
        stream2_color = tp[1]
        color_list = [mhd_color, adv_color, diff0_color, diff1_color, diff2_color, stream0_color, stream1_color, stream2_color]
    elif compare == 'transport_relative':
        mhd_color = palettable.wesanderson.Darjeeling2_5.mpl_colors[-1]
        tp = palettable.scientific.diverging.Berlin_12.mpl_colors
        diff0_color = tp[8]
        diff1_color = tp[9]
        diff2_color = tp[10]

        stream0_color = tp[3]
        stream1_color = tp[2]
        stream2_color = tp[1]
        ap = palettable.scientific.diverging.Tofino_12.mpl_colors
        adv1_color = ap[8]
        adv2_color = ap[9]
        adv3_color = ap[10]
        color_list = [mhd_color, adv1_color, adv2_color, adv3_color, diff0_color, 
                      diff1_color, diff2_color, stream0_color, stream1_color, stream2_color]
    return color_list


def constant_T_line(T):
    p_list = []
    entropy_list = []
    rho_list = np.linspace(-30, -24, 100)
    rho_list = np.power(10, rho_list)
    for rho in rho_list:
        p_list.append(calculate_pressure(rho, T))
        entropy_list.append(calculate_entropy(rho, T))
    return np.array(p_list), np.array(entropy_list)

def constant_rho_line(rho):
    rho = YTQuantity(rho, 'g/cm**3')
    p_list = []
    entropy_list = []
    T_list = np.linspace(3, 8, 100)
    T_list = np.power(10, T_list)
    for T in T_list:
        p_list.append(calculate_pressure(rho, T))
        entropy_list.append(calculate_entropy(rho, T))
    return np.array(p_list), np.array(entropy_list)

def constant_entropy_line(e):
    e = YTQuantity(e, 'keV * cm**2')
    rho_list = []
    T_list = []
    p_list = np.linspace(-15.5, -11, 100)
    p_list = np.power(10, p_list)
    for p in p_list:
        p = YTQuantity(p, 'dyn / cm**2')
        rho_list.append(calculate_density(e, p))
        T_list.append(calculate_temperature(e, p))
    return np.array(rho_list), np.array(T_list)

def constant_pressure_line(p):
    p = YTQuantity(p, 'dyn / cm**2')
    rho_list = []
    T_list = []
    e_list = np.linspace(-2, 3, 100)
    e_list = np.power(10, e_list)
    for e in e_list:
        e = YTQuantity(e, 'keV * cm**2')
        rho_list.append(calculate_density(e, p))
        T_list.append(calculate_temperature(e, p))
    return np.array(rho_list), np.array(T_list)
