import yt
from yt import YTArray
from yt import YTQuantity
import sys
import glob
import os
import numpy as np
import matplotlib.pylab as plt
import palettable
import seaborn as sns

import plotting_tools as pt

def sim_has_all_outputs(sim_loc, min_output = 40, max_output = 60):
    all_outputs = True
    for i in range(min_output, max_output+1):
        if not os.path.isfile('%s/DD%04d/DD%04d'%(sim_loc, i, i)):
            all_outputs = False
    return all_outputs
        
def skip_sim(sim_loc):
    skip = False
    skip_list = ['tctf_1.0', 'tctf_3.0', 'tctf_10.0']
    for s in skip_list:
        if sim_loc.__contains__(s):
            skip = True
    return skip

def get_marker(sim_name):
    if sim_name.__contains__('cr'):
        marker = 'o'
        if sim_name.__contains__('tdiff'):
            marker = 's'
        elif sim_name.__contains__('stream'):
            marker = 'v'
    return marker

def get_sim_cell_size(sim_fam, sim_name):    
    dx = 0.685212  #kpc
    if sim_fam == 'production/high_res':
        dx = 0.342606
    elif sim_fam == 'production/low_res':
        dx = 1.370424
    if sim_name.__contains__('isothermal'):
        dx *= 1.5
    return dx

def get_clump_size(sim_name, sim_fam = 'production', work_dir = '../../simulations',
                  T_cold = 3.33333e5, zstart = 0.8, zend = 1.2):

    sim_location = '%s/%s/%s'%(work_dir, sim_fam, sim_name)
    time_list, clump_data = pt.get_time_data('clump', sim_location = sim_location, use_mpi = True,\
                                           zstart = zstart, zend = zend, load = True, save = True,
                                           work_dir = work_dir, sim_fam = sim_fam)
    nclumps_list, clump_size_list, clump_std_list = clump_data
    nclumps_ave = np.nanmean(nclumps_list[40:61])
    clump_size_ave = np.nanmean(clump_size_list[40:61])
    clump_size_std = np.nanstd(clump_std_list[40:61])

    dx = get_sim_cell_size(sim_fam, sim_name)
    clump_size_ave *= np.power(dx, 1./3.)
    clump_size_std *= np.power(dx, 1./3.)    
    return np.log10(nclumps_ave), np.log10(clump_size_ave), np.log10(clump_size_std)
        
def get_creta(sim_name, sim_fam = 'production', work_dir = '../../simulations',
                        T_cold = 3.33333e5, zstart = 0.8, zend = 1.2):
    log_T, log_creta, mass_list = pt.get_2d_hist_data('temperature', 'cr_eta', sim_name,
                                                        zstart = zstart, zend = zend,
                                                    work_dir = work_dir, sim_fam = sim_fam)
    cold_mask = log_T <= np.log10(T_cold)
    creta_cold_ave = np.nanmean(log_creta[cold_mask])
    err_creta = np.nanstd(log_creta[cold_mask])
    return creta_cold_ave, err_creta

def get_plot_data(profile = 'isocool', sim_fam = 'production', work_dir = '../../simulations', 
                  T_cold = 3.33333e5, zstart = 0.8, zend = 1.2):
    sim_list = sorted(glob.glob('%s/%s/%s*'%(work_dir, sim_fam, profile)))
    nclumps_list    = np.array([])
    clump_size_err_list = np.array([])
    clump_size_list = np.array([])
    creta_list     = np.array([])
    creta_err_list = np.array([])
    marker_list = np.array([])
    for sim_loc in sim_list:
        sim_name = os.path.basename(sim_loc)
        if sim_has_all_outputs(sim_loc) and not skip_sim(sim_name):
            if sim_name.__contains__('cr'):
                nclumps, clump_size, clump_err = get_clump_size(sim_name, sim_fam = sim_fam, 
                      work_dir = work_dir, T_cold = T_cold, zstart = zstart, zend = zend)
                creta, creta_err = get_creta(sim_name, sim_fam = sim_fam,
                      work_dir = work_dir, T_cold = T_cold, zstart = zstart, zend = zend)
                marker = get_marker(sim_name)
                
                nclumps_list  = np.append(nclumps_list,  nclumps)
                clump_size_list = np.append(clump_size_list,   clump_size)
                clump_size_err_list = np.append(clump_size_err_list, clump_err)
                creta_list     = np.append(creta_list,     creta)
                creta_err_list = np.append(creta_err_list, creta_err)
                marker_list    = np.append(marker_list, marker)
                
    return  nclumps_list, clump_size_list, clump_size_err_list, creta_list, creta_err_list, marker_list
                                    
def plot_density_fluctuation(plot_type = 'mean', compare = 'transport', profile = 'isocool', output = 50,
                              T_cold = 3.3333333e5, zstart = 0.8, zend = 1.2, relative = 0, 
                              work_dir = '../../simulations/', grid_rank = 3, load = True, save = True):

    
    fig, ax = plt.subplots(nrows=1, ncols = 2, figsize = (8.8, 4), sharex = True, sharey = False)
    ax[0].set_xlim(-2, 2.5)
    ax[0].set_ylim(0.4, 1.8)
    ax[1].set_ylim(0, 2)
    ax[0].set_xlabel('Log $(P_c / P_g$)')
    ax[0].set_xlabel('Log $(P_c / P_g$)')
    ax[0].set_ylabel('Log Average Clump Size (kpc)')
    ax[1].set_ylabel('Log Average Number of Clumps')

    sim_fam_list = ['production', 'production/high_res', 'production/low_res',
                    'production/Tmin1e4']

    pal = sns.cubehelix_palette(8, start=.5, rot=-.75)

    color_list = [pal[6], pal[4], pal[2], 'orange']
    label_list = ['Fiducial', 'High-res', 'Low-res', 'Tmin = 1e4']
    for i in range(len(sim_fam_list)):
        nclumps_list, clump_size_list, clump_size_err_list, creta_list, creta_err_list, marker_list = \
                           get_plot_data(profile = profile, 
                            sim_fam = sim_fam_list[i], work_dir = work_dir, 
                            zstart = zstart, zend = zend, T_cold = T_cold)
#        ax.axhline(np.mean(rho_cont_mhd), color = color_list[i], linestyle = 'dashed')
        print(sim_fam_list[i], creta_list, clump_size_list)
        for j in range(len(creta_list)):
            ax[0].errorbar(creta_list[j], clump_size_list[j], xerr = creta_err_list[j], 
                           yerr = clump_size_err_list[j], fmt = 'none', color = color_list[i], alpha = 0.3, zorder = 1)
            ax[0].scatter(creta_list[j], clump_size_list[j], color = color_list[i], facecolors = 'none',
                          label = None, marker = marker_list[j], alpha = 1, zorder = 2)
            ax[1].scatter(creta_list[j], nclumps_list[j], color = color_list[i], facecolors = 'none',
                          label = None, marker = marker_list[j], alpha = 1, zorder = 2)

        # just for the label:
        ax[0].scatter(-100, -100, color = color_list[i], label = label_list[i], marker = 'o', alpha = 0.9)
    ax[0].scatter(-100, -100, color = color_list[0], facecolors = 'none', marker = 'o', label = 'Advection')
    ax[0].scatter(-100, -100, color = color_list[0], facecolors = 'none', marker = 's', label = 'Diffusion')
    ax[0].scatter(-100, -100, color = color_list[0], facecolors = 'none', marker = 'v', label = 'Streaming')

    ax[0].legend(fontsize = 8, ncol = 2)
    fig.tight_layout()
    fig_basename = 'clump_size_creta_%s'%(profile)
    figname = '../../plots/%s/%s.png'%(sim_fam, fig_basename)
    print(figname)
    plt.savefig(figname, dpi = 300)

sim_fam = 'production'
work_dir = '../../simulations'

plot_type = 'mean'
plot_density_fluctuation(plot_type, work_dir = work_dir, profile = 'iso')

