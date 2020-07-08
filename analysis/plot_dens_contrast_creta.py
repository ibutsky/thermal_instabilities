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
    skip_list = ['tctf_3.0', 'tctf_10.0']
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

def get_density_contrast(sim_name, sim_fam = 'production', work_dir = '../../simulations',
                  T_cold = 3.33333e5, zstart = 0.8, zend = 1.2):

    log_rho, log_T, mass_list = pt.get_2d_hist_data('density', 'temperature', sim_name,
                                                            zstart = zstart, zend = zend,
                                                            work_dir = work_dir, sim_fam = sim_fam)
    if len(log_rho) == 0:
        print("%s has len(rho) = 0"%sim_loc)

    cold_mask = log_T <= np.log10(T_cold)
    hot_mask  = log_T > np.log10(T_cold)

    rho_cold_ave = np.nanmean(log_rho[cold_mask])
    rho_hot_ave = np.nanmean(log_rho[hot_mask])
    err_cold = np.nanstd(log_rho[cold_mask])
    err_hot = np.nanstd(log_rho[hot_mask])

    err = np.sqrt(err_cold**2 + err_hot**2)            
    # in log space, rho_c / rho_h
    rho_cont = rho_cold_ave - rho_hot_ave
    err = np.sqrt(err_cold**2 + err_hot**2)
    
    return rho_cont, err
        
def get_creta(sim_name, sim_fam = 'production', work_dir = '../../simulations',
                        T_cold = 3.33333e5, zstart = 0.8, zend = 1.2, warm = False):
    log_T, log_creta, mass_list = pt.get_2d_hist_data('temperature', 'cr_eta', sim_name,
                                                        zstart = zstart, zend = zend,
                                                    work_dir = work_dir, sim_fam = sim_fam)
    temp_mask = log_T <= np.log10(T_cold)
    if warm:
        temp_mask = log_T > np.log10(T_cold)
        
    creta_ave = np.nanmean(log_creta[temp_mask])
    err_creta = np.nanstd(log_creta[temp_mask])
    return creta_ave, err_creta

def get_plot_data(profile = 'isocool', sim_fam = 'production', work_dir = '../../simulations', 
                  T_cold = 3.33333e5, zstart = 0.8, zend = 1.2, warm = False):
    sim_list = glob.glob('%s/%s/%s*'%(work_dir, sim_fam, profile))
    rho_cont_mhd_list = np.array([])
    rho_err_mhd_list = np.array([])
    rho_cont_list  = np.array([])
    rho_err_list   = np.array([])
    creta_list     = np.array([])
    creta_err_list = np.array([])
    marker_list = np.array([])
    for sim_loc in sim_list:
        sim_name = os.path.basename(sim_loc)
        if sim_has_all_outputs(sim_loc) and not skip_sim(sim_name):
            if sim_name.__contains__('cr'):
                rho_cont, rho_err = get_density_contrast(sim_name, sim_fam = sim_fam, 
                      work_dir = work_dir, T_cold = T_cold, zstart = zstart, zend = zend)
                creta, creta_err = get_creta(sim_name, sim_fam = sim_fam, warm = warm,
                      work_dir = work_dir, T_cold = T_cold, zstart = zstart, zend = zend)
                marker = get_marker(sim_name)
                
                rho_cont_list  = np.append(rho_cont_list,  rho_cont)
                rho_err_list   = np.append(rho_err_list,   rho_err)
                creta_list     = np.append(creta_list,     creta)
                creta_err_list = np.append(creta_err_list, creta_err)
                marker_list    = np.append(marker_list, marker)

            elif sim_name.__contains__('beta_100.0') and sim_name.__contains__('tctf_0.1'):
                rho_cont, rho_err = get_density_contrast(sim_name, sim_fam = sim_fam,
                      work_dir = work_dir, T_cold = T_cold, zstart = zstart, zend = zend)
                rho_cont_mhd_list = np.append(rho_cont_mhd_list, rho_cont)
                rho_err_mhd_list  = np.append(rho_err_mhd_list, rho_err)

    return  rho_cont_mhd_list, rho_err_mhd_list, rho_cont_list, rho_err_list, creta_list, creta_err_list, marker_list
                                    
def plot_density_fluctuation(plot_type = 'mean', compare = 'transport', profile = 'isocool', output = 50,
                              T_cold = 3.3333333e5, zstart = 0.8, zend = 1.2, relative = 0, warm = False,
                              work_dir = '../../simulations/', grid_rank = 3, load = True, save = True):

    
    fig, ax = plt.subplots(nrows=1, ncols = 1, figsize = (4.4, 4), sharex = True, sharey = False)
    ax.set_xlim(-2, 2.5)
    ax.set_ylim(-0.1, 1.6)
    if warm:
        ax.set_xlabel('Log $(P_c / P_g)_{\mathrm{hot}}$')
    else:
        ax.set_xlabel('Log $(P_c / P_g)_{\mathrm{cold}}$')
    ax.set_ylabel('Log Average Density Contrast')

#    color_list  = pt.get_color_list('tctf')
    color_list = pt.get_color_list('transport')
    marker = 'o'
    sim_fam_list = ['production', 'production/high_res', 'production/low_res',
                    'production/Tmin1e4']
    pal = sns.cubehelix_palette(8, start=.5, rot=-.75)
#    color_list = ['black', 'gray']
    color_list = [pal[6], pal[4], pal[2], 'orange']
    label_list = ['Fiducial', 'High-res', 'Low-res', 'Tmin = 1e4']
    for i in range(len(sim_fam_list)):
        rho_cont_mhd, rho_err_mhd, rho_contrast, rho_err, creta, creta_err, marker_list =  \
                           get_plot_data(profile = profile, 
                            sim_fam = sim_fam_list[i], work_dir = work_dir, 
                                         zstart = zstart, zend = zend, T_cold = T_cold, warm = warm)
        ax.axhline(np.mean(rho_cont_mhd), color = color_list[i], linestyle = 'dashed', zorder = 0)
        
        for j in range(len(creta)):
            ax.errorbar(creta[j], rho_contrast[j], xerr = creta_err[j], yerr = rho_err[j], fmt = 'none', 
                                     color = color_list[i], alpha = 0.3, zorder = 1)
            ax.scatter(creta[j], rho_contrast[j], color = color_list[i], facecolors = 'none',
                       label = None, marker = marker_list[j], alpha = 1, zorder = 2)

        # just for the label:
        ax.scatter(-100, -100, color = color_list[i], label = label_list[i], marker = marker, alpha = 0.9)

    #just for the label
    ax.scatter(-100, -100, color = color_list[0], facecolors = 'none', marker = 'o', label = 'Advection')
    ax.scatter(-100, -100, color = color_list[0], facecolors = 'none', marker = 's', label = 'Diffusion')
    ax.scatter(-100, -100, color = color_list[0], facecolors = 'none', marker = 'v', label = 'Streaming')

    # analytic line
    x = np.linspace(0, 3, 20)
#    y = np.power((1 + 1/x), 3./2.)
    y = -3./2. *  x + 1.5
    if warm:
        y = 3./2 * np.log10(1 + 1/np.power(10, x))
    ax.plot(x, y, color = 'black', linestyle = 'dotted') 

    ax.legend(fontsize = 8, ncol = 2, loc = 3)
    fig.tight_layout()
    fig_basename = 'density_contrast_%s'%(profile)
    if warm:
        fig_basename += '_warm'
    figname = '../../plots/%s/%s.png'%(sim_fam, fig_basename)
    print(figname)
    plt.savefig(figname, dpi = 300)

sim_fam = 'production'
work_dir = '../../simulations'

plot_type = 'mean'
plot_density_fluctuation(plot_type, work_dir = work_dir, profile = 'iso', warm = False)

