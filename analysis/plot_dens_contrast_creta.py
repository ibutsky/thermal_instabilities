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
    sim_list = glob.glob('%s/%s/%s*'%(work_dir, sim_fam, profile))
    rho_cont_mhd_list = np.array([])
    rho_err_mhd_list = np.array([])
    rho_cont_list  = np.array([])
    rho_err_list   = np.array([])
    creta_list     = np.array([])
    creta_err_list = np.array([])
                          
    for sim_loc in sim_list:
        sim_name = os.path.basename(sim_loc)
        if sim_has_all_outputs(sim_loc) and not skip_sim(sim_name):
            if sim_name.__contains__('cr'):
                rho_cont, rho_err = get_density_contrast(sim_name, sim_fam = sim_fam, 
                      work_dir = work_dir, T_cold = T_cold, zstart = zstart, zend = zend)
                creta, creta_err = get_creta(sim_name, sim_fam = sim_fam,
                      work_dir = work_dir, T_cold = T_cold, zstart = zstart, zend = zend)

                rho_cont_list  = np.append(rho_cont_list,  rho_cont)
                rho_err_list   = np.append(rho_err_list,   rho_err)
                creta_list     = np.append(creta_list,     creta)
                creta_err_list = np.append(creta_err_list, creta_err)
            elif sim_name.__contains__('beta_100.0'):
                rho_cont, rho_err = get_density_contrast(sim_name, sim_fam = sim_fam,
                      work_dir = work_dir, T_cold = T_cold, zstart = zstart, zend = zend)
                rho_cont_mhd_list = np.append(rho_cont_mhd_list, rho_cont)
                rho_err_mhd_list  = np.append(rho_err_mhd_list, rho_err)
    return  rho_cont_mhd_list, rho_err_mhd_list, rho_cont_list, rho_err_list, creta_list, creta_err_list
                                    
def plot_density_fluctuation(plot_type = 'mean', compare = 'transport', profile = 'isocool', output = 50,
                              T_cold = 3.3333333e5, zstart = 0.8, zend = 1.2, relative = 0, 
                              work_dir = '../../simulations/', grid_rank = 3, load = True, save = True):

    
    fig, ax = plt.subplots(nrows=1, ncols = 1, figsize = (4.4, 4), sharex = True, sharey = False)
    ax.set_xlim(-2, 2.5)
    ax.set_ylim(-0.1, 1.6)
    ax.set_xlabel('Log $(P_c / P_g$)')
    ax.set_ylabel('Log Average Density Contrast')

#    color_list  = pt.get_color_list('tctf')
    color_list = pt.get_color_list('transport')
    marker = 'o'
    sim_fam_list = ['production', 'production/high_res', 'production/low_res']
    pal = sns.cubehelix_palette(8, start=.5, rot=-.75)
#    color_list = ['black', 'gray']
    color_list = [pal[6], pal[4], pal[2]]
    label_list = ['Fiducial', 'High-res', 'Low-res']
    for i in range(len(sim_fam_list)):
        rho_cont_mhd, rho_err_mhd, rho_contrast, rho_err, creta, creta_err =  \
                           get_plot_data(profile = profile, 
                            sim_fam = sim_fam_list[i], work_dir = work_dir, 
                            zstart = zstart, zend = zend, T_cold = T_cold)
        ax.axhline(np.mean(rho_cont_mhd), color = color_list[i], linestyle = 'dashed')
        
        ax.scatter(creta, rho_contrast, color = color_list[i], label = None, marker = marker, alpha = 0.9)
        ax.errorbar(creta, rho_contrast, xerr = creta_err, yerr = rho_err, fmt = 'none', 
                                     color = color_list[i], alpha = 0.3)

        # just for the label:
        ax.scatter(-100, -100, color = color_list[i], label = label_list[i], marker = marker, alpha = 0.9)

    ax.legend(fontsize = 8)
    fig.tight_layout()
    fig_basename = 'density_contrast_%s'%(profile)
    figname = '../../plots/%s/%s.png'%(sim_fam, fig_basename)
    print(figname)
    plt.savefig(figname, dpi = 300)

sim_fam = 'production'
work_dir = '../../simulations'

plot_type = 'mean'
plot_density_fluctuation(plot_type, work_dir = work_dir, profile = 'iso')

