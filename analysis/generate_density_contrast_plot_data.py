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

    cold_mass = np.sum(mass_list[cold_mask])
    hot_mass = np.sum(mass_list[hot_mask])

    # if the cold mass fraction is very low, ignore
    if cold_mass / hot_mass < 0.01:
        rho_cont = -99
        rho_err = 0
    else:
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
    transport_list = np.array([])
    for sim_loc in sim_list:
        sim_name = os.path.basename(sim_loc)
        if sim_has_all_outputs(sim_loc) and not skip_sim(sim_name):
            if sim_name.__contains__('cr'):
                rho_cont, rho_err = get_density_contrast(sim_name, sim_fam = sim_fam, 
                      work_dir = work_dir, T_cold = T_cold, zstart = zstart, zend = zend)
                creta, creta_err = get_creta(sim_name, sim_fam = sim_fam, warm = warm,
                      work_dir = work_dir, T_cold = T_cold, zstart = zstart, zend = zend)
                
                rho_cont_list  = np.append(rho_cont_list,  rho_cont)
                rho_err_list   = np.append(rho_err_list,   rho_err)
                creta_list     = np.append(creta_list,     creta)
                creta_err_list = np.append(creta_err_list, creta_err)
                transport = 1
                if sim_name.__contains__('diff'):
                    transport = 2
                elif sim_name.__contains__('stream'):
                    transport = 3

                transport_list  = np.append(transport_list, transport)

            elif sim_name.__contains__('beta_100.0') and sim_name.__contains__('tctf_0.1'):
                rho_cont, rho_err = get_density_contrast(sim_name, sim_fam = sim_fam,
                      work_dir = work_dir, T_cold = T_cold, zstart = zstart, zend = zend)
                rho_cont_mhd_list = np.append(rho_cont_mhd_list, rho_cont)
                rho_err_mhd_list  = np.append(rho_err_mhd_list, rho_err)

    return  rho_cont_mhd_list, rho_err_mhd_list, rho_cont_list, rho_err_list, creta_list, creta_err_list, transport_list
                                    
def generate_density_contrast_plot_data(plot_type = 'mean', compare = 'transport', profile = 'isocool', output = 50,
                        T_cold = 3.3333333e5, zstart = 0.8, zend = 1.2, relative = 0, warm = False, 
                        sim_fam = 'production',work_dir = '../../simulations/', grid_rank = 3, load = True, save = True):

    rho_cont_mhd, rho_err_mhd, rho_contrast, rho_err, creta, creta_err, transport_list =  \
                           get_plot_data(profile = profile, 
                            sim_fam = sim_fam, work_dir = work_dir, 
                                         zstart = zstart, zend = zend, T_cold = T_cold, warm = warm)
    data = np.vstack((rho_contrast, rho_err, creta, creta_err, transport_list))
    np.save('../../data/%s/dens_contrast_plot_data'%sim_fam, data)
        
    mhd_data = np.vstack((rho_cont_mhd, rho_err_mhd))
    np.save('../../data/%s/dens_contrast_mhd_data'%sim_fam, mhd_data)
        

sim_fam = 'production'
work_dir = '../../simulations'

plot_type = 'mean'
generate_density_contrast_plot_data(plot_type, work_dir = work_dir, profile = 'isocool', warm = False)

