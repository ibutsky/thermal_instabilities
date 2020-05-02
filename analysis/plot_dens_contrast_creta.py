import yt
from yt import YTArray
from yt import YTQuantity
import sys
import os
import numpy as np
import matplotlib.pylab as plt
import palettable
import seaborn as sns

import plotting_tools as pt


def plot_density_fluctuation(plot_type = 'mean', compare = 'transport', sim = 'isocool', output = 50,
                              T_cold = 3.3333333e5, zstart = 0.8, zend = 1.2, relative = 0, 
                              work_dir = '../../simulations/', grid_rank = 3, load = True, save = True):

    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list\
                        = pt.generate_lists(compare, 0.1, crdiff = 0, cr = 1, beta = 100, cr_only = 1)
    tctf_list = [0.1, 0.3, 1]
    all_cr_list = [0.01, 0.1, 1, 10]
    print(tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list)
    
    fig, ax = plt.subplots(nrows=1, ncols = 1, figsize = (4.4, 4), sharex = True, sharey = False)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-2, 3e2)
    ax.set_ylim(0.5, 30)
    ax.set_xlabel('$ P_c / P_g $')
    ax.set_ylabel('Average Density Contrast')


#    color_list  = pt.get_color_list('tctf')
    color_list = pt.get_color_list('transport')
    marker = 'o'

    for acr in all_cr_list:
        cr_list = len(beta_list)*[acr]
        for i in range(len(cr_list)):
            for j, tctf in enumerate(tctf_list):
                sim_name = pt.get_sim_name(sim, tctf, beta_list[i], cr_list[i], diff = diff_list[i],
                                               stream = stream_list[i], heat = heat_list[i])
                log_rho, log_T, mass_list = pt.get_2d_hist_data('density', 'temperature', sim_name, 
                                                                zstart = zstart, zend = zend,
                                                          work_dir = work_dir, sim_fam = sim_fam)
                log_rho, log_cre, mass_list = pt.get_2d_hist_data('density', 'cr_eta', sim_name,
                                                                zstart = zstart, zend = zend,
                                                          work_dir = work_dir, sim_fam = sim_fam)
                if len(log_rho) == 0:
                    continue
                rho = log_rho
                T = log_T

                cold_mask = T <= np.log10(T_cold)
                hot_mask  = T > np.log10(T_cold)

                if plot_type == 'mean':
                    rho_cold = np.nanmean(rho[cold_mask])
                    rho_hot = np.nanmean(rho[hot_mask])
                elif plot_type == 'median':
                    rho_cold = np.median(rho[cold_mask])
                    rho_hot = np.mean(rho[hot_mask])

                err_cold = np.nanstd(rho[cold_mask])
                err_hot = np.nanstd(rho[hot_mask])

                err = np.sqrt((err_cold / rho_cold)**2 + (err_hot/rho_hot)**2)*(rho_cold / rho_hot)
                err = np.sqrt(err_cold**2 + err_hot**2)
                
                if cr_list[i] == 10:
                    print(err_cold, len(rho[cold_mask]),  err_hot ,len(rho[hot_mask]))
                creta_time_list, creta_list = pt.get_time_data('cold_creta', sim, tctf, beta_list[i], cr_list[i], \
                                           diff = diff_list[i], stream = stream_list[i], heat = heat_list[i],
                                           T_min = T_cold, zstart = zstart, zend = zend, grid_rank = grid_rank,
                                           load = load, save = save, work_dir = work_dir, sim_fam = sim_fam)
                if len(creta_list) > 0:                    
                    creta = np.nan_to_num(creta_list[output-10:output+10])                    
                    
                if acr == 1 and j == 0:
                    if i == 0:
                        label = 'Advection'
                    else:
                        label = pt.get_label_name('transport', tctf, beta_list[i], cr_list[i], 
                                              crdiff = diff_list[i], crstream = stream_list[i], 
                                              crheat = heat_list[i], counter = i)
                else:
                    label = None
                
                ax.scatter(np.mean(creta), np.power(10, rho_cold - rho_hot), 
                           color = color_list[i+1], label = label, 
                                        marker = marker, alpha = 0.7)

                ax.errorbar(np.mean(creta), np.power(10, rho_cold - rho_hot), xerr = np.std(creta), yerr = np.power(10, err),
                                     color = color_list[i+1], alpha = 0.7)
#                ax.errorbar(np.mean(creta), rho_hot, xerr = np.std(creta), yerr = np.std(rho_hot),
#                                     color = color_list[j])

    ax.legend(fontsize = 8)
    fig.tight_layout()
    fig_basename = 'density_contrast_%s_%s'%(plot_type, sim)
    figname = '../../plots/%s/%s.png'%(sim_fam, fig_basename)
    print(figname)
    plt.savefig(figname, dpi = 300)

sim_fam = 'production'
work_dir = '../../simulations'

for plot_type in ['mean', 'median']:
    plot_density_fluctuation(plot_type, work_dir = work_dir)

