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


def plot_density_fluctuation(output, sim, compare, tctf, beta, cr, diff = 0, stream = 0, heat = 0,
                              T_cold = 3.3333333e5, zstart = 0.8, zend = 1.2, relative = 0, fs = 12, 
                              work_dir = '../../simulations/', grid_rank = 3):

    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list\
                        = pt.generate_lists(compare, tctf, crdiff = diff, cr = cr, beta = beta)
    tctf_list = [0.1, 0.3, 1, 3]#, 10]
    print(tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list)
    
    fig, ax = plt.subplots(nrows=1, ncols = 3, figsize = (4*3, 3.8), sharex = True, sharey = False)
    for col in range(3):
        ax[col].set_xscale('log')
        ax[col].set_yscale('log')
        ax[col].set_xlim(.09, 10)
        ax[col].set_xlabel('$t_{cool} / t_{ff}$', fontsize = fs)

    if relative == 0:
        ax[0].set_ylim(1e-2, 5)
        ax[1].set_ylim(5e-3, 4)
        ax[2].set_ylim(2e-2, 10)
        ax[0].set_ylabel('Density Fluctuation', fontsize = fs)
        ax[1].set_ylabel('Cold Mass Fraction', fontsize = fs)
        ax[2].set_ylabel('Cold Mass Flux', fontsize = fs)
    else:
        ax[0].set_ylim(1e-2, 3)
        ax[1].set_ylim(1e-3, 100)
        ax[2].set_ylim(1e-2, 10)
        ax[0].set_ylabel('Relative Density Fluctuation', fontsize = fs)
        ax[1].set_ylabel('Relative Cold Mass Fraction', fontsize = fs)
        ax[2].set_ylabel('Relative Cold Mass Flux', fontsize = fs)               

    color_list = pt.get_color_list(compare)
    print(beta_list, cr_list, stream_list, heat_list, diff_list)
    for i in range(len(cr_list)):
        for col, plot_type in enumerate(['density_fluctuation', 'cold_fraction', 'cold_flux']):
            for sim_fam in sim_fam_list:
                x_list = []
                y_list = []
                err_list = []

                x_rel = []
                y_rel = []
                err_rel_list = []
                for tctf in tctf_list:
                    time_list, data_list = pt.get_time_data(plot_type, sim, tctf, beta_list[i], cr_list[i], \
                                           diff = diff_list[i], stream = stream_list[i], heat = heat_list[i],
                                           T_min = T_cold, zstart = zstart, zend = zend, grid_rank = grid_rank,
                                           load = load, save = save, work_dir = work_dir, sim_fam = sim_fam)
                    if len(data_list) > 0:
                        data = np.nan_to_num(data_list[output-10:output+10])
                        mean = np.nanmean(data)
                        err  = np.nanstd(data)
                        x_list.append(tctf)
                        y_list.append(mean)
                        err_list.append(err)

                    if relative:
                        time_nocr, data_nocr = pt.get_time_data(plot_type, sim, tctf, beta_list[i], cr = 0, \
                                           diff = 0, stream = 0, heat = 0,
                                           T_min = T_cold, zstart = zstart, zend = zend, grid_rank = grid_rank,
                                           load = load, save = save, work_dir = work_dir, sim_fam = sim_fam)
                        if len(data_list) > 0 and len(data_nocr) > 0:
                            data_nocr = np.nan_to_num(data_nocr[output-10:output+10])
                            mean_nocr = np.mean(data_nocr)
                            err_nocr  = np.std(data_nocr)
                            data =  np.nan_to_num(data_list[output-10:output+10])
                            mean_cr = np.mean(data)
                            err_cr = np.std(data)
                            mean = mean_cr / mean_nocr
                            err = mean * np.sqrt( np.power(err_nocr / mean_nocr, 2) + np.power(err_cr / mean, 2)) 

                            x_rel.append(tctf)
                            y_rel.append(mean)
                            err_rel_list.append(err)
                    
                if sim_fam == 'production':
                    label = pt.get_label_name(compare, tctf, beta_list[i], cr_list[i], crdiff = diff_list[i], \
                                              crstream = stream_list[i], crheat = heat_list[i], counter = i)
                    linestyle = 'solid'
                    if compare == 'stream' and heat_list[i] == 0:
                        linestyle = 'dashed'
                        #label += ', no heating'
 #                   if compare == 'transport' and i == 1:
#                        linestyle = 'dashed'
                else:
                    label =  None
                    linestyle = 'dashed'

                color = color_list[i]
                marker = 'o'
                
                x_list = np.array(x_list)
                y_list = np.array(y_list)
                err_list = np.array(err_list)

                if relative == 0:    
                    mask = y_list > 0
                    x_list   =   x_list[mask]
                    y_list   =   y_list[mask]
                    err_list = err_list[mask]
                    ax[col].plot(x_list, y_list, color = color_list[i], label = label, 
                                 linewidth = 2, marker = marker, linestyle = linestyle)
                    ax[col].errorbar(x_list, y_list, err_list, color = color_list[i], linestyle = '')
                elif relative and cr_list[i] > 0:
                    ax[col].plot(x_rel, y_rel, color = color_list[i], label = label,
                            linewidth = 2, marker = marker, linestyle = linestyle)
                    ax[col].errorbar(x_rel, y_rel, err_rel_list, color = color_list[i], linestyle = '')
                    ax[col].axhline(y = 1, linestyle = 'dashed', color = 'gray', linewidth = 1)
                    ax[col].tick_params(labelsize = fs)
    ax[0].legend(fontsize = 7, ncol = 2, loc = 3)
    fig.tight_layout()
    fig_basename = 'dens_cfrac_cflux_tctf'
    if relative:
        fig_basename += '_relative'
    figname = pt.get_fig_name(fig_basename, sim, compare, \
                              tctf, beta, cr, crdiff = diff, crstream = stream, \
                              crheat = heat, time = output, sim_fam = sim_fam,\
                              loc = '../../plots')
    print(figname)
    plt.savefig(figname, dpi = 300)

sim_fam_list = ['production']#, 'production/high_res']
work_dir = '../../simulations'
load = True
save = True

def make_all_plots(output, compare, beta = 100, cr = 0.1,\
                   tctf = 0.1, crdiff = 0, crstream = 0, crheat = 0):
    for sim in ['isocool']:
        plot_density_fluctuation(output, sim, compare, tctf, beta, cr, work_dir = work_dir)

                                                                

sim = 'isocool'
tctf = .1
output = 50
stream = 0
heat = 0
diff = 0
beta = 100
relative = 0

compare = sys.argv[1]
cr_list = [0.01, 0.1, 1, 10]
if compare == 'transport' or compare == 'cr' or compare == 'stream':
    if compare == 'transport':
#        relative = 1
        compare = 'transport_relative'
        cr_list = [0.01, 0.1, 1, 10]
    else:
        cr_list = [0, 0.01, 0.1, 1, 10]
    for cr in cr_list:
            plot_density_fluctuation(output, sim, compare, tctf, beta, cr, diff = diff, stream = stream, heat = heat, \
                                     work_dir = work_dir, relative = relative)

elif compare == 'beta':
    cr = 0
    relative = 0
    plot_density_fluctuation(output, sim, compare, tctf, beta, cr, diff = diff, stream = stream, heat = heat, \
                                     work_dir = work_dir, relative = relative)
