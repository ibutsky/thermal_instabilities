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
                              T_cold = 3.3333333e5, zstart = 0.8, zend = 1.2, relative = 0, 
                              work_dir = '../../simulations/', grid_rank = 3):

    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list\
                        = pt.generate_lists(compare, tctf, crdiff = diff, cr = cr, beta = beta)
    mask = cr_list > 0
    tctf_list = tctf_list[mask]
    beta_list = beta_list[mask]
    cr_list = cr_list[mask]
    diff_list = diff_list[mask]
    stream_list = stream_list[mask]
    heat_list = heat_list[mask]
        
    tctf_list = [0.1, 0.3, 1, 3, 10]

    print(tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list)
    
    fig, ax = plt.subplots(nrows=1, ncols = 1, figsize = (4.4, 4), sharex = True, sharey = False)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(.09, 10)
    ax.set_xlabel('$t_{cool} / t_{ff}$')


    ax.set_ylim(1e-2, 5e2)
    ax.set_ylabel('CR Pressure Ratio')

    color_list = pt.get_color_list(compare)    
    for i in range(len(cr_list)):
        for col, plot_type in enumerate(['cold_creta']):
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
                    mean = np.mean(data)
                    err  = np.std(data)
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
                    
            label = pt.get_label_name(compare, tctf, beta_list[i], cr_list[i], crdiff = diff_list[i], \
                                      crstream = stream_list[i], crheat = heat_list[i], counter = i)
        
            color = color_list[i]
            marker = 'o'
            if label is None:
                marker = None

            linestyle = pt.get_linestyle(compare, tctf, beta_list[i], cr_list[i], crdiff = diff_list[i], \
                                              crstream = stream_list[i], crheat = heat_list[i], counter = i)

            x_list = np.array(x_list)
            y_list = np.array(y_list)
            err_list = np.array(err_list)


            if relative == 0:    
                mask = y_list > 0
                x_list   =   x_list[mask]
                y_list   =   y_list[mask]
                err_list = err_list[mask]

                ax.plot(x_list, y_list, color = color_list[i], label = label, 
                            linewidth = 2, marker = marker, linestyle = linestyle)
                ax.errorbar(x_list, y_list, err_list, color = color_list[i])
            elif relative and cr_list[i] > 0:
                ax.plot(x_rel, y_rel, color = color_list[i], label = label,
                            linewidth = 2, marker = marker, linestyle = linestyle)
                ax.errorbar(x_rel, y_rel, err_rel_list, color = color_list[i])
                ax.axhline(y = 1, linestyle = 'dashed', color = 'gray', linewidth = 1)
    ax.legend(fontsize = 8)
    fig.tight_layout()
    fig_basename = 'creta_tctf'
    if relative:
        fig_basename += '_relative'
    figname = pt.get_fig_name(fig_basename, sim, compare, \
                              tctf, beta, cr, crdiff = diff, crstream = stream, \
                              crheat = heat, time = output, sim_fam = sim_fam,\
                              loc = '../../plots')
    print(figname)
    plt.savefig(figname, dpi = 300)

sim_fam = 'production'
work_dir = '../../simulations'
load = True
save = True

def make_all_plots(output, compare, beta = 100, cr = 0.1,\
                   tctf = 0.1, crdiff = 0, crstream = 0, crheat = 0):
    for sim in ['isocool']:
        plot_density_fluctuation(output, sim, compare, tctf, beta, cr, work_dir = work_dir)

                                                                

sim = 'isocool'
tctf = .1

stream = 0
heat = 0
diff = 0
beta = 100

relative = 0
#for compare in ['stream', 'diff']:
#    for cr in [.1, 1, 10]:
#        for output in [40]:

for compare in ['transport']:
    for cr in [0.01, 0.1, 1, 10]:
#    for cr in [1.0]:
        for output in [50]:
            plot_density_fluctuation(output, sim, compare, tctf, beta, cr, diff = diff, stream = stream, heat = heat, \
                                     work_dir = work_dir, relative = relative)

