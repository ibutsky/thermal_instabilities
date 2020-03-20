import yt
from yt import YTArray
from yt import YTQuantity
import sys
import os
import numpy as np
import matplotlib.pylab as plt
import palettable

import plotting_tools as pt
import yt_functions as ytf

def plot_density_fluctuation_growth(sim, compare, tctf, beta, cr, diff = 0, stream = 0, heat = 0,
                                    zstart = 0.8, zend = 1.2, 
                                    field = 'density', work_dir = '../../simulations/', grid_rank = 3):


    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list \
                        = pt.generate_lists(compare, tctf, beta = beta, crdiff = crdiff, cr = cr)

    if field == 'cr_pressure':
        mask = cr_list > 0
        tctf_list = tctf_list[mask]
        beta_list = beta_list[mask]
        cr_list = cr_list[mask]
        diff_list = diff_list[mask]
        stream_list = stream_list[mask]
        heat_list = heat_list[mask]
    print(tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list)

    fig, ax = plt.subplots(figsize = (4.4, 4))
    ax.set_yscale('log')
    ax.set_ylim(5e-3, 5)
    ax.set_xlim(0, 10)
    
    gamma = 5./3.
    time_list = np.arange(0, 12, 1)
    wcool = 1.0 / (gamma * 1.0)
    pi = (5./3.) * wcool
    linecolor = 'black'
    if pt.dark_mode:
        linecolor = 'white'
    ax.plot(time_list, 0.02*np.exp(pi*time_list), color = linecolor,\
            linestyle = 'dashed', label = 'Linear Theory', linewidth = 3)

    if compare == 'tctf':
        cpal = palettable.cmocean.sequential.Tempo_5_r.mpl_colors
    else:
        cpal = palettable.scientific.sequential.Batlow_8.mpl_colors


    for i, tctf in enumerate(tctf_list):
        time_list, dzfield_rms_list = pt.get_time_data('rms_fluctuation', sim, tctf, beta_list[i], cr_list[i], \
                                           diff = diff_list[i], stream = stream_list[i], heat = heat_list[i], 
                                           field = field, zstart = zstart, zend = zend, grid_rank = grid_rank, 
                                           load = load, save = save, work_dir = work_dir, sim_fam = sim_fam)

        label = pt.get_label_name(compare, tctf, beta_list[i], cr_list[i], crdiff = diff_list[i], \
                       crstream = stream_list[i], crheat = heat_list[i])
        ax.plot(time_list/tctf, dzfield_rms_list, linewidth = 3, label = label, color = cpal[i])

        if (resolution_compare):
            low_res_sim_location = pt.get_sim_location(sim, tctf, beta_list[i], cr_list[i], \
                                    diff = diff_list[i], stream = stream_list[i],
                                    heat = heat_list[i], work_dir = work_dir+'/low_res')

            low_res_time_list, low_res_dzfield_rms_list = \
                    calculate_rms_fluctuation(low_res_sim_location, output_list, \
                    field = field, grid_rank = grid_rank, data_loc = '../../data/low_res')
            ax.plot(low_res_time_list/tctf, low_res_dzfield_rms_list, linewidth = 3, \
                    linestyle = 'dashed', label = None, color = cpal[i])


            high_res_sim_location = pt.get_sim_location(sim, tctf, beta_list[i], cr_list[i], \
                                    diff = diff_list[i], stream = stream_list[i],
                                    heat = heat_list[i], work_dir = work_dir+'/high_res')

            if not os.path.isdir(high_res_sim_location):
                continue

            high_res_time_list, high_res_dzfield_rms_list = \
                    calculate_rms_fluctuation(high_res_sim_location, output_list, \
                    field = field, grid_rank = grid_rank, data_loc = '../../data/high_res')
            ax.plot(high_res_time_list/tctf, high_res_dzfield_rms_list, linewidth = 3, \
                    linestyle = 'dotted', label = None, color = cpal[i])

     
    ax.set_xlabel('t/t$_{cool}$')
    if field == 'density':
        ax.set_ylabel('$ \\langle \\delta \\rho / \\rho\\rangle_{\\mathrm{rms}}$ ')
    elif field == 'temperature':
        ax.set_ylabel('RMS Temperature Fluctuation')
    elif field == 'cr_pressure':
        ax.set_ylabel('RMS CR Pressure Fluctuation')
    ax.legend()
    fig.tight_layout()
    figname = pt.get_fig_name('%s_fluctuation_growth'%field, sim, compare, \
                              tctf, beta, cr, diff,  sim_fam = sim_fam)
    plt.savefig(figname, dpi = 300)


def make_all_plots(compare, beta = 100, cr = 0, field = 'density'):
    all_tctf = [.1, 0.3, 1, 3]
    all_cr = [0.01, .1, 1, 10]
    for sim in ['isocool']:
        if compare == 'diff' or compare == 'stream' or compare == 'transport':
            for tctf in all_tctf:
                for cr in all_cr:
                    plot_density_fluctuation_growth(sim, compare, tctf, beta, cr, work_dir = work_dir, field = field)

        elif compare == 'cr':
            for tctf in all_tctf:
                plot_density_fluctuation_growth(sim, compare, tctf, beta, cr, work_dir = work_dir, field = field)
                
        elif compare == 'tctf':
            tctf = 0.1
            plot_density_fluctuation_growth(sim, compare, tctf, beta, cr, work_dir = work_dir, field = field)
            
        elif compare == 'beta':
            cr = 0
            for tctf in all_tctf:
                plot_density_fluctuation_growth(sim, compare, tctf, beta, cr, work_dir = work_dir, field = field)


sim_fam = 'production'
work_dir = '../../simulations'
save = True
load = True
resolution_compare = 0

crdiff = 0
compare = sys.argv[1]

#field = 'cr_pressure'
field = 'density'

make_all_plots(compare, field = field)
