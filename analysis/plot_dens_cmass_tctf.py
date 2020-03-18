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
                              T_cold = 3.3333333e5, zstart = 0.8, zend = 1.2,
                              work_dir = '../../simulations/', grid_rank = 3):

    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list\
                        = pt.generate_lists(compare, tctf, crdiff = diff, cr = cr, beta = beta)
    print(tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list)
    
    fig, ax = plt.subplots(nrows=1, ncols = 2, figsize = (8.8, 4), sharex = True, sharey = False)
    for col in range(2):
        ax[col].set_xscale('log')
        ax[col].set_yscale('log')
        ax[col].set_xlim(.09, 5)
        ax[col].set_xlabel('$t_{cool} / t_{ff}$')

    ax[0].set_ylim(1e-2, 5)
    ax[1].set_ylim(5e-3, 4)
    ax[0].set_ylabel('Density Fluctuation')
    ax[1].set_ylabel('Cold Mass Fraction')
    
    if compare == 'cr':
        cpal = palettable.scientific.sequential.Batlow_6.mpl_colors
        color_list = cpal
    else:
        cpal = palettable.scientific.sequential.Bamako_6.mpl_colors
        color_list = cpal
        if compare == 'stream':
            cpal = palettable.scientific.sequential.Bamako_4.mpl_colors
            color_list = [cpal[0], cpal[0], cpal[0], cpal[1], cpal[1], cpal[1], cpal[2], cpal[2], cpal[2]]
            
    
    for i in range(len(cr_list)):
        out_loc = '../../data' 
        out_name = 'tctf_%s_beta_%.1f_cr_%.2f'%(sim, beta_list[i], cr_list[i])
        if diff_list[i] > 0:
            out_name +='_tdiff_%.1f'%diff_list[i]
        if stream_list[i] > 0:
            out_name +='_stream'
            if heat_list[i] > 0:
                out_name += '_heat'
        out_name += '_%i'%output

        for col, plot_type in enumerate(['density_fluctuation', 'cold_fraction']):
            tctf_list = []
            y_list = []
            data_file = '%s/%s_%s'%(out_loc, plot_type, out_name)
            if os.path.isfile(data_file) and load == True:
                tctf_list, y_list = np.loadtxt(data_file, unpack=True)
            else:
                for tctf in [.1, .3, 1, 3]:
                    sim_location = pt.get_sim_location(sim, tctf, beta_list[i], cr_list[i], \
                                                       diff = diff_list[i], stream = stream_list[i],
                                                       heat = heat_list[i], work_dir = work_dir)

                    sim_base = os.path.basename(sim_location)
                    ds_loc = '%s/DD%04d/DD%04d'%(sim_location, output, output)
                    if os.path.isfile(ds_loc):
                        ds = yt.load(ds_loc)
                        if plot_type == 'density_fluctuation':
                            y_data = pt.calculate_rms_fluctuation(ds)
                        elif plot_type == 'cold_fraction':
                            y_data = pt.calculate_cold_fraction(ds)
                        y_list.append(y_data)
                        tctf_list.append(tctf)
                if save and len(y_list) > 0:
                    outf = open(data_file, 'w')
                    for tctf_out, y_out in zip(tctf_list, y_list):
                        outf.write("%e %e\n"%(tctf_out, y_out))
                    outf.close()

            label = pt.get_label_name(compare, tctf, beta_list[i], cr_list[i], crdiff = diff_list[i], \
                                      crstream = stream_list[i], crheat = heat_list[i])
        
            color = color_list[i]
            linestyle = 'solid'
            marker = 'o'

            if compare != 'cr':
                if cr_list[i] == 0:
                    linestyle = 'dotted'
                    marker = None
                    if i > 0:
                        label = None
                else:
                    if diff_list[i] == 0 and stream_list[i] == 0:
                        linestyle = 'dashed'
                        marker = None
                        if i > 1:
                            label = None
                
            ax[col].plot(tctf_list, y_list, color = color_list[i], label = label, 
                            linewidth = 2, marker = marker, linestyle = linestyle)


    ax[0].legend(fontsize = 8)
    fig.tight_layout()
    fig_basename = 'dens_cfrac_tctf'
    figname = pt.get_fig_name(fig_basename, sim, compare, \
                              tctf, beta, cr, crdiff = diff, crstream = stream, \
                              crheat = heat, time = output, \
                              loc = '../../plots/production')
    print(figname)
    plt.savefig(figname, dpi = 300)


work_dir = '../../simulations/production'
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

#for compare in ['stream', 'diff']:
#    for cr in [.1, 1, 10]:
#        for output in [40]:

for compare in ['stream']:
    for cr in [0.01, 0.1, 1, 10]:
#    for cr in [1.0]:
        for output in [40]:
            plot_density_fluctuation(output, sim, compare, tctf, beta, cr, diff = diff, stream = stream, heat = heat, \
                         work_dir = work_dir)

