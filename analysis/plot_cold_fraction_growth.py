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

def plot_cold_fraction_growth(sim, compare, tctf, beta, cr, diff = 0, stream = 0, heat = 0,
                              T_min = 3.3333333e5, zstart = 0.8, zend = 1.2,
                              work_dir = '../../simulations/', grid_rank = 3):

    tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list\
                        = pt.generate_lists(compare, tctf, crdiff = crdiff, cr = cr)

    fig, ax = plt.subplots(figsize = (4, 4))
    ax.set_yscale('log')
    ax.set_ylim(5e-3, 5)
    ax.set_xlim(0, 10)
    
    cpal = palettable.scientific.sequential.Batlow_6.mpl_colors

    output_list = np.linspace(0, 100, 10)
    for i, tctf in enumerate(tctf_list):
        sim_location = pt.get_sim_location(sim, tctf, beta_list[i], cr_list[i], \
                                           diff = diff_list[i], stream = stream_list[i], 
                                           heat = heat_list[i], work_dir = work_dir)

        out_name = '../../data/cold_fraction_growth_%s.dat'%(os.path.basename(sim_location))
        if os.path.isfile(out_name) and load == True:
            time_list, cold_fraction_list = np.loadtxt(out_name, unpack=True)

        else:
            if not os.path.isdir(sim_location):
                continue

            time_list = []
            cold_fraction_list = []
            for output in output_list:
                ds_loc = '%s/DD%04d/DD%04d'%(sim_location, output, output)
                if os.path.isfile(ds_loc):
                    ds = yt.load(ds_loc)
                    cold_frac = pt.calculate_cold_fraction(ds, T_min = T_min, z_min = zstart, z_max = zend, grid_rank = grid_rank)
                    time_list.append(ds.current_time/tctf)
                    cold_fraction_list.append(cold_frac)

        label = pt.get_label_name(compare, tctf, beta_list[i], cr_list[i], crdiff = diff_list[i], \
                           crstream = stream_list[i], crheat = heat_list[i])
        ax.plot(time_list, cold_fraction_list, linewidth = 2, label = label, color = cpal[i])

        if save and len(cold_fraction_list) > 0:
            outf = open(out_name, 'w')
            for i in range(len(time_list)):
                outf.write("%e %e\n"%(time_list[i], cold_fraction_list[i]))
            outf.close()
     
    ax.set_xlabel('t/t$_{cool}$')
    ax.set_ylabel('Cold Mass Fraction')

    ax.legend()
    fig.tight_layout()
    figname = pt.get_fig_name('cold_fraction_growth', sim, compare, tctf, beta, cr, diff
                              loc = '../../plots/production')
    plt.savefig(figname, dpi = 300)

def make_all_plots(compare, beta = 100, cr = 0.1):
    all_tctf = [.1, 0.3, 1, 3]
    all_cr = [0.01, .1, 1, 10]
    for sim in ['isothermal', 'isocool']:
        if compare == 'diff' or compare == 'stream':
            for tctf in all_tctf:
                for cr in all_cr:
                    plot_cold_fraction_growth(sim, compare, tctf, beta, cr, work_dir = work_dir)

        elif compare == 'cr':
            for tctf in all_tctf:
                    plot_cold_fraction_growth(sim, compare, tctf, beta, cr, work_dir = work_dir)

        

work_dir = '../../simulations/production'
load = True
save = True

crdiff = 0

compare = sys.argv[1]

make_all_plots(compare)
