import yt
from yt import YTArray
from yt import YTQuantity
import sys
import os
import numpy as np
import matplotlib.pylab as plt
import palettable

import plotting_tools as pt

def calculate_rms_fluctuation(sim_folder, output_list, field = 'density'):
    time_list     = np.array([])
    dzfield_rms_list = np.array([])
    if not os.path.isdir(sim_folder):
        return time_list, dzfield_rms_list
    sim_base = os.path.basename(sim_folder)
    print(sim_base)
    
    out_name = '../../data/fluctuation_growth_%s_%s'%(sim_base, field)
    if os.path.isfile(out_name) and load == True:
        time_list, dzfield_rms_list = np.loadtxt(out_name, unpack=True)
    
    else:

        for output in output_list:
            ds_path = "%s/DD%04d/DD%04d"%(sim_folder, output, output)

            if os.path.isfile(ds_path):
                ds = yt.load(ds_path)
                region1 = ds.r[0, 0, 0.9:1.1]
                region2 = ds.r[0, 0, -1.1:-0.9]
                zlist    = np.append(region1[('gas', 'z')].in_units('kpc'),\
                                     region2[('gas', 'z')].in_units('kpc'))
                
                dzfield = np.array([])
                for z in zlist:
                    zslice  = ds.r[:, :, YTQuantity(z, 'kpc')]
                    zfield     = zslice[('gas', field)]
                    zfield_ave = np.mean(zfield)
                    dzfield    = np.append(dzfield, (zfield - zfield_ave) / zfield_ave)

                dzfield_rms = np.sqrt(np.mean(dzfield**2))
                time_list     = np.append(time_list,    ds.current_time)
                dzfield_rms_list = np.append(dzfield_rms_list, dzfield_rms)

        if save:
            outf = open(out_name, 'w')
            for i in range(len(time_list)):
                outf.write("%e %e\n"%(time_list[i], dzfield_rms_list[i]))
            outf.close()

    
    return time_list, dzfield_rms_list




def plot_density_fluctuation_growth(sim, beta = 'inf', tctf_list = None, cr_list = None, \
                                    field = 'density', beta_list = None, work_dir = '../../simulations'):
    if tctf_list == None:
        tctf_list = [0.1, 0.3, 1.0, 10.0]
    if beta_list == None:
        beta_list = len(tctf_list)*['inf']
    if cr_list == None:
        cr_list = len(tctf_list) * [0]
    fig, ax = plt.subplots(figsize = (6, 6))
    ax.set_yscale('log')
    ax.set_ylim(5e-3, 5)
    if nocool:
        ax.set_ylim(1e-3, 1e-1)
    ax.set_xlim(0, 10)
    
    gamma = 5./3.
    time_list = np.arange(0, 12, 1)
    wcool = 1.0 / (gamma * 1.0)
    if 0:
        pi_high_eta = 4./9 * wcool
        pi_low_eta = 16./15 * wcool
        ax.plot(time_list, 0.02*np.exp(pi_high_eta*time_list), color = 'black',\
            linestyle = 'dashed', label = 'Linear Theory, high $\\eta$', linewidth = 3)
        ax.plot(time_list, 0.02*np.exp(pi_low_eta*time_list), color = 'black',\
            linestyle = 'dotted', label = 'Linear Theory, low $\\eta$', linewidth = 3)
    elif not nocool:
        pi = (5./3.) * wcool
        ax.plot(time_list, 0.02*np.exp(pi*time_list), color = 'black',\
            linestyle = 'dashed', label = 'Linear Theory', linewidth = 3)

    cpal = palettable.cmocean.sequential.Tempo_7_r.mpl_colors
    cpal = palettable.scientific.sequential.Batlow_11.mpl_colors
    cpal = palettable.scientific.sequential.Batlow_6.mpl_colors

    #output_list = np.linspace(0, 100, 10)
    output_list = np.arange(0, 210, 10)
    for i, tctf in enumerate(tctf_list):
        if beta_list[i] == 'inf':
            sim_location = '%s/%s_tctf_%.1f'%(work_dir, sim, tctf)
        else:
            sim_location =  '%s/%s_tctf_%.1f_beta_%.1f'%(work_dir, sim, tctf, beta_list[i])
        cr = cr_list[i]
        if cr > 0:
            if cr < 0.1:
                sim_location += '_cr_%.2f'%(cr)
            else:
                sim_location += '_cr_%0.1f'%(cr)
        if nocool:
            sim_location += '_nocool'
        print(sim_location)
        if not os.path.isdir(sim_location):
            continue

        label = '$t_{cool}/t_{ff}$ = %.1f'%tctf
        if cr_compare:
            label = 'P$_c$ / P$_g$ = %.1f'%cr
            if cr < 0.1:
                label = 'P$_c$ / P$_g$ = %.2f'%cr
        if beta_compare:
            if beta_list[i] == 'inf':
                label = 'Hydro'
            else:
                label = '$\\beta = $%.1f'%beta_list[i]
        time_list, dzfield_rms_list = calculate_rms_fluctuation(sim_location, output_list, field = field)
        ax.plot(time_list/tctf, dzfield_rms_list, linewidth = 3, label = label, color = cpal[i])

        if tctf == 1.0 and beta == 'inf' and beta_compare == 0 and cr_compare == 0:
            for res, linestyle in zip([64, 256], ['dashed', 'dotted']):
                sim_location = '%s/%s_%i'%(work_dir, sim, res)
                time_list, dzfield_rms_list = calculate_rms_fluctuation(sim_location, output_list, field = field)
                ax.plot(time_list/tctf, dzfield_rms_list, linewidth = 1, alpha = 0.7, color = cpal[i], \
                        linestyle = linestyle, label = label + ', res = %i$^3$'%res)
    if do_test:
        time_list, dzfield_rms_list = calculate_rms_fluctuation(test_sim_location, output_list, field = field)
        ax.plot(time_list/tctf_test, dzfield_rms_list, linewidth = 3, color = 'red', label = 'test')
    
    ax.set_xlabel('t/t$_{cool}$')
    if field == 'density':
        ax.set_ylabel('$ \\langle \\delta \\rho / \\rho\\rangle_{\\mathrm{rms}}$ ')
    elif field == 'temperature':
        ax.set_ylabel('RMS Temperature Fluctuation')
    ax.legend()
    fig.tight_layout()
    if beta == 'inf':
        plot_name = '../../plots/%s_fluctuation_growth_%s'%(field, sim)
    else:
        plot_name = '../../plots/%s_fluctuation_growth_%s_beta_%.1f'%(field, sim, beta)
    if beta_compare:
        plot_name = '../../plots/%s_fluctuation_growth_%s_tctf_%.1f_beta_compare'%(field, sim,tctf_list[0])
    if cr_compare:
        plot_name += '_tctf_%.2f_beta_%.1f_cr'%(tctf_list[0], beta_list[0])
    plot_name += '.png'
    print(plot_name)

    plt.savefig(plot_name, dpi = 300)


do_test = False
test_sim_location = '../../simulations/isocool_tctf_0.1_beta_3.0_constB'
tctf_test = 0.1

tctf_list = [0.1, 0.3, 1.0, 3.0, 10]
cr_list = None

#cr_list = [0.01, 0.1, 0.3, 0.6, 0.9, 1.0, 1.1, 2.0, 3.0, 10.0, 30.0]
cr_list = [0.1, 0.3, 1.0, 3.0, 10.0]
cr_list = [0.01, 0.1, 1.0, 10.0, 100.0]
tctf_list = len(cr_list) * [3.0]
beta_list = len(cr_list)* [10.0]

tctf_list = 6*[3.0]
cr_list = [0, 0, 0, 0, 0, 0]
beta_list = [3, 10, 30, 100, 300, 'inf']
#beta_list = [10, 'inf']

sim = 'isocool_isochoric'
sim = 'isocool'
#sim = 'isothermal'
nocool = False

beta_compare = 1
cr_compare = 0
field = 'density'
#field = 'temperature'
save = False
load = True 


plot_density_fluctuation_growth(sim, tctf_list = tctf_list, beta_list = beta_list, cr_list = cr_list, field = field)
