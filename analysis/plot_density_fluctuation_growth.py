import yt
from yt import YTArray
from yt import YTQuantity
import sys
import os
import numpy as np
import matplotlib.pylab as plt
import palettable

import plotting_tools as pt

def calculate_rms_fluctuation(sim_folder, output_list, field = 'density', grid_rank = 3):
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
                print(zstart, zend)

                if (grid_rank == 3):
                    region1 = ds.r[0, 0, zstart:zend]
                    region2 = ds.r[0, 0, -zend:-zstart]
                    zlist    = np.append(region1[('gas', 'z')].in_units('kpc'),\
                                         region2[('gas', 'z')].in_units('kpc'))
                
                    dzfield = np.array([])
                    for z in zlist:
                        zslice  = ds.r[:, :, YTQuantity(z, 'kpc')]
                        zfield     = zslice[('gas', field)]
                        zfield_ave = np.mean(zfield)
                        dzfield    = np.append(dzfield, (zfield - zfield_ave) / zfield_ave)
                elif (grid_rank == 2):
                    all_y = ds.ortho_ray('y', (0, 0))[('gas', 'y')].in_units('kpc')
                    ymask = np.abs(all_y / ds.length_units('kpc') > zstart) & np.abs(all_y / ds.length_units('kpc') < zend)
                    ylist = all_y[ymask]
                    for y in ylist:
                        yray = ds.ortho_ray('x', (y, 0))
                        zfield = yray[('gas', field)]
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




def plot_density_fluctuation_growth(sim, beta = 'inf', tctf_list = None, cr_list = None, diff_list = None, \
                                    field = 'density', beta_list = None, work_dir = '../../simulations/'):
    if tctf_list == None:
        tctf_list = [0.1, 0.3, 1.0, 10.0]
    if beta_list == None:
        beta_list = len(tctf_list)*['inf']
    if cr_list == None:
        cr_list = len(tctf_list) * [0]
    if diff_list == None:
        diff_list = len(tctf_list)*[0]
    fig, ax = plt.subplots(figsize = (6, 6))
    ax.set_yscale('log')
    ax.set_ylim(5e-3, 5)
    ax.set_xlim(0, 10)
    
    gamma = 5./3.
    time_list = np.arange(0, 12, 1)
    wcool = 1.0 / (gamma * 1.0)
    pi = (5./3.) * wcool
    ax.plot(time_list, 0.02*np.exp(pi*time_list), color = 'black',\
            linestyle = 'dashed', label = 'Linear Theory', linewidth = 3)

    cpal = palettable.cmocean.sequential.Tempo_7_r.mpl_colors
    cpal = palettable.scientific.sequential.Batlow_11.mpl_colors
    cpal = palettable.scientific.sequential.Batlow_6.mpl_colors

    #output_list = np.linspace(0, 100, 10)
    output_list = np.arange(0, 210, 10)
    for i, tctf in enumerate(tctf_list):
        sim_location = pt.get_sim_location(sim, tctf, beta_list[i], cr_list[i], diff = diff_list[i], work_dir = work_dir)

        print(sim_location)
        if not os.path.isdir(sim_location):
            continue

        label = pt.get_label_name(compare, tctf, beta_list[i], cr_list[i])

        time_list, dzfield_rms_list = calculate_rms_fluctuation(sim_location, output_list, field = field)
        ax.plot(time_list/tctf, dzfield_rms_list, linewidth = 3, label = label, color = cpal[i])
        
        if tctf > 0:
#             sim_location = '%s/%s_tctf_%.1f_2d_128'%(work_dir, sim, tctf)
#             time_list, dzfield_rms_list = calculate_rms_fluctuation(sim_location, output_list, field = field)
#             ax.plot(time_list/tctf, dzfield_rms_list, linewidth = 2, alpha = 0.7, color = cpal[i], \
#                        linestyle = 'solid', label = label + ', 2D 128')

             sim_location = '%s/%s_tctf_%.1f_2d_256'%(work_dir, sim, tctf)
             time_list, dzfield_rms_list = calculate_rms_fluctuation(sim_location, output_list, field = field)
             ax.plot(time_list/tctf, dzfield_rms_list, linewidth = 2, alpha = 0.7, color = cpal[i], \
                        linestyle = 'dashed', label = label + ', 2D 256')
     
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
    if compare == 'beta':
        plot_name = '../../plots/%s_fluctuation_growth_%s_tctf_%.1f_beta_compare'%(field, sim,tctf_list[0])
    elif compare == 'cr':
        plot_name += '_tctf_%.2f_beta_%.1f_cr'%(tctf_list[0], beta_list[0])
        if diff_list[-1] > 0:
            plot_name += '_diff_%.1f'%diff_list[-1]
    plot_name += '.png'
    print(plot_name)

    plt.savefig(plot_name, dpi = 300)

        

field = 'density'
zstart = 0.8
zend = 1.2
save = False
load = False

crdiff = 0

sim = sys.argv[1]
compare = sys.argv[2]

if sys.argv[3]:
    tctf = float(sys.argv[3])

print(compare)
tctf_list, beta_list, cr_list, diff_list = pt.generate_lists(compare, tctf, crdiff = crdiff)
beta_list = len(tctf_list)*[10.0]
print(tctf_list, beta_list, cr_list, diff_list)

plot_density_fluctuation_growth(sim, tctf_list = tctf_list, beta_list = beta_list, \
                                cr_list = cr_list, diff_list = diff_list, field = field)
