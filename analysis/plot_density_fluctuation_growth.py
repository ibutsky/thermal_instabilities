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
                else:
                    all_y = ds.ortho_ray('y', (0, 0))[('gas', 'y')].in_units('kpc')
                    ymask = np.abs(all_y / ds.length_unit.in_units('kpc') > zstart) & np.abs(all_y / ds.length_unit.in_units('kpc') < zend)
                    ylist = all_y[ymask]
                    dzfield = np.array([])
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




def plot_density_fluctuation_growth(sim, tctf_list = [1.0], cr_list = [0], diff_list = [0], \
                                    stream_list = [0], heat_list = [0], field = 'density', \
                                    beta_list = ['inf'], work_dir = '../../simulations/', grid_rank = 3):

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
        sim_location = pt.get_sim_location(sim, tctf, beta_list[i], cr_list[i], \
                                           diff = diff_list[i], stream = stream_list[i], 
                                           heat = heat_list[i], work_dir = work_dir)

        if not os.path.isdir(sim_location):
            continue

        label = pt.get_label_name(compare, tctf, beta_list[i], cr_list[i], crdiff = diff_list[i], \
                       crstream = stream_list[i], crheat = heat_list[i])

        time_list, dzfield_rms_list = calculate_rms_fluctuation(sim_location, output_list, field = field, grid_rank = grid_rank)
        ax.plot(time_list/tctf, dzfield_rms_list, linewidth = 3, label = label, color = cpal[i])
        
        if tctf < 0:
            for res in [256]:
                sim_location = '%s/%s_tctf_%.1f'%(work_dir, sim, tctf)
                if beta_list[i] != 'inf':
                    sim_location += '_beta_%.1f'%beta_list[i]
                sim_location += '_2d_%i'%(res)
                time_list, dzfield_rms_list = calculate_rms_fluctuation(sim_location, output_list, field = field, grid_rank = 2)
                print(time_list, dzfield_rms_list)
                ax.plot(time_list/tctf, dzfield_rms_list, linewidth = 2, alpha = 0.7, color = cpal[i], \
                        linestyle = 'dashed', label = label + ', 2D 256')
     
    ax.set_xlabel('t/t$_{cool}$')
    if field == 'density':
        ax.set_ylabel('$ \\langle \\delta \\rho / \\rho\\rangle_{\\mathrm{rms}}$ ')
    elif field == 'temperature':
        ax.set_ylabel('RMS Temperature Fluctuation')
    ax.legend()
    fig.tight_layout()
    figname = pt.get_fig_name('density_fluctuation_growth', sim, compare, tctf_list[0], beta_list[0], \
                              cr_list[0], diff_list[0])
    plt.savefig(figname, dpi = 300)

        

work_dir = '../../simulations/'
grid_rank = 3

field = 'density'
zstart = 0.8
zend = 1.2
save = True
load = False

crdiff = 0

sim = sys.argv[1]
compare = sys.argv[2]

if sys.argv[3]:
    tctf = float(sys.argv[3])


print(compare)
tctf_list, beta_list, cr_list, diff_list, stream_list, heat_list\
    = pt.generate_lists(compare, tctf, crdiff = crdiff)

print(tctf_list, beta_list, cr_list, diff_list)

plot_density_fluctuation_growth(sim, tctf_list = tctf_list, beta_list = beta_list, \
                                grid_rank = grid_rank, stream_list = stream_list, \
                                heat_list = heat_list, cr_list = cr_list, diff_list = diff_list, \
                                field = field, work_dir = work_dir)
