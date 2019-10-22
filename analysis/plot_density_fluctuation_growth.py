import yt
from yt import YTArray
from yt import YTQuantity
import sys
import numpy as np
import matplotlib.pylab as plt
import palettable


def calculate_drho_rms(sim_folder, output_list):
    time_list     = np.array([])
    drho_rms_list = np.array([])

    for output in output_list:
        ds      = yt.load("%s/DD%04d/DD%04d"%(sim_folder, output, output))
        region1 = ds.r[0, 0, 0.9:1.1]
        region2 = ds.r[0, 0, -1.1:-0.9]
        zlist    = np.append(region1[('gas', 'z')].in_units('kpc'),\
                         region2[('gas', 'z')].in_units('kpc'))
                     
        drho = np.array([])
        for z in zlist:
            zslice  = ds.r[:, :, YTQuantity(z, 'kpc')]
            rho     = zslice[('gas', 'density')]
            rho_ave = np.mean(rho)
            drho    = np.append(drho, (rho - rho_ave) / rho_ave)

            drho_rms = np.sqrt(np.mean(drho**2))

        time_list     = np.append(time_list,    ds.current_time)
        drho_rms_list = np.append(drho_rms_list, drho_rms)
    
    return time_list, drho_rms_list



def plot_density_fluctuation_growth(sim, tctf_list = None, work_dir = '../../simulations'):
    if tctf_list == None:
        tctf_list = [0.1, 0.3, 1.0, 10.0]
    fig, ax = plt.subplots(figsize = (6, 6))
    ax.set_yscale('log')
    ax.set_ylim(5e-3, 5)
    ax.set_xlim(0, 10)
    
    gamma = 5./3.
    time_list = np.linspace(0, 10, 10)
    wcool = 1.0 / (gamma * 1.0)
    pi = (5./3.) * wcool
    ax.plot(time_list, 0.02*np.exp(pi*time_list), color = 'black',\
            linestyle = 'dashed', label = 'Linear Theory', linewidth = 3)

    cpal = palettable.wesanderson.Darjeeling4_5.mpl_colors
    output_list = np.linspace(0, 100, 10)
    for i, tctf in enumerate(tctf_list):
        sim_location = '%s/%s_tctf_%.1f'%(work_dir, sim, tctf)
        label = '$t_{cool}/t_{ff}$ = %0.1f'%tctf
        time_list, drho_rms_list = calculate_drho_rms(sim_location, output_list)
        ax.plot(time_list/tctf, drho_rms_list, linewidth = 3, label = label, color = cpal[i])
        if tctf == 1.0:
            for res, linestyle in zip([64, 256], ['dashed', 'dotted']):
                sim_location = '%s/%s_%i'%(work_dir, sim, res)
                time_list, drho_rms_list = calculate_drho_rms(sim_location, output_list)
                ax.plot(time_list/tctf, drho_rms_list, linewidth = 3, color = cpal[i], \
                        linestyle = linestyle, label = label + ', res = %i$^3$'%res)
    
    ax.set_xlabel('t/t$_{cool}$')
    ax.set_ylabel('$ \\langle \\delta \\rho / \\rho\\rangle_{\\mathrm{rms}}$ ')
    ax.legend()
    fig.tight_layout()
    plt.savefig('../../plots/density_fluctuation_growth_%s.png'%sim, dpi = 300)


tctf_list = [0.1, 0.3, 1.0, 10]
    
sim = sys.argv[1]
plot_density_fluctuation_growth(sim)
