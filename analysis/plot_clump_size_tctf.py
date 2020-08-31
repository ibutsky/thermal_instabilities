import yt
from yt.data_objects.level_sets.api import *

import numpy as np
import matplotlib.pylab as plt

import palettable
import plotting_tools as pt


def _neg_log_T(field, data):
    log_rho = -np.log10(data[('gas', 'temperature')])
    return log_rho

def find_cold_clumps(ds, field = 'neg_log_T'):
    ds.add_field(('gas', 'neg_log_T'), function = _neg_log_T, units = '')
    ad = ds.all_data()
    ad_slice = ds.r[:, :, 0.8:1.2]
    master_clump = Clump(ad, ('gas', field))

    master_clump.add_validator("min_cells", 3)
    master_clump.add_info_item("total_cells")

    c_min = -5
    c_max = ad["gas", field].max()
    step = -0.1
    find_clumps(master_clump, c_min, c_max, step)
    return master_clump

def find_radius(leaf_clumps):
    radius_list = []
    n_clumps = 0
    H = 43.8535772
    z_min = 0.8 * H 
    z_max = 1.2 * H
    for i, clump in enumerate(leaf_clumps):
        z = np.abs(clump[('gas', 'z')].in_units('kpc'))
        z_ave = np.mean(z)
        if (z_ave >= z_min) and (z_ave <= z_max):
            n_clumps += 1
            ncells = len(z)
            radius = np.power(ncells, 1./3.)
            print(i, ncells, radius)
            radius_list.append(radius)

    return n_clumps, np.mean(radius_list), np.std(radius_list)

def plot_clump_size(output = 40, beta = 100.0, cr = 0, work_dir = '../../simulations'):
    tctf_list = [0.1, 0.3, 1]#, 3]

    fig, ax = plt.subplots(nrows=1, ncols = 2, figsize = (8, 3.8), sharex = True)
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    ax[0].set_xlim(.09, 1.5)
    ax[0].set_ylim(1, 20)
    ax[0].set_xlabel('$t_{cool} / t_{ff}$')
    ax[0].set_ylabel('Clump Radius')

    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[1].set_xlim(.09, 1.5)
    ax[1].set_ylim(1, 300)
    ax[1].set_xlabel('$t_{cool} / t_{ff}$')
    ax[1].set_ylabel('Number of Clumps')

    sim_fam_list = ['production', 'production/high_res', 'production']
    cr_list = [0, 0, 0, 0, 1]
    beta_list = [100, 100, 10, 3]
    color_list = ['blue', 'blue', 'green',  'orange']
    linestyle_list = ['solid', 'dashed', 'solid', 'solid']
    label_list = ['No CR', 'High Res', 'Beta = 10','P$_c$ / P$_g$ = 1']
    for i in range(len(sim_fam_list)):
        n_clump_list = []
        r_clump_list = []
        std_clump_list = []
        for tctf in tctf_list:
            sim_loc = pt.get_sim_location('isocool', tctf, beta_list[i], cr = cr_list[i], sim_fam = sim_fam_list[i])
            print(sim_loc)
            ds = yt.load('%s/DD%04d/DD%04d'%(sim_loc, output, output))
            master_clump = find_cold_clumps(ds)
            
            leaves = master_clump.leaves
            N_clump, r_clump, std_clump = find_radius(leaves)
            n_clump_list.append(N_clump)
            r_clump_list.append(r_clump)
            std_clump_list.append(std_clump)


        print(tctf, n_clump_list, r_clump_list, std_clump_list)
        ax[0].plot(tctf_list, r_clump_list, color = color_list[i], linestyle = linestyle_list[i], marker = 'o', label = label_list[i])
        ax[0].errorbar(tctf_list, std_clump_list, color = color_list[i], linestyle = "")
        ax[1].plot(tctf_list, n_clump_list, color = color_list[i], linestyle = linestyle_list[i], marker = 'o')
    ax[0].legend()
    fig.tight_layout()
    plt.savefig('../../plots/production/clump_test_%i.png'%output, dpi = 300)

plot_clump_size(output = 30)
plot_clump_size(output = 50)
plot_clump_size(output = 60)
