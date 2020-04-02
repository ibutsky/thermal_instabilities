import yt

import matplotlib.pylab as plt
from matplotlib.colors import LogNorm, ListedColormap

import numpy as np
import seaborn as sns

import palettable
import plotting_tools as pt

def get_log_T_rho_data(ds, z_min = 0.8, z_max = 1.2):
    ad = ds.all_data()
    z_abs_code = np.abs(ad[('gas', 'z')] / ds.length_unit.in_units('kpc'))
    zmask = (z_abs_code >= z_min) & (z_abs_code <= z_max)

    temperature = ad[('gas', 'temperature')][zmask]# / 1e6
    density = ad[('gas', 'density')][zmask]# / 1e-27
    mass = ad[('gas', 'cell_mass')].in_units('Msun')[zmask]

    logt = np.log10(temperature)
    logd = np.log10(density)

    return logt, logd, mass

def make_plot(sim, weighted = True, normed = True, z_min = 0.8, z_max = 1.2, 
              work_dir = '../../simulations', plot_dir = '../../plots', sim_fam = 'production'):
    logt = np.array([])
    logd = np.array([])
    mass = np.array([])
    
    output_list = np.arange(40, 61, 1)
#    output_list = [40]
    for output in output_list:
        ds = yt.load('%s/%s/%s/DD%04d/DD%04d'%(work_dir, sim_fam, sim, output, output))
        t, d, m = get_log_T_rho_data(ds, z_min = z_min, z_max = z_max)
        logt = np.append(logt, t)
        logd = np.append(logd, d)
        mass = np.append(mass, m)
#    xedges = np.linspace(-.45, .15, 100)
#    yedges = np.linspace(-1.3, .5, 100)
    xedges = np.linspace(-28, -26, 100)
    yedges = np.linspace(4.3, 7, 100)
    if weighted:
        weights = mass
    else:
        weights = None
    H, xedges, yedges = np.histogram2d(logd, logt, bins = (xedges, yedges), 
                                       weights = weights, normed = normed)
    H = H.T
    X, Y = np.meshgrid(xedges, yedges)


    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    fig = plt.figure(figsize=(6, 6))

    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True)
    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(direction='in', labelbottom=False)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False)
    

#    fig, ax = plt.subplots(figsize=(4, 4))
#    ax.set_ylim(-1.3, .5)
#    ax.set_xlim(-.45, .15)
    ax_scatter.set_xlim(-28, -26)
    ax_scatter.set_ylim(4.3, 7)
    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())
    ax_scatter.set_xlabel('Log Density (g cm$^{-3}$)')
    ax_scatter.set_ylabel('Log Temperature (K)')


    palette = sns.color_palette("Blues")
#    palette = sns.cubehelix_palette(start=2.8, rot=.1, n_colors = 8)
#    palette = palettable.cartocolors.sequential.Mint_7
#    palette = palettable.cmocean.sequential.Ice_20_r.mpl_colors
    cmap = ListedColormap(palette)
    background = palette[0]
#    cmap = palette.mpl_colormap
#    background = palette.mpl_colors[0]
    ax_scatter.pcolormesh(X, Y, H, norm = LogNorm(), cmap = cmap)
    ax_scatter.set_facecolor(background)
    ax_histx.hist(logd, bins=xedges, weights = weights, normed = normed, color = palette[2])
    ax_histx.hist(logd, bins=xedges, weights = weights, normed = normed, 
                  histtype = 'step', color = palette[-1])
    ax_histy.hist(logt, bins=yedges, weights = weights, normed = normed, 
                  orientation = 'horizontal', color = palette[2])
    ax_histy.hist(logt, bins=yedges, weights = weights, normed = normed, 
                  orientation = 'horizontal', histtype = 'step', color = palette[-1])


    #fig.tight_layout()
    figname = '%s/%s/new_phase_%s'%(plot_dir, sim_fam, sim)
    if weighted:
        figname += '_weighted'
    plt.savefig('%s.png'%figname, dpi = 300)


compare = 'transport'
for cr in [0.1, 1]:
    for tctf in [0.3, 1]:
        sim_list = pt.generate_sim_list(compare, tctf = tctf, cr = cr)
        for sim in sim_list:
            make_plot(sim)

