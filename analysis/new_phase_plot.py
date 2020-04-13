import yt

import matplotlib.pylab as plt
from matplotlib.colors import LogNorm, ListedColormap

import numpy as np
import seaborn as sns

import palettable
import plotting_tools as pt
import yt_functions as ytf

def get_log_phase_data(ds, xfield = 'density', yfield = 'temperature', z_min = 0.8, z_max = 1.2):
    ad = ds.all_data()
    z_abs_code = np.abs(ad[('gas', 'z')] / ds.length_unit.in_units('kpc'))
    zmask = (z_abs_code >= z_min) & (z_abs_code <= z_max)

    xfield = ad[('gas', xfield)][zmask]# / 1e27
    yfield = ad[('gas', yfield)][zmask]# / 1e-27
    mass = ad[('gas', 'cell_mass')].in_units('Msun')[zmask]

    logx = np.log10(xfield)
    logy = np.log10(yfield)

    return logx, logy, mass


def constant_T_line(T):
    p_list = []
    entropy_list = []
    rho_list = np.linspace(-30, -24, 100)
    rho_list = np.power(10, rho_list)
    for rho in rho_list:
        p_list.append(pt.calculate_pressure(rho, T))
        entropy_list.append(pt.calculate_entropy(rho, T))
    return np.array(p_list), np.array(entropy_list)
                      
def constant_rho_line(rho):
    p_list = []
    entropy_list = []
    T_list = np.linspace(3, 8, 100)
    T_list = np.power(10, T_list)
    for T in T_list:
        p_list.append(pt.calculate_pressure(rho, T))
        entropy_list.append(pt.calculate_entropy(rho, T))
    return np.array(p_list), np.array(entropy_list)

def make_plot(sim, xfield = 'density', yfield = 'temperature', weighted = True, normed = True, 
              label = None,  z_min = 0.8, z_max = 1.2, 
              work_dir = '../../simulations', plot_dir = '../../plots', sim_fam = 'production'):
    logx_list = np.array([])
    logy_list = np.array([])
    mass_list = np.array([])
    
    output_list = np.arange(40, 61, 1)
#    output_list = [40]
    for output in output_list:
        ds = ytf.load('%s/%s/%s/DD%04d/DD%04d'%(work_dir, sim_fam, sim, output, output))
        logx, logy, mass = get_log_phase_data(ds, xfield = xfield, yfield = yfield, z_min = z_min, z_max = z_max)
        logx_list = np.append(logx_list, logx)
        logy_list = np.append(logy_list, logy)
        mass_list = np.append(mass_list, mass)
#    xedges = np.linspace(-.45, .15, 100)
#    yedges = np.linspace(-1.3, .5, 100)

    if xfield == 'density':
        xlims = (-28, -26)
    elif xfield == 'pressure':
        xlims = (-14.55, -13.01)
    if yfield == 'temperature':
        ylims = (4.3, 7)
    elif yfield == 'entropy':
        ylims = (-1.2, 2.4)

    xedges = np.linspace(xlims[0], xlims[1], 100)
    yedges = np.linspace(ylims[0], ylims[1], 100)
    if weighted:
        weights = mass_list
    else:
        weights = None
    H, xedges, yedges = np.histogram2d(logx_list, logy_list, bins = (xedges, yedges), 
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

    ax_scatter.set_xlim(xlims[0], xlims[1])
    ax_scatter.set_ylim(ylims[0], ylims[1])
    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())
    if xfield == 'density':
        ax_scatter.set_xlabel('Log Density (g cm$^{-3}$)')
    elif xfield == 'pressure':
        ax_scatter.set_xlabel('Log Pressure ($\\frac{\\mathrm{dyn}}{\\mathrm{cm}^2}$)')
    if yfield == 'temperature':
        ax_scatter.set_ylabel('Log Temperature (K)')
    elif yfield == 'entropy':
        ax_scatter.set_ylabel('Log Entropy ($\\mathrm{cm}^2 \cdot \\mathrm{keV}$)')


    palette = sns.color_palette("Blues")
    cmap = ListedColormap(palette)
    background = palette[0]

    ax_scatter.pcolormesh(X, Y, H, norm = LogNorm(), cmap = cmap)
    ax_scatter.set_facecolor(background)
    

    # add constant T line:
    T_list = np.linspace(4, 7, 6)
    rho_list = np.linspace(-28, -26, 6)
    for T in T_list:
        T = np.power(10, T)
        pT, eT = constant_T_line(T)
        ax_scatter.plot(np.log10(pT), np.log10(eT),linestyle = (0, (5, 5)), color = 'black', alpha = 0.4)
    for rho in rho_list:
        rho = np.power(10, rho)
        pT, eT = constant_rho_line(rho)
        ax_scatter.plot(np.log10(pT), np.log10(eT),linestyle = (0, (5, 5)), color = 'black', alpha = 0.4)

    if label is not None:
        xtext = 0.08*(xlims[1] - xlims[0]) + xlims[0]
        ytext = 0.9*(ylims[1] - ylims[0]) + ylims[0]
        ax_scatter.text(xtext, ytext, label, fontsize = 16, color = 'black')
    ax_histx.hist(logx_list, bins=xedges, weights = weights, normed = normed, color = palette[2])
    ax_histx.hist(logx_list, bins=xedges, weights = weights, normed = normed, 
                  histtype = 'step', color = palette[-1])
    ax_histy.hist(logy_list, bins=yedges, weights = weights, normed = normed, 
                  orientation = 'horizontal', color = palette[2])
    ax_histy.hist(logy_list, bins=yedges, weights = weights, normed = normed, 
                  orientation = 'horizontal', histtype = 'step', color = palette[-1])


    #fig.tight_layout()
    figname = '%s/%s/new_phase_%s_%s_%s'%(plot_dir, sim_fam, xfield, yfield, sim)
    if weighted:
        figname += '_weighted'
    print('%s.png'%figname)
    plt.savefig('%s.png'%figname, dpi = 300)


compare = 'transport'
xfield = 'density'
yfield = 'temperature'

xfield = 'pressure'
yfield = 'entropy'
for cr in [1, 0.1]:
    for tctf in [0.3, 1]:
        sim_list = pt.generate_sim_list(compare, tctf = tctf, cr = cr)
        label_list = pt.generate_label_list(compare, tctf = tctf, cr = cr)
        for sim, label in zip(sim_list, label_list):
            make_plot(sim, xfield = xfield, yfield = yfield, label = label)

