import yt
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm, ListedColormap
import astropy.constants as const

import numpy as np
import seaborn as sns

import palettable
import plotting_tools as pt
import yt_functions as ytf

def make_plot(sim, xfield = 'density', yfield = 'temperature', weighted = True, normed = True, 
              label = None,  z_min = 0.8, z_max = 1.2, number_density = True,
              work_dir = '../../simulations', plot_dir = '../../plots', sim_fam = 'production'):


    logx_list, logy_list, mass_list = pt.get_2d_hist_data(xfield, yfield, sim, zstart = z_min, zend = z_max, 
                                                          work_dir = work_dir, sim_fam = sim_fam)

    mu = 1.22
    mh = const.m_p.cgs.value
    log_mumh = np.log10(mu*mh)

    if xfield == 'density':
        xlims = (-28, -26)
        if sim_fam == 'production/Tmin1e4':
            xlims = (-28, -24)
        if number_density:
            xlims -= log_mumh
            logx_list -= log_mumh

    elif xfield == 'pressure':
        xlims = (-14.55, -13.01)
    if yfield == 'temperature':
        ylims = (4.3, 6.9)
        if sim_fam == 'production/Tmin1e4':
            ylims = (3.5, 6.9)
    elif yfield == 'entropy':
        ylims = (-1.2, 2.4)

    xedges = np.linspace(xlims[0], xlims[1], 100)
    yedges = np.linspace(ylims[0], ylims[1], 100)
    if weighted:
        weights = mass_list
    else:
        weights = None
    H, xedges, yedges = np.histogram2d(logx_list, logy_list, bins = (xedges, yedges), weights = weights, normed = normed)
    H = H.T
    X, Y = np.meshgrid(xedges, yedges)


    # definitions for the axes
    left, width = 0.12, 0.65
    bottom, height = 0.11, 0.65
    spacing = 0.005
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    fig = plt.figure(figsize=(7, 5.9))
    fs = 14
    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True, labelsize = fs)
    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(direction='in', labelbottom=False, labelsize = fs)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False, labelsize = fs)

    ax_scatter.set_xlim(xlims[0], xlims[1])
    ax_scatter.set_ylim(ylims[0], ylims[1])
    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())
    if xfield == 'density':
        ax_scatter.set_xlabel('Log Density (g cm$^{-3}$)', fontsize = fs)
        if number_density:
            ax_scatter.set_xlabel('Log Number Density (cm$^{-3}$)', fontsize = fs)

    elif xfield == 'pressure':
        ax_scatter.set_xlabel('Log Pressure ($\\frac{\\mathrm{dyn}}{\\mathrm{cm}^2}$)', fontsize = fs)
    if yfield == 'temperature':
        ax_scatter.set_ylabel('Log Temperature (K)', fontsize = fs)
    elif yfield == 'entropy':
        ax_scatter.set_ylabel('Log Entropy ($\\mathrm{cm}^2 \cdot \\mathrm{keV}$)', fontsize = fs)


    palette = sns.color_palette("Blues")
    cmap = ListedColormap(palette)
    background = palette[0]

    ax_scatter.pcolormesh(X, Y, H, norm = LogNorm(), cmap = cmap)
    ax_scatter.set_facecolor(background)
    

    # add constant T line:
    if xfield == 'pressure':
        T_list = np.linspace(4, 7, 6)
        rho_list = np.linspace(-28, -26, 6)
        for T in T_list:
            T = np.power(10, T)
            pT, eT = pt.constant_T_line(T)
            ax_scatter.plot(np.log10(pT), np.log10(eT),linestyle = (0, (5, 5)), color = 'black', alpha = 0.4)
        for rho in rho_list:
            rho = np.power(10, rho)
            pT, eT = pt.constant_rho_line(rho)
            ax_scatter.plot(np.log10(pT), np.log10(eT),linestyle = (0, (5, 5)), color = 'black', alpha = 0.4)
    elif xfield == 'density':
        e_list = np.linspace(-1.5, 3.0, 6)
        p_list = np.linspace(-15.5, -11, 6)
        color_list = ['green', 'blue', 'purple', 'black', 'orange', 'red']
        e_label_pos = [(-4, 4.5), (-4, 4.5), (-6, 4.7), (-6, 5.5), (-6, 6.5), (-6, 6.5)]
        p_label_pos = [(-6.2, 4.7), (-5.1, 5.1), (-4, 5.3), (-3.5, 5.5), (-3, 6), (-3, 6.5)]
        e_index = [0, 23, 43, 43, 43, 0]
        p_index = [49, 59, 59, 59, 69, 0]
        for e, pos, color, i in zip(e_list, e_label_pos, color_list, e_index):
            e = np.power(10, e)
            rho, T = pt.constant_entropy_line(e)
            log_rho = np.log10(rho)
            log_T = np.log10(T)
            if number_density:
                log_rho -= log_mumh
            ax_scatter.plot(log_rho, log_T, linestyle = (0, (5, 5)), color = 'black', alpha = 0.4)
            angle = 24 #27
            ax_scatter.text(log_rho[i]+.06, log_T[i]-.06, '$\\kappa = 10^{%.1f}\mathrm{cm}^2\mathrm{keV}$'%np.log10(e), 
                            rotation = angle, rotation_mode = 'anchor', 
                            color = 'black', alpha = 0.4, fontsize = fs-6)
        for p, pos, color, i in zip(p_list, p_label_pos, color_list, p_index):
            p = np.power(10, p)
            rho, T = pt.constant_pressure_line(p)
            log_rho = np.log10(rho)
            log_T = np.log10(T)
            if number_density:
                log_rho -= log_mumh
            ax_scatter.plot(log_rho, log_T,linestyle = (0, (5, 5)), color = 'black', alpha = 0.4)
            # label lines of constant pressure
            angle = -33 #30
            ax_scatter.text(log_rho[i]+.01, log_T[i]+.01, '$\\mathrm{P} = 10^{%.1f}\mathrm{dyn/cm}^2$'%np.log10(p), 
                            rotation = angle, rotation_mode = 'anchor', 
                            color = 'black', alpha = 0.4, fontsize = fs-6)
            

    if label is not None:
        xtext = 0.08*(xlims[1] - xlims[0]) + xlims[0]
        ytext = 0.9*(ylims[1] - ylims[0]) + ylims[0]
        ax_scatter.text(xtext, ytext, label, fontsize = fs+2, color = 'black')
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

sim_fam = 'production'

compare = 'transport_pdf'
xfield = 'density'
yfield = 'temperature'


for cr in [0]:
    for tctf in [0.3, 1.0]:
        sim_list = pt.generate_sim_list(compare, tctf = tctf, cr = cr)
        label_list = pt.generate_label_list(compare, tctf = tctf, cr = cr)
        for sim, label in zip(sim_list, label_list):
            make_plot(sim, xfield = xfield, yfield = yfield, label = label, sim_fam = sim_fam)

