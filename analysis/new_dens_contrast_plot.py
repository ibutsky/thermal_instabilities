import numpy as np
import matplotlib.pylab as plt
import palettable
import seaborn as sns

import plotting_tools as pt

def get_marker(transport):
    if transport == 1:
        marker = 'o'
    elif transport == 2:
        marker = 's'
    elif transport == 3:
        marker = 'v'
    return marker

def plot_density_fluctuation(sim_fam = 'production'):

    fig, ax = plt.subplots(nrows=1, ncols = 1, figsize = (4.4, 4), sharex = True, sharey = False)
    ax.set_xlim(-2, 2.5)
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlabel('Log $(P_c / P_g)_{\mathrm{cold}}$')
    ax.set_ylabel('Log Average Density Contrast')

    pal = sns.cubehelix_palette(8, start=.5, rot=-.75)
    pal = pt.get_color_list('transport_relative')
    color_list = [pal[0], pal[3], pal[6]]
    label_list = ['Fiducial', 'High-res', 'Low-res', 'Tmin = 1e4']
    marker = 'o'

    all_data = np.load('../../data/%s/dens_contrast_plot_data.npy'%sim_fam)#, allow_pickle = True)
    rho_contrast, rho_err, creta, creta_err, transport_list = all_data

    for transport in range(3):
        i = transport # defining redundantly for now
        mask = transport_list == transport
        
        ax.errorbar(creta[mask], rho_contrast[mask], xerr = creta_err[mask], yerr = rho_err[mask], fmt = 'none', 
                                     color = color_list[i],linewidth = 0.5, alpha = 0.2, zorder = 1)
        ax.scatter(creta[mask], rho_contrast[mask], color = color_list[i], #facecolors = 'none',
                       label = None, marker = get_marker(transport), alpha = 1, zorder = 2)

        # just for the label:
        #ax.scatter(-100, -100, color = color_list[i], label = label_list[i], marker = marker, alpha = 0.9)


    #just for the label
    ax.scatter(-100, -100, color = color_list[0], facecolors = 'none', marker = 'o', label = 'Advection')
    ax.scatter(-100, -100, color = color_list[0], facecolors = 'none', marker = 's', label = 'Diffusion')
    ax.scatter(-100, -100, color = color_list[0], facecolors = 'none', marker = 'v', label = 'Streaming')

    # analytic line
    x = np.linspace(0, 3, 20)

    y1 = -3./4. *  x + 2
    y2 = -3./4. *  x + 3.5
    if warm:
        y = 3./2 * np.log10(1 + 1/np.power(10, x))
#    ax.plot(x, y, color = 'black', linestyle = 'dotted') 
    ax.fill_between(x, y1, y2 = y2, color = 'gray', alpha = 0.4, zorder = 0)

    ax.legend(fontsize = 8, ncol = 2, loc = 3)
    fig.tight_layout()
    fig_basename = 'density_contrast'

    figname = '../../plots/%s/%s.png'%(sim_fam, fig_basename)
    print(figname)
    plt.savefig(figname, dpi = 300)

sim_fam = 'production'
plot_density_fluctuation(sim_fam = sim_fam)

