from scipy import optimize
from sklearn.neighbors import KernelDensity
import numpy as np
import matplotlib.pylab as plt
import palettable
import seaborn as sns

import plotting_tools as pt

def get_marker(transport):
    if transport  < 4:
        marker = 'o'
    elif transport < 7:
        marker = 's'
    elif transport < 13:
        marker = 'v'
    return marker

def get_color(transport):
    pal = pt.get_color_list('transport_relative')
    if transport < 10:
        color = pal[transport - 1]
    else:
        color = pal[transport - 4]
    return color

def solve_dcontrast(a = 4./3., X = 20, inv_beta = 0, b = 4./3.):
    etas = np.logspace(-3, 3, 1000)
    #X = 20 # Thot/Tcold
    func = lambda d,eta: X/d + eta/d**(a) + inv_beta/d**(b) - (1+eta+inv_beta)
    d_solutions = np.array([optimize.root_scalar(func,args=(eta), bracket=[1,1000]).root for eta in etas])
    return np.log10(etas), np.log10(d_solutions)
        

def calculate_median_profile(x, y,  xmin = None, xmax = None, nbins = 50, centered = True):
    if xmin == None:
        xmin = x.min()
    if xmax == None:
        xmax = x.max()

    xbins           = np.linspace(xmin, xmax, nbins+1)
    centered_x_bins = (xbins + (xmax/nbins/2.0))
    median          = np.zeros(nbins)
    mean            = np.zeros(nbins)
    std             = np.zeros(nbins)

    for i in range(nbins):
        mask = (x >= xbins[i]) & (x < xbins[i+1])
        sample = y[mask]
        median[i] = np.median(sample)
        mean[i]   = np.mean(sample)
        std[i] = np.std(sample)

    if centered:
        xbins = centered_x_bins
    return xbins[:-1], median, mean, std


def generate_sample_data(x, y, n = 10000, bandwidth = 0.1):
    data = np.vstack([x, y]).T
    kde = KernelDensity(kernel = 'gaussian', bandwidth = bandwidth).fit(data)
    sample = kde.sample(n)
    return sample

def generate_fit(x, y, n = 100000, bandwidth = 0.1, nbins = 15, xmin = -1, xmax = 2.5):
    data = np.vstack([x, y]).T
    kde = KernelDensity(kernel = 'gaussian', bandwidth = bandwidth).fit(data)
    sample = kde.sample(n)

    x_fit, median, mean, std = calculate_median_profile(sample[:, 0], sample[:, 1],
                                                        xmin = xmin, xmax = xmax, nbins = nbins)
    return sample, x_fit, median, mean, std


def plot_fit(ax, creta, rho_contrast, transport_list):
    #advection
    mask = (transport_list < 4) & (rho_contrast > -99)                                                                                                                      
    sample, xfit, median, mean, std = generate_fit(creta[mask], rho_contrast[mask],
                                                   bandwidth = .1, nbins = 5, xmin = -2, xmax = 3)
  #  std *= 0.5
    ax.fill_between(xfit, median-std, median+std, color = get_color(2), alpha =0.4)
    ax.plot(xfit, median, color = get_color(2), zorder = 10)
    ax.fill_between(xfit, median-std, median+std, color = get_color(2), alpha =0.4)

    # diffusion                                                                                                                                                             
    mask = (transport_list > 3) & (transport_list < 7) & (rho_contrast > -99)
    sample, xfit, median, mean, std = generate_fit(creta[mask], rho_contrast[mask],
                                                   bandwidth = .1, nbins = 5, xmin = -2, xmax = 3)
 #   std *= 0.5
    ax.fill_between(xfit, median-std, median+std, color = get_color(5), alpha =0.4)
    ax.plot(xfit, median, color = get_color(5), zorder = 10)
    ax.fill_between(xfit, median-std, median+std, color = get_color(5), alpha =0.4)

    #streaming
    mask = (transport_list > 6) & (rho_contrast > -99)
    sample, xfit, median, mean, std = generate_fit(creta[mask], rho_contrast[mask],
                                                   bandwidth = .1, nbins = 5, xmin = -2, xmax = 3)
#    std *= 0.5
    ax.fill_between(xfit, median-std, median+std, color = get_color(5), alpha =0.4)
    ax.plot(xfit, median, color = get_color(8), zorder = 10)
    ax.fill_between(xfit, median-std, median+std, color = get_color(8), alpha =0.4)

def load_plot_data(sim_fam_list = ['production', 'high_res']):
    rho_contrast = np.array([])
    rho_err = np.array([])
    creta = np.array([])
    creta_err = np.array([])
    transport_list = np.array([])
    
    for sim_fam in sim_fam_list:
        all_data = np.load('../data/%s/dens_contrast_plot_data.npy'%sim_fam)
        rho_contrast = np.append(rho_contrast, all_data[0])
        rho_err = np.append(rho_err, all_data[1])
        creta = np.append(creta, all_data[2])
        creta_err = np.append(creta_err, all_data[3])
        transport_list = np.append(transport_list, all_data[4])

    return rho_contrast, rho_err, creta, creta_err, transport_list
    
def plot_density_fluctuation(sim_fam = 'production'):

    fig, ax = plt.subplots(nrows=1, ncols = 1, figsize = (6, 4), sharex = True, sharey = False)
    ax.set_xlim(-2, 2.5)
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlabel('Log $(P_c / P_g)_{\mathrm{cold}}$')
    ax.set_ylabel('Log Average Density Contrast')

    sim_fam_list = ['production', 'high_res']
    rho_contrast, rho_err, creta, creta_err, transport_list = load_plot_data(sim_fam_list = sim_fam_list)

    for transport in range(1, 13):
        i = transport -1
        mask = (transport_list == transport) & (rho_contrast > -99) 
        
        ax.errorbar(creta[mask], rho_contrast[mask], xerr = creta_err[mask], yerr = rho_err[mask], fmt = 'none', 
                                     color = get_color(transport),linewidth = 0.8, alpha = 0.2, zorder = 1)
        ax.scatter(creta[mask], rho_contrast[mask], color = get_color(transport), #facecolors = 'none',
                       label = None, marker = get_marker(transport), alpha = 0.8, zorder = 2)


    #just for the label                                                                                           
    ax.scatter(-100, -100, color = get_color(1), marker = 'o', label = 'Advection, $ \\beta = 100 $')
    ax.scatter(-100, -100, color = get_color(2), marker = 'o', label = 'Advection, $ \\beta = 10 $')
    ax.scatter(-100, -100, color = get_color(3), marker = 'o', label = 'Advection, $ \\beta = 3 $')

    ax.scatter(-100, -100, color = get_color(4), marker = 's', label = 'Diffusion, $t_{\\rm diff}/t_{\\rm ff} = 10 $')
    ax.scatter(-100, -100, color = get_color(5), marker = 's', label = 'Diffusion, $t_{\\rm diff}/t_{\\rm ff} = 3 $')
    ax.scatter(-100, -100, color = get_color(6), marker = 's', label = 'Diffusion, $t_{\\rm diff}/t_{\\rm ff} = 1 $')
    
    ax.scatter(-100, -100, color = get_color(7), marker = 'v', label = 'Streaming, $ \\beta = 100 $')
    ax.scatter(-100, -100, color = get_color(8), marker = 'v', label = 'Streaming, $ \\beta = 10 $')
    ax.scatter(-100, -100, color = get_color(9), marker = 'v', label = 'Streaming, $ \\beta = 3 $')

  #  plot_fit(ax, creta, rho_contrast, transport_list)

#    X = 10
#    x, y = solve_dcontrast(a = 4./3., X = 20, inv_beta = 0.25)
#    plt.plot(x-0.3, y, color = 'black', label = '$P_c \\propto \\rho^{4/3}$')

    inv_beta = 0.333
#    X = 12.3
    X = 20
    x, y = solve_dcontrast(a = 4./3., X = X, inv_beta = inv_beta)
    plt.plot(x, y, color = 'black', linestyle = 'dashed', label = 'Analytic, $P_c \\propto \\rho^{4/3}$')

    x, y = solve_dcontrast(a = 1./3., X = X, inv_beta = inv_beta)
    plt.plot(x, y, color = 'black', linestyle = 'dashdot', label = 'Analytic, $P_c \\propto \\rho^{1/3}$')

    x, y = solve_dcontrast(a = 1./12., X = X, inv_beta = inv_beta)
    plt.plot(x, y, color = 'black', linestyle = 'dotted', label = 'Analytic, $P_c \\propto \\rho^{1/12}$')

    box = ax.get_position()
    # 1.1, 0.6
#    ax.legend(fontsize = 6, ncol = 1, loc = 'center left', bbox_to_anchor = (1.1, .6),
#              facecolor = 'white', framealpha = 1)
    ax.legend(fontsize = 5, ncol = 1, loc = 'center left', bbox_to_anchor = (0.004, 0.24),#(0.008, 0.28),
              facecolor = 'white', framealpha = 1)
    fig = plt.gcf()
    fig.set_size_inches(4.2, 4)
    #ax.legend(fontsize = 8, ncol = 3, loc = 'upper center', bbox_to_anchor = (0.5, 1.1), facecolor = 'white', framealpha = 1)

    fig.tight_layout()
    fig_basename = 'density_contrast'

    #figname = '../../plots/%s/%s.png'%(sim_fam, fig_basename)
    figname = '../plots/%s.png'%(fig_basename)

    print(figname)
    plt.savefig(figname, dpi = 300, bbox_inches = 'tight', pad_inches = 0.05)

sim_fam = 'production'
plot_density_fluctuation(sim_fam = sim_fam)
