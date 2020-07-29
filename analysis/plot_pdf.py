import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import astropy.constants as const
import os
import palettable
import plotting_tools as pt

sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

mu = 1.22
mh = const.m_p.cgs.value
log_mumh = np.log10(mu*mh)

# Define and use a simple function to label the plot in axes coordinates
def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)


def create_pdf_plot(x, y, num_cols, bw = 'scott', xlabel = 'x', height = 2, aspect = 7, pal = None, use_label = True):    
    df = pd.DataFrame(dict(x=x, y=y))

    # Initialize the FacetGrid object
    if pal is None:
        pal = sns.cubehelix_palette(num_cols, rot=-.25, light=.7)
    g = sns.FacetGrid(df, row="y", hue="y", aspect=aspect, height=height, palette=pal, )

    # Draw the densities in a few steps
    g.map(sns.kdeplot, "x", clip_on=False, shade=True, alpha=1, lw=1.5, bw=bw)
    g.map(sns.kdeplot, "x", clip_on=False, color="w", lw=2, bw=bw)
    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    if use_label:
        g.map(label, "x")

    # Set the subplots to overlap
    g.fig.subplots_adjust(hspace=-.25)

    # Remove axes details that don't play well with overlap
  #  g.set(xlabel = xlabel)
    g.set_titles("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)
    return(g)
   
def format_data_for_pdf(field, sim_list, label_list, nbins = 100, weighted = True, work_dir = '../../simulations', sim_fam = 'production'):
    # Create the data
    x = np.array([])
    y = np.array([])
    for i, sim in enumerate(sim_list):
        print(i, label_list[i])
        data, ignore, mass_list = pt.get_2d_hist_data(field, field, sim, work_dir = work_dir, sim_fam = sim_fam)
        print(len(data))
        if weighted:
            weights = mass_list
        else:
            weights = np.array(len(mass_list)*[1.0])
        if field == 'density':
            data -= log_mumh
        hist, bin_edges = np.histogram(data, weights = weights, density = True, bins = nbins)
        # normalize hist values:
        print("hist: %e"%np.sum(hist*np.diff(bin_edges)))
        max_hist = max(hist)
        if max_hist > 0:
            hist *= (500 / max_hist)
        new_array = np.array([])
#        print(max(hist), min(hist), np.median(hist))
        for j in range(len(hist)):
            bin_edge = 0.5*(bin_edges[j] + bin_edges[j+1])
            # add hist[i] copies of the value "bin_edge" to artificially recreate the data of the weighted histogram
            new_array = np.append(new_array, (int(hist[j])*[bin_edge]))
        x = np.append(x, new_array)
        print("new array: %e \n"%np.sum(new_array))
        y = np.append(y, len(new_array)*[label_list[i]])
    return x, y

def make_plot(field, compare, tctf = 0.3, cr = 1, weighted = True, nbins = 100, 
              work_dir = '../../simulations', plot_dir = '../../plots', sim_fam = 'production'):

    sim_list = pt.generate_sim_list(compare, tctf = tctf, cr = cr)
    label_list = pt.generate_label_list(compare, tctf = tctf, cr = cr)

    if field == 'cr_eta':
        sim_list = sim_list[1:]
        label_list = label_list[1:]
    #color_list = pt.get_color_list(compare)
    if field == 'density':
        pal = sns.cubehelix_palette(len(sim_list), rot=-.25, light=.7)
    elif field == 'temperature':
        pal = palettable.scientific.sequential.LaJolla_13.mpl_colors[2:-1]
        if compare == 'transport_pdf':
            pal = palettable.scientific.sequential.LaJolla_16.mpl_colors[2:-1]

    elif field == 'cr_eta':
        pal = sns.cubehelix_palette(len(sim_list)+1)[1:]


    x, y = format_data_for_pdf(field, sim_list, label_list, weighted = weighted, 
                               nbins = nbins, work_dir = work_dir, sim_fam = sim_fam)

    if field == 'density':
        xlabel = 'Log Number Density (cm$^{-3}$)'
#        xlims = (-28.5, -25.8)
        xlims = (-28.7, -25.6)
        xlims -= log_mumh
        aspect = 8
        use_label = True
        ylims = (0, 1)
    elif field == 'temperature':
        xlabel = 'Log Temperature (K)'
 #       xlims = (3.5, 7)
        xlims = (4.3, 6.8)
        aspect = 8.137 #7.5
        use_label = False
        ylims = (0, 1)
    elif field == 'cr_eta':
        xlabel = 'Log ($P_c / P_g$)'
  #      xlims = (np.log10(cr) - 2, np.log10(cr) + 2)
        xlims = (np.log10(cr) - 0.8, np.log10(cr) + 1.7)
        aspect = 8
        use_label = False
        ylims = (0, 1)
    ylims = (0, 1)

    g = create_pdf_plot(x, y, len(sim_list), height = 1, aspect = aspect, pal = pal, use_label = use_label)
    ax = g.axes

    ax[-1][0].set_xlabel(xlabel, color = 'black', fontsize = 16)
    ax[-1][0].tick_params(axis='x', colors='black', bottom = True, labelsize = 'large')
    g.set(xlim = xlims)
    g.set(ylim = ylims)
    if field == 'cr_eta':
        field = 'creta'
    fig_basename = 'pdf_%s'%field
    if weighted:
        fig_basename += '_weighted'
    figname = pt.get_fig_name(fig_basename, 'isocool', compare, \
                              tctf, 100, cr, sim_fam = sim_fam,\
                              loc = plot_dir)
    g.savefig(figname, dpi = 300, bbox_inches = 'tight', pad_inches = 0)




sim_fam = 'production'#/Tmin1e4'

cr = 1
tctf = 0.3
compare = 'transport_pdf'
#compare = 'cr'
cr_list = [0.01, 0.1, 1, 10]
tctf_list = [0.1, 0.3, 1, 3]

cr_list = [0.1, 10] 
tctf_list = [0.3, 1.0]
for cr in cr_list:
    for tctf in tctf_list:
        for field in ['temperature', 'cr_eta']:#, 'density', 'temperature']:
            make_plot(field, compare, tctf, cr, weighted = True, sim_fam = sim_fam)

