import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

import plotting_tools as pt

sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

# Define and use a simple function to label the plot in axes coordinates
def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)


def create_pdf_plot(x, y, num_cols, bw = 'scott', xlabel = 'x', height = 2, aspect = 4):    
    df = pd.DataFrame(dict(x=x, y=y))

    # Initialize the FacetGrid object
    pal = sns.cubehelix_palette(num_cols, rot=-.25, light=.7)
    g = sns.FacetGrid(df, row="y", hue="y", aspect=aspect, height=height, palette=pal, )

    # Draw the densities in a few steps
    g.map(sns.kdeplot, "x", clip_on=False, shade=True, alpha=1, lw=1.5, bw=bw)
    g.map(sns.kdeplot, "x", clip_on=False, color="w", lw=2, bw=bw)
    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    g.map(label, "x")

    # Set the subplots to overlap
    g.fig.subplots_adjust(hspace=-.25)

    # Remove axes details that don't play well with overlap
  #  g.set(xlabel = xlabel)
    g.set_titles("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)
    return(g)
   
def get_masked_data(field, output, sim_loc, workdir = '../../simulations/production'):
    ds = yt.load('%s/%s/DD%04d/DD%04d'%(workdir, sim_loc, output, output))
    ad = ds.all_data()
    z_abs_code = np.abs(ad[('gas', 'z')] / ds.length_unit.in_units('kpc'))
    z_max = 1.2
    z_min = 0.8
    zmask = (z_abs_code >= z_min) & (z_abs_code <= z_max)

#    cell_mass = ad[('gas', 'cell_mass')][zmask]
#    total_mass = np.sum(cell_mass)
#    weighted_field = ad[('gas', field)][zmask] * cell_mass / total_mass   
    temp = ad[('gas', 'temperature')][zmask]
    print(len(temp[temp < 3.333e5]) / len(temp))
    return np.log10(ad[('gas', field)][zmask].d)

def format_data_for_pdf(field, output, sim_list, label_list, nbins = 100, work_dir = '../../simulations/production'):
    # Create the data
    x = np.array([])
    y = np.array([])
    for i, sim in enumerate(sim_list):
        print(i, label_list[i])
        data, ignore, mass_list = pt.get_2d_hist_data(field, field, sim, work_dir = work_dir, sim_fam = sim_fam)
        hist, bin_edges = np.histogram(data, weights = mass_list, density = False, bins = nbins)
        # normalize hist values:
        max_hist = max(hist)
        hist *= (1000 / max_hist)
        new_array = np.array([])
#        print(max(hist), min(hist), np.median(hist))
        for j in range(len(hist)):
            bin_edge = 0.5*(bin_edges[j] + bin_edges[j+1])
            # add hist[i] copies of the value "bin_edge" to artificially recreate the data of the weighted histogram
            new_array = np.append(new_array, (int(hist[j])*[bin_edge]))
#        if len(data) > 250000:
#            data = signal.resample(data, 250000)
        x = np.append(x, new_array)
        y = np.append(y, len(new_array)*[label_list[i]])
    return x, y



sim_fam = 'production'


sim_list = ['isocool_tctf_0.3_beta_100.0', 'isocool_tctf_0.3_beta_100.0_cr_0.01', \
            'isocool_tctf_0.3_beta_100.0_cr_0.1', 'isocool_tctf_0.3_beta_100.0_cr_1.0',\
            'isocool_tctf_0.3_beta_100.0_cr_10.0']
label_list = ['No CR', 'P_c/P_g = .01', 'P_c/P_g = .1', 'P_c/P_g = 1', 'P_c/P_g = 10']
  

x, y = format_data_for_pdf('temperature', 50, sim_list, label_list)
g = create_pdf_plot(x, y, len(sim_list))  
plt.savefig('test.png', dpi = 300)
