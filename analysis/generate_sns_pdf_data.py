import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
%matplotlib inline
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
   

def generate_sns_histogram(field, sim_loc, output_list, sim_loc, weighted = True, use_log = True, \
                    nbins = 100, cold_mask = False, z_max = 1.2, z_min = 0.8):
    combined_dens_array = np.array([])
    combined_temp_array = np.array([])
    for output in output_list:
        ds = ytf.load('%s/%s/DD%04d/DD%04d'%(workdir, sim_loc, output, output))
        ad = ds.all_data()
        z_abs_code = np.abs(ad[('gas', 'z')] / ds.length_unit.in_units('kpc'))
        zmask = (z_abs_code >= z_min) & (z_abs_code <= z_max)

        if cold_mask:
            zmask = zmask & (ad[('gas', 'temperature')] <= 3.333333e5)
        
        dens = ad[('gas', 'density')][zmask].d
        temp = ad[('gas', 'temperature')][zmask].d

        if weighted:
            weights = np.array(ad[('gas', 'cell_mass')][zmask].in_units('Msun') / 1e4)
        else:
            weights = np.ones(len(dens))
    
        if use_log:
            dens = np.log10(dens)
            temp = np.log10(temp)
        combined_dens_array = np.append(combined_dens_array, dens)
        combined_temp_array = np.append(combined_temp_array, temp)

    dens_hist, dens_bin_edges = np.histogram(combined_dens_array, weights = weights, density = False, bins = nbins)
    temp_hist, temp_bin_edges = np.histogram(combined_temp_array, weights = weights, density = False, bins = nbins)
    
    new_dens_array = np.array([])
    new_temp_array = np.array([])
    for i in range(len(dens_hist)):
        dens_bin_edge = 0.5*(dens_bin_edges[i] + dens_bin_edges[i+1])
        new_dens_array = np.append(new_dens_array, (int(dens_hist[i]) * [dens_bin_edge]))
        temp_bin_edge = 0.5*(temp_bin_edges[i] + temp_bin_edges[i+1])
        new_temp_array = np.append(new_temp_array, (int(temp_hist[i]) * [temp_bin_edge]))
    # save thisi to an h5 File
    return new_array


def format_data_for_pdf(field, output, sim_list, label_list, nbins = 100, use_log = True, weighted = False, \
                            cold_mask = False, workdir = '../../simulations/production'):
    # Create the data
    x = []
    y = []
    for i, sim_loc in enumerate(sim_list):
        data = get_masked_data(field, output, sim_loc, nbins = nbins, weighted = weighted, use_log = use_log, 
                                   cold_mask = cold_mask, workdir = workdir)
        x = np.append(x, data)
        y = np.append(y, len(data)*[label_list[i]])
    return x, y
