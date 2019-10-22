import yt
from yt import YTArray
from yt import YTQuantity
import sys
import numpy as np
import matplotlib.pylab as plt
import palettable

import plotting_tools as pt


def plot_cold_fraction(model, beta_list = None, tctf_list = None, output_list = [50], \
                                    work_dir = '../../simulations', Tmin = 3.3333333e5):
    if tctf_list == None:
        tctf_list = np.array([0.1, 0.3, 1.0, 10.0])

    fig, ax = plt.subplots(figsize = (6, 6))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.08, 15)
    ax.set_ylim(1e-2, 1)
    
    cold_fraction_list = np.array([])
    for tctf in tctf_list:        
        sim_location = '%s/%s_tctf_%.1f'%(work_dir, model, tctf)
        cold_fraction_list = np.append(cold_fraction_list, \
                          pt.calculate_averaged_cold_fraction(output_list, sim_location, Tmin = Tmin))

    # add in a couple measurements by hand
    sim_location = '%s/%s_64'%(work_dir, model)
    cold_fraction = pt.calculate_averaged_cold_fraction(output_list, sim_location, Tmin = Tmin)
    ax.scatter(1, cold_fraction, marker = 'v', color = 'black', label = 'res = 64$^3$')

    sim_location = '%s/%s_256'%(work_dir, model)
    cold_fraction = pt.calculate_averaged_cold_fraction(output_list, sim_location, Tmin = Tmin)
    ax.scatter(1, cold_fraction, marker = '^', color = 'black', label= 'res = 256$^3$')

    mask = cold_fraction_list > 0
    ax.plot(tctf_list[mask], cold_fraction_list[mask], color = 'black', label = 'Hydro', marker = 'o')

    ax.set_xlabel('t$_{cool}$ / t$_{ff}$')
    ax.set_ylabel('Cold Fraction')
    ax.legend()
    fig.tight_layout()
    plt.savefig('../../plots/cold_fraction_%s.png'%model, dpi = 300)


tctf_list = [0.1, 0.3, 1.0, 10]
                                                                
Tmin = 1e6 / 3.0
model = sys.argv[1]
output_list = [90, 95, 100]
plot_cold_fraction(model, output_list = output_list, Tmin = Tmin)
