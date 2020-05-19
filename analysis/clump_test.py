import yt
from yt import YTQuantity
from yt.data_objects.level_sets.api import *
import numpy as np
import sys

def _negative_temp(field, data):
    return -data[('gas', 'temperature')]

sim_name = 'isocool_tctf_0.1_beta_100.0'
sim_name = sys.argv[1]
ds = yt.load('../../simulations/production/%s/DD0040/DD0040'%sim_name) 
ds.add_field(('gas', 'negative_temperature'), function = _negative_temp, sampling_type = 'cell', 
             units = 'K')
ad = ds.all_data()
ad_slice = ds.r[:, :, 0.8:1.2]
cold_mask = ad_slice[('gas', 'temperature')] < 3.33e5
rho_cold = ad_slice[('gas', 'density')][cold_mask]
rho_cold_ave = np.mean(rho_cold)
rho_cold_std = np.std(rho_cold)

field = 'density'
master_clump = Clump(ad, ('gas', field))

master_clump.add_validator("min_cells", 20)
# notes for tomorrow:
# try finding clumps using inverse temperature or negative!
# maybe try adding validator for temperature or density cutoffs
# idea: contour around 1 sigma of average cold temperature
# idea: filter out leafs (especially larger ones... what's going on there...
#        maybe i should change step back down to something  small t oallow smaller clumps
# quantify  shape of clump
#master_clump.add_validator("max_temperature")
#idea
master_clump.add_info_item("total_cells")


#c_min = YTQuantity(1e-27, 'g/cm**3')#
#c_max = ad["gas", field].max()
#c_max = YTQuantity(-4e4, 'K')
#c_min = YTQuantity(-3e5, 'K')#ad["gas", field].max()
c_min = rho_cold_ave - 0.5*rho_cold_std
c_max = ad["gas", field].max()  
print(c_min, c_max)

step = 1000.0
find_clumps(master_clump, c_min, c_max, step)

leaf_clumps = master_clump.leaves

p = yt.ProjectionPlot(ds, 'x', ('gas', 'temperature'), weight_field = 'density')   
#p = yt.SlicePlot(ds, 'x', ('gas', 'temperature'))
p.set_zlim(('gas', 'temperature'), 1e5, 3e6)
p.annotate_clumps(leaf_clumps) 
p.save('../../plots/clump_plots/%s.png'%sim_name) 

radius_list = []
for i, clump in enumerate(leaf_clumps):
    z = np.abs(clump[('gas', 'z')].in_units('kpc'))
    danger_zone = (z < 5) | (z > 83)
    if len(z[danger_zone]) > 0:
        print(i, "danger zone", len(z[danger_zone]) / len(z))
        continue
    else:
        ncells = len(z)
        radius = np.power(ncells, 1./3.)
        print(i, ncells, radius)
        radius_list.append(radius)

print(np.mean(radius_list), np.median(radius_list), np.std(radius_list))
