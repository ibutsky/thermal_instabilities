import yt
from yt import YTQuantity
from yt.data_objects.level_sets.api import *


ds = yt.load('simulations/transfer/isocool_tctf_0.3_beta_100.0/DD0050/DD0050') 
ad = ds.all_data()

field = 'density'

master_clump = Clump(ad, ('gas', field))

master_clump.add_validator("min_cells", 20)
master_clump.add_info_item("total_cells")


c_min = YTQuantity(1e-27, 'g/cm**3')#ad["gas", field].min()
c_max = ad["gas", field].max()
step = 100.0
find_clumps(master_clump, c_min, c_max, step)

leaf_clumps = master_clump.leaves

p = yt.ProjectionPlot(ds, 'x', ('gas', 'temperature'), weight_field = 'density')   
p.annotate_clumps(leaf_clumps) 
p.save('temp_clumps.png') 
