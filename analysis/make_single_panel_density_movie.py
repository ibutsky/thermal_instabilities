import yt

import numpy as np
import glob
import os
import sys

import matplotlib.pylab as plt
from matplotlib.colors import LogNorm, SymLogNorm
import palettable
import astropy.constants as const

import multiprocessing as mp

import plotting_tools as pt

half_range = 1
rho0 = 1e-27
T0 = 1e6
mu = 1.22
mh = const.m_p.cgs.value
kb = const.k_B.cgs.value
p0 = (rho0 / mu / mh) * kb*T0

workdir = '../../simulations/2d_4096'
plot_folder = '../../movies/2d_4096/temp'

sim = sys.argv[1]
field = sys.argv[2]
orient = sys.argv[3]
last_output = int(sys.argv[4])

def plot_density_slices(output, folder = '.'):
    cmap_list = [palettable.cmocean.sequential.Tempo_20.mpl_colormap]
    cmap_list = [palettable.scientific.sequential.Oslo_20.mpl_colormap]

    if field == 'temperature':
        cmap_list = [palettable.scientific.sequential.LaJolla_20_r.mpl_colormap]

    fig, ax = plt.subplots(ncols = 1, nrows = 1, figsize=(5, 5))
    print('%s/%s/DD%04d/DD%04d'%(workdir, sim, output, output))
    ds = yt.load('%s/%s/DD%04d/DD%04d'%(workdir, sim, output, output))

    s = yt.SlicePlot(ds, orient, ('gas', field))
    s.set_cmap(('gas', field), cmap_list[0])
    if field == 'density':
        s.set_zlim(('gas', field), 1e-28, 1e-26)
    elif field == 'temperature':
        s.set_zlim(('gas', field), 5e4, 1e7)
    s.hide_axes()
    s.hide_colorbar()
    s.set_buff_size(4096)
    figname = '%s/%04d.png'%(plot_folder, output)
    s.save(figname, mpl_kwargs={'dpi':300})


output_list = np.arange(0, 201, 1) 
output_list = np.arange(0, last_output, 1)

#output_list = np.arange(30, 50, 1)
pool = mp.Pool(mp.cpu_count())
print("Number of processors: ", mp.cpu_count())

pool.map(plot_density_slices, [output for output in output_list])
pool.close()


cwd = os.getcwd()
os.chdir(plot_folder)
#os.system('ffmpeg -r 10 -f image2 -s 1920x1080 -i %04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p density.mov')
os.system('ffmpeg -framerate 12 -pattern_type glob -i "*.png" -c:v mpeg4 -pix_fmt yuv420p -q:v 3 %s.mov'%field)
os.rename('%s.mov'%field, '../%s_%s_only.mov'%(sim, field))
png_files = glob.glob('*.png')
#for pic in png_files:
#    os.remove(pic)
os.chdir(cwd)

