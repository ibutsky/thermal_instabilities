import yt
from yt import YTArray
from yt import YTQuantity

import numpy as np
import sys
import os
import glob

import matplotlib.pylab as plt
from matplotlib.colors import SymLogNorm, LogNorm

import multiprocessing as mp

import palettable
import plotting_tools as pt
import yt_functions as ytf


work_dir = '../../simulations/'
plot_folder = '../../movies/temp'

def fft_comp(ds, field, level = 0 ):
    cube = ds.covering_grid(level, left_edge=ds.domain_left_edge,
                            dims=ds.domain_dimensions)

    rho = cube[('gas', 'density')].d
    if field == 'rho':
        fft_field = rho
    elif field == 'drho':
        drho = np.ndarray(shape = rho.shape)
        for i in range(len(rho)):
            rho_slice = rho[:, :, i]
            rho_ave = np.mean(rho_slice)
            drho[:, :, i]  = (rho_slice - rho_ave) / rho_ave
        fft_field = drho

    nx, ny, nz = rho.shape

    # do the FFTs -- note that since our data is real, there will be     
    # too much information here.  fftn puts the positive freq terms in
    # the first half of the axes -- that's what we keep.  Our 
    # normalization has an '8' to account for this clipping to one
    # octant. 

    ru = np.fft.fftn(fft_field)[0:nx//2+1,0:ny//2+1,0:nz//2+1]
    ru = 8.0*ru/(nx*ny*nz)

    return np.abs(ru)**2

def make_power_spectrum(ds, field = 'drho'):
    dims = ds.domain_dimensions
    nx, ny, nz = dims

    Kk = np.zeros( (nx//2+1, ny//2+1, nz//2+1))
    Kk = fft_comp(ds, field)

    # wavenumbers in units of box length
    L = np.array([1.0, 1.0, 1.0])
    
    kx = np.fft.rfftfreq(nx)*nx/L[0]
    ky = np.fft.rfftfreq(ny)*ny/L[1]
    kz = np.fft.rfftfreq(nz)*nz/L[2]

    # physical limits to the wavenumbers
    kmin = np.min(1.0/L)
    kmax = np.min(0.5*np.array(dims)/L)
    
    
    kbins = np.arange(kmin, kmax, kmin)
    N = len(kbins)

    # bin the Fourier KE into radial kbins
    kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")
    k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

    whichbin = np.digitize(k.flat, kbins)
    ncount = np.bincount(whichbin)

    P_spectrum = np.zeros(len(ncount)-1)

    for n in range(1,len(ncount)):
        P_spectrum[n-1] = np.sum(Kk.flat[whichbin==n])

    k = kbins[1:N]
    P_spectrum = P_spectrum[1:N]
    plt.loglog(k, P_spectrum)
    plt.show()
    return k, P_spectrum


def plot_power_spectrum(output, field = 'drho'):
    fig, ax = plt.subplots(figsize = (6, 6))
    for tctf, beta, cr, label, color, k in zip(tctf_list, beta_list, cr_list, label_list, color_list, k_list):
        # load the simulation
        sim_loc = '%s/%s_tctf_%.1f'%(work_dir, model, tctf)
        if beta != 'inf':
            sim_loc += '_beta_%.1f'%(beta)
            if beta == 256:
                sim_loc += '_k_%.1f'%(k)
        if cr > 0:
            sim_loc += '_cr_%.1f'%cr
        ds_path = '%s/DD%04d/DD%04d'%(sim_loc, output, output)
        if not os.path.isfile(ds_path):
            return
        ds = yt.load(ds_path)
        
        k, P_k = make_power_spectrum(ds, field)
        ax.loglog(k, P_k, label = label, linewidth = 3, color = color)
     
    ax.set_xlabel('k')
    ax.set_ylabel('P(k)')
    ax.set_xlim(1, 1e2)
    ax.set_ylim(1e-6, 1e2)
    ax.legend(loc = 1)
    fig.tight_layout()
    figname = '%s/%04d.png'%(plot_folder, output)
    plt.savefig(figname, dpi = 300)

    

model = sys.argv[1]
compare = sys.argv[2]
#tctf = float(sys.argv[3])
#beta = 100.0 


if compare == 'tctf':
    color_list  = palettable.scientific.sequential.Batlow_6.mpl_colors
    tctf_list = [0.1, 0.3, 1.0, 3.0, 10.0]
    beta_list = 5*[beta]
    cr_list = 5*[0]
    label_list = []
    for tctf in tctf_list:
        label_list.append("t$_{cool}$/t$_{ff}$ = %.1f"%tctf)

elif compare == 'beta':
    color_list  = palettable.scientific.sequential.Batlow_6.mpl_colors
    tctf_list = 3*[tctf]
    beta_list = [4, 100, 'inf']
    cr_list = 3*[0]
    label_list = []
    for beta in beta_list:
        if beta == 'inf':
            label_list.append('Hydro')
        else:
            label_list.append("$\\beta =$ %.1f"%beta)

elif compare == 'ktest':
    color_list = palettable.scientific.sequential.Batlow_6.mpl_colors
    tctf_list  = [3, 3]
    beta_list  = [256, 256]
    k_list = [ 4, 32]
    cr_list = [0, 0]
    label_list = ['k = 4', 'k = 32']
    

elif compare == 'cr':
    color_list  = palettable.scientific.sequential.Batlow_6.mpl_colors
    tctf_list = 6*[1.0]
    beta_list = 6*[10.0]
    cr_list = [0, 0.1, 0.3, 1.0, 3.0, 10.0]
    k_list = 6*[0]
    label_list = []
    for cr in cr_list:
        label_list.append('$\\eta = $%0.1f'%cr)

    


pool = mp.Pool(mp.cpu_count())
print("Number of processors: ", mp.cpu_count())

output_list = np.arange(0, 101, 1)
pool.map(plot_power_spectrum, [output for output in output_list])
pool.close()

cwd = os.getcwd()
os.chdir(plot_folder)
os.system('ffmpeg -r 10 -f image2 -s 1920x1080 -i %04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p power_spectrum.mov')
if compare == 'tctf':
    if beta == 'inf':
        os.rename('power_spectrum.mov', '../power_spectrum_%s.mov'%(model))
    else:
        os.rename('power_spectrum.mov', '../power_spectrum_%s_beta_%.1f.mov'%(model, beta))
elif compare == 'beta':
    os.rename('power_spectrum.mov', '../power_spectrum_%s_tctf_%.1f.mov'%(model, tctf))
elif compare == 'ktest':
    os.rename('power_spectrum.mov', '../power_spectrum_%s_ktest.mov'%model)
elif compare == 'cr':
    os.rename('power_spectrum.mov', '../power_spectrum_%s_cr.mov'%model)

png_files = glob.glob('*.png')
for pic in png_files:
    os.remove(pic)
os.chdir(cwd)

