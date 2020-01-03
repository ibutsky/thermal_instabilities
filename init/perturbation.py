"""
Turbulent velocity perturbation generator
-- Drummond Fielding
"""

import h5py
import numpy as np
import numpy.fft
from math import *
from optparse import OptionParser
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt

from constants_and_parameters import *

def init_perturbations(n, kmin, kmax, alpha, dtype, skinny_ratio = 1):
    print('skinny_ratio = %i'%skinny_ratio)
    kx = np.zeros(n, dtype=dtype)
    ky = np.zeros(n, dtype=dtype)
    kz = np.zeros(n, dtype=dtype)
    # perform fft k-ordering convention shifts
    for j in range(0,n[1]):
        for k in range(0,n[2]):
            kx[:,j,k] = n[0]*np.fft.fftfreq(n[0])
    for i in range(0,n[0]):
        for k in range(0,n[2]):
            ky[i,:,k] = n[1]*np.fft.fftfreq(n[1])
    for i in range(0,n[0]):
        for j in range(0,n[1]):
            kz[i,j,:] = n[2]*np.fft.fftfreq(n[2])
            
    kx = np.array(kx, dtype=dtype)
    kx *= skinny_ratio
    ky = np.array(ky, dtype=dtype)
    ky *= skinny_ratio
    kz = np.array(kz, dtype=dtype)

    k = np.sqrt(np.array(kx**2+ky**2+kz**2, dtype=dtype))

    # only use the positive frequencies
    inds = np.where(np.logical_and(k**2 >= kmin**2, k**2 < (kmax+1)**2))
    nr = len(inds[0])

    phasex = np.zeros(n, dtype=dtype)
    phasex[inds] = 2.*pi*np.random.uniform(size=nr)
    fx = np.zeros(n, dtype=dtype)
    fx[inds] = np.random.normal(size=nr)
    
    phasey = np.zeros(n, dtype=dtype)
    phasey[inds] = 2.*pi*np.random.uniform(size=nr)
    fy = np.zeros(n, dtype=dtype)
    fy[inds] = np.random.normal(size=nr)
    
    phasez = np.zeros(n, dtype=dtype)
    phasez[inds] = 2.*pi*np.random.uniform(size=nr)
    fz = np.zeros(n, dtype=dtype)
    fz[inds] = np.random.normal(size=nr)

    # rescale perturbation amplitude so that low number statistics
    # at low k do not throw off the desired power law scaling.
    for i in range(kmin, kmax+1):
        slice_inds = np.where(np.logical_and(k >= i, k < i+1))
        rescale = sqrt(np.sum(np.abs(fx[slice_inds])**2 + np.abs(fy[slice_inds])**2 + np.abs(fz[slice_inds])**2))
        fx[slice_inds] = fx[slice_inds]/rescale
        fy[slice_inds] = fy[slice_inds]/rescale
        fz[slice_inds] = fz[slice_inds]/rescale

    # set the power law behavior
    # wave number bins
    fx[inds] = fx[inds]*k[inds]**-(0.5*alpha)
    fy[inds] = fy[inds]*k[inds]**-(0.5*alpha)
    fz[inds] = fz[inds]*k[inds]**-(0.5*alpha)

    # add in phases
    fx = np.cos(phasex)*fx + 1j*np.sin(phasex)*fx
    fy = np.cos(phasey)*fy + 1j*np.sin(phasey)*fy
    fz = np.cos(phasez)*fz + 1j*np.sin(phasez)*fz

    return fx, fy, fz, kx, ky, kz


def normalize(n, fx, fy, fz):
    norm = np.sqrt(np.sum(fx**2 + fy**2 + fz**2)/np.product(n))
    fx = fx/norm
    fy = fy/norm
    fz = fz/norm
    return fx, fy, fz


def make_perturbations(n, kmin, kmax, alpha, f_solenoidal, dtype=np.float64, skinny_ratio = 1):
    fx, fy, fz, kx, ky, kz = init_perturbations(n, kmin, kmax, alpha, dtype, skinny_ratio = skinny_ratio)
    if f_solenoidal != None:
        k2 = kx**2+ky**2+kz**2
        # solenoidal part
        fxs = 0.; fys =0.; fzs = 0.
        if f_solenoidal != 0.0:
            fxs = np.real(fx - kx*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16))
            fys = np.real(fy - ky*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16))
            fzs = np.real(fz - kz*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16))
            ind = np.where(k2 == 0)
            fxs[ind] = 0.; fys[ind] = 0.; fzs[ind] = 0.
            # need to normalize this before applying relative weighting of solenoidal / compressive components
            norm = np.sqrt(np.sum(fxs**2+fys**2+fzs**2))
            fxs = fxs/norm
            fys = fys/norm
            fzs = fzs/norm
        # compressive part
        # get a different random cube for the compressive part
        # so that we can target the RMS solenoidal fraction,
        # instead of setting a constant solenoidal fraction everywhere.
        fx, fy, fz, kx, ky, kz = init_perturbations(n, kmin, kmax, alpha, dtype, skinny_ratio = skinny_ratio)
        fxc = 0.; fyc =0.; fzc = 0.
        if f_solenoidal != 1.0:
            fxc = np.real(kx*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16))
            fyc = np.real(ky*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16))
            fzc = np.real(kz*(kx*fx+ky*fy+kz*fz)/np.maximum(k2,1e-16))
            ind = np.where(k2 == 0)
            fxc[ind] = 0.; fyc[ind] = 0.; fzc[ind] = 0.
            # need to normalize this before applying relative weighting of solenoidal / compressive components
            norm = np.sqrt(np.sum(fxc**2+fyc**2+fzc**2))
            fxc = fxc/norm
            fyc = fyc/norm
            fzc = fzc/norm
        # back to real space
        pertx = np.real(np.fft.ifftn(f_solenoidal*fxs + (1.-f_solenoidal)*fxc))
        perty = np.real(np.fft.ifftn(f_solenoidal*fys + (1.-f_solenoidal)*fyc))
        pertz = np.real(np.fft.ifftn(f_solenoidal*fzs + (1.-f_solenoidal)*fzc))
    else:
        # just convert to real space
        pertx = np.real(np.fft.ifftn(fx))
        perty = np.real(np.fft.ifftn(fy))
        pertz = np.real(np.fft.ifftn(fz))
    
    # subtract off COM (assuming uniform density)
    pertx = pertx-np.average(pertx)
    perty = perty-np.average(perty)
    pertz = pertz-np.average(pertz)
    # scale RMS of perturbation cube to unity
    pertx, perty, pertz = normalize(n, pertx, perty, pertz)
    return pertx, perty, pertz


def get_erot_ke_ratio(n, pertx, perty, pertz):
    x, y, z = np.mgrid[0:n[0], 0:n[1], 0:n[2]]
    x = x - (n[0]-1)/2.
    y = y - (n[1]-1)/2.
    z = z - (n[2]-1)/2.
    r2 = x**2+y**2+z**2
    erot_ke_ratio = (np.sum(y*pertz-z*perty)**2 +
                     np.sum(z*pertx-x*pertz)**2 +
                     np.sum(x*perty-y*pertx)**2)/(np.sum(r2)*np.product(n)) 
    return erot_ke_ratio


def plot_spectrum1D(n, pertx, perty, pertz, kmin, kmax, dtype):
    # plot the 1D power to check the scaling.
    fx = np.abs(np.fft.fftn(pertx))
    fy = np.abs(np.fft.fftn(perty))
    fz = np.abs(np.fft.fftn(pertz))
    fx = np.abs(fx)
    fy = np.abs(fy)
    fz = np.abs(fz)
    kx = np.zeros(n, dtype=dtype)
    ky = np.zeros(n, dtype=dtype)
    kz = np.zeros(n, dtype=dtype)
    # perform fft k-ordering convention shifts
    for j in range(0,n[1]):
        for k in range(0,n[2]):
            kx[:,j,k] = n[0]*np.fft.fftfreq(n[0])
    for i in range(0,n[0]):
        for k in range(0,n[2]):
            ky[i,:,k] = n[1]*np.fft.fftfreq(n[1])
    for i in range(0,n[0]):
        for j in range(0,n[1]):
            kz[i,j,:] = n[2]*np.fft.fftfreq(n[2])
    k = np.sqrt(np.array(kx**2+ky**2+kz**2,dtype=dtype))
    k1d = []
    power = []
    for i in range(kmin,kmax+1):
        slice_inds = np.where(np.logical_and(k >= i, k < i+1))
        k1d.append(i+0.5)
        power.append(np.sum(fx[slice_inds]**2 + fy[slice_inds]**2 + fz[slice_inds]**2))
    plt.loglog(k1d, power)
    plt.savefig("spectrum1D.png", bbox_inches='tight', dpi=200)
    plt.clf()

def plot_slice(pertx,name):
    plt.pcolormesh(pertx[:,:,0])
    ax = plt.gca()
    ax.set_aspect(1.0)
    plt.savefig("pert_slice_%s.png" %name, bbox_inches='tight', dpi=200)
    plt.clf()


def read_parameters_from_command_line():
    parser = OptionParser()
    parser.add_option('--kmin', dest='kmin', 
                  help='minimum wavenumber.', 
                  default=default_kmin)
    parser.add_option('--kmax', dest='kmax', 
                  help='maximum wavenumber.', 
                  default=default_kmax)
    parser.add_option('--size', dest='size', 
                  help='size of each direction of data cube.  default=256', 
                  default=default_n)
    parser.add_option('--size_x', dest='size_x',
                  help='size to use in x-direction of data cube.  default=256',
                  default=0)
    parser.add_option('--size_y', dest='size_y',
                  help='size to use in y-direction of data cube.  default=256',
                  default=0)
    parser.add_option('--size_z', dest='size_z',
                  help='size to use in z-direction of data cube.  default=256',
                  default=0)
    parser.add_option('--alpha', dest='alpha', 
                  help='negative of power law slope.  (Power ~ k^-alpha) '+
                  'supersonic turbulence is near alpha=2.  '+
                  'driving over a narrow band of two modes is often done with alpha=0', 
                  default = default_alpha)
    parser.add_option('--seed', dest='seed', 
                  help='seed for random # generation.  default=0', 
                  default = default_seed)
    parser.add_option('--f_solenoidal', dest='f_solenoidal', 
                  help='volume RMS fraction of solenoidal componet of the perturbations relative to the total.  ' + 
                       'If --f_solenoidal=None, the motions are purely random.  For low wave numbers ' +
                       'the relative imporance of solenoidal to compressive may be sensitve to the ' +
                       'choice of radom seed.  It has been suggested (Federrath 2008) that ' +
                       'f_solenoidal=2/3 is the most natural driving mode and this is currently' + 
                       'the suggested best-practice.', 
                  default = default_f_solenoidal)

    (options, args) = parser.parse_args()

    # size of the data domain
    if int(options.size_x) > 0:
        n = [int(options.size_x), int(options.size_y), int(options.size_z)]
    else:
        n = [int(options.size), int(options.size), int(options.size)]

    # range of perturbation length scale in units of the smallest side of the domain
    kmin = int(options.kmin)
    kmax = int(options.kmax)
    if kmin > kmax or kmin < 0 or kmax < 0:
        print(kmin, kmax)
        print("kmin must be < kmax, with kmin > 0, kmax > 0.  See --help.")
        sys.exit(0)
#    if kmax > floor(np.min(n))/2:
#        print("kmax must be <= floor(size/2).  See --help.")
#        sys.exit(0)
    f_solenoidal = options.f_solenoidal
    if f_solenoidal == "None" or f_solenoidal == "none":
        f_solenoidal = None
    else:
        f_solenoidal = float(options.f_solenoidal)
        if f_solenoidal > 1. or f_solenoidal < 0.:
            print("You must choose f_solenoidal.  See --help.")
            sys.exit(0)
    alpha = options.alpha
    if alpha==None:
        print("You must choose a power law slope, alpha.  See --help.")
        sys.exit(0)
    alpha = float(options.alpha)
    if alpha < 0.:
        print("alpha is less than zero. Thats probably not what you want.  See --help.")
        sys.exit(0)
    seed = int(options.seed)
    # data precision
    dtype = np.float64
    # ratio of solenoidal to compressive components
    if options.f_solenoidal=="None" or options.f_solenoidal==None:
        f_solenoidal = None
    else:
        f_solenoidal = min(max(float(options.f_solenoidal), 0.), 1.)
    return n, kmin, kmax, alpha, f_solenoidal, seed, dtype



def generate_perturbation_infile(skinny_ratio = 1):
    n, kmin, kmax, alpha, f_solenoidal, seed, dtype = read_parameters_from_command_line()
    n = [32, 32, 128]
    kmin = 4
    kmax = 64

    np.random.seed(seed=seed)      
    pertx, perty, pertz = make_perturbations(n, kmin, kmax, alpha, f_solenoidal, skinny_ratio = skinny_ratio)
#    erot_ke_ratio = get_erot_ke_ratio(n, pertx, perty, pertz)
#    print("erot_ke_ratio = ", erot_ke_ratio)
    plot_spectrum1D(n, pertx, perty, pertz, kmin, kmax, dtype = dtype)
    plot_slice(pertx, "kmin%i_kmax%i_alpha%i"%(kmin, kmax, alpha))


    # Write the perturbations to a file.
    outf = open('perturbation.in', 'w')

    # Write out a header.
    outf.write('# density perturbation field\n')
    outf.write('# k_min: %d\n' % kmin)
    outf.write('# k_max: %d\n' % kmax)
    outf.write('# alpha: %d\n' % alpha)
    outf.write('# seed: %d\n'  % seed)
    outf.write('# f_solenoidal: %d\n' % f_solenoidal)
    outf.write('# Dimensions: %i\n' % len(n))
    outf.write('# Grid size:')

    ncopies = [32, 32, 128]
    for component in ncopies:
        outf.write(' %i' % component)
    outf.write('\n')

    # First non commented line is the dimensions
    outf.write('%i\n' % len(n))

    # Second non commented line is the grid size
    for component in ncopies:
        outf.write('%i ' % component)
    outf.write('\n')

    # Third non commented line is the standard deviation
    pertx_1d = pertx.flatten()
    std = np.std(pertx_1d)
    outf.write('%f\n' %std)

    # Write out each index
    for i in range(n[0]):
        for j in range(n[1]):
            for k in range(n[2]):
                outf.write('%i %i %i %e\n'%(i, j, k, pertx[i,j,k]))
               
    outf.close()


generate_perturbation_infile()
