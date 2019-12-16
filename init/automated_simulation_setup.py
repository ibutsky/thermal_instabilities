import os
import shutil

sim_dir = '../../simulations'

def get_folder_name(tctf, halo_prof, perturb, beta, cr, grid_rank = 3, resolution = 128):
    if halo_prof == 1:
        sim = 'isothermal'
    elif halo_prof == 2:
        sim = 'isentropic'
    elif halo_prof == 3:
        sim = 'isocool'
    if perturb == 2:
        sim += '_isochoric'

    folder_name = sim
    folder_name += '_tctf_%.1f'%(tctf)

    if beta != 'inf':
        folder_name += '_beta_%.1f'%beta

    if cr > 0:
        if cr >= 0.1:
            folder_name += '_cr_%.1f'%cr
        elif cr < 0.01:
            folder_name += '_cr_%.3f'%cr
        else:
            folder_name += '_cr_%.2f'%cr

        if cr_diffusion:
            folder_name += '_diff_%.1f'%(cr_kappa/1e29)

        if cr_streaming:
            folder_name += '_stream'

        if cr_heating:
            folder_name += '_heat'
    if grid_rank == 2:
        folder_name += '_2d_%i'%resolution
    return folder_name
    


def create_constants_and_parameters_file(fn, tcool_tff_ratio, halo_profile, perturb_type, beta, eta, \
                                         grid_rank = 3, ndim = 128, wall_time = '16:00:00'):
    f = open(fn, 'w')
    f.write('from yt import YTQuantity\n')
    f.write('from yt import YTArray\n\n')

    f.write('tcool_tff_ratio = %f    # ratio of t_cool / t_free_fall\n'%tcool_tff_ratio)
    f.write('halo_profile    = %i    # 1 = isothermal, 2 = isentropic, 3 = iso-tcool\n'%halo_profile)
    f.write('perturb_type    = %i \n'%perturb_type)
    
    if beta == 'inf':
        f.write('magnetic_pressure_ratio = 0.0 \n')
    else:
        f.write('magnetic_pressure_ratio = %f \n'%(1.0 / beta))
    f.write('cr_pressure_ratio       = %f\n\n' %eta)

    f.write('cr_diffusion = %i\n'%cr_diffusion)
    f.write('cr_kappa     = %e\n'%cr_kappa)
    f.write('cr_streaming = %i\n'%cr_streaming)
    f.write('cr_heating = %i\n'%cr_heating)
    
    
    f.write('########### These parameters don\'t really need to change ###########\n')
    f.write('bfield_direction    = [1, 0, 0]\n')
    f.write('resolution          = %i  # resolution along z-axis; scaled for x-y axes\n'%ndim)
    f.write('grid_rank           = %i\n'%grid_rank)

    f.write('# cooling function parameters\n')
    f.write('# Tmin should be ~ T0/20    \n')
    f.write('T_min = YTQuantity(5e4, \'K\')\n')
    f.write('T_max = YTQuantity(1e9, \'K\')\n')
    f.write('T_power_law_index = (-2./3.)\n')
    f.write('smooth_factor = 0.02\n\n')

    f.write('####### gas parameters ######\n')
    f.write('# changing these will just rescale the units\n')
    f.write('T0              = YTQuantity(1e6,   \'K\')    \n')
    f.write('rho0            = YTQuantity(1e-27, \'g/cm**3\') \n')
    f.write('g0              = YTQuantity(5e-10, \'cm/s**2\')\n')
    f.write('gsoft_scale     = 0.1    # units of scale height\n\n')

    f.write('###### perturbatio default parameters\n')
    f.write('perturbation_amplitude = 0.02\n')
    f.write('default_n              = resolution\n')
    f.write('default_kmin           = 4\n')
    f.write('default_kmax           = int(resolution / 2)\n')
    f.write('default_f_solenoidal   = 2./3.\n')
    f.write('default_alpha          = 0\n')
    f.write('default_seed           = 4085281318\n')

    f.write('###### computing parameters\n')
    f.write('nodes              = %i\n'%nodes)
    f.write('tasks_per_node     = 32\n')
    f.write('num_cooling_cycles = 10  # simulation stop time = 10 * t_cool\n')
    f.write('num_outputs        = 100\n')
    f.write('wall_time          = \'%s\'\n'%(wall_time))

    f.close()


grid_rank = 2
ndim = 256
nodes = 2

num = 5
#tctf_list = num*[1.0]
tctf_list = [0.1, 0.3, 1.0, 3.0, 10.0]
halo_prof_list = num*[3]
perturb_list = num*[1]
beta_list = num*[10]
#beta_list = ['inf', 300, 100, 30, 10, 3]
#cr_list = [0.01, 0.1, 1.0, 10.0, 100.0]
cr_list = num*[0]
wall_time = '0:20:00'

cr_diffusion = 0
cr_kappa     = 0 #3e29 #1e28
cr_streaming = 0
cr_heating   = 0

for i in range(len(tctf_list)):
    folder_basename = get_folder_name(tctf_list[i], halo_prof_list[i], perturb_list[i], \
                                      beta_list[i], cr_list[i], grid_rank = grid_rank, resolution = ndim)
    folder_path = '%s/%s'%(sim_dir, folder_basename)

    if not os.path.isdir(folder_path):
        os.mkdir(folder_path)

    else:
        print("%s already exists; skipping\n"%folder_basename)
        continue

    cwd = os.getcwd()
    os.chdir(folder_path)
    shutil.copyfile('%s/perturbation_%i.in'%(cwd, ndim), 'perturbation.in')
    shutil.copyfile('%s/generate_initial_conditions.py'%cwd, 'generate_initial_conditions.py')

    create_constants_and_parameters_file('constants_and_parameters.py', tctf_list[i],\
                                         halo_prof_list[i], perturb_list[i], beta_list[i], cr_list[i], \
                                         grid_rank = grid_rank, ndim = ndim, wall_time = wall_time)

    os.system('python generate_initial_conditions.py')
    os.system('sbatch submit.sbatch')
    os.chdir(cwd)
