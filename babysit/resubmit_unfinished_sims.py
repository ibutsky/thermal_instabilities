import os
import glob
import numpy as np

def generate_restart_sbatch_file(sim_loc, output, nodes = 1, tasks_per_node = 32, wall_time = '4:00:00'):
    cwd = os.getcwd()
    sim_basename = os.path.basename(sim_loc)
    outf = open('%s/resubmit.sbatch'%sim_loc, 'w')
    outf.write('#!/bin/bash\n')
    outf.write('#SBATCH --job-name=%s\n'%sim_basename)
    outf.write('#SBATCH --nodes=%i\n'%nodes)
    outf.write('#SBATCH --ntasks-per-node=%i\n'%(tasks_per_node))
    outf.write('#SBATCH --time=%s\n\n'%wall_time)

    outf.write('module purge\n')
    outf.write('module load slurm gcc lib/hdf5 openmpi\n\n')

    outf.write('cd %s\n'%sim_loc)

    outf.write('c=0\nwhile true\ndo\n')
    outf.write('    if [ ! -e estd.out.$c ]; then\n        break\n    fi\n')
    outf.write('    c=$[$c+1]\ndone\nmv estd.out estd.out.$c\n\n')

    procs = nodes * tasks_per_node
    outf.write('mpirun -n %i ~/enzo-dev-cr/src/enzo/enzo.exe -d -r %s/%s >&estd.out\n'\
               %(procs, output, output))
    outf.close()


def find_last_output(sim_loc):
    DD_list = sorted(glob.glob('%s/DD*'%sim_loc))
    # sometimes the last output is corrupt
    return os.path.basename(DD_list[-2])


def resubmit_unfinished_sims(nodes = 1, tasks_per_node = 32, wall_time = '4:00:00'):
    unfinished_sim_list = np.loadtxt('unfinished_sims.dat', dtype = 'str')
    workdir = unfinished_sim_list[0]
    unfinished_sim_list = unfinished_sim_list[1:]

    for unfinished_sim in unfinished_sim_list:
        sim_loc = '%s/%s'%(workdir, unfinished_sim)
        last_output = find_last_output(sim_loc)
        generate_restart_sbatch_file(sim_loc, last_output, nodes = nodes, \
                                     tasks_per_node = tasks_per_node, wall_time = wall_time)

        os.system('sbatch %s/resubmit.sbatch'%sim_loc)
resubmit_unfinished_sims()

    