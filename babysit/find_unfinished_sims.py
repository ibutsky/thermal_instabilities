import os
import glob

def find_last_output(sim_loc):
    DD_list = sorted(glob.glob('%s/DD*'%sim_loc))
    # sometimes the last output is corrupt
    return os.path.basename(DD_list[-1])

workdir = '../../simulations/production'

sim_loc_list = glob.glob('%s/*'%workdir)
unfinished_sim_list = []

for sim_loc in sim_loc_list:
    if os.path.isdir('%s/DD0000'%sim_loc) and not os.path.isdir('%s/DD0100'%sim_loc):
        unfinished_sim_list.append(os.path.basename(sim_loc))
                                   

outf = open('unfinished_sims.dat', 'w')
outf.write('%s\n'%workdir)
for unfinished_sim in unfinished_sim_list:
    last_output = find_last_output('%s/%s'%(workdir, unfinished_sim))
    outf.write('%s, %s\n'%(last_output, unfinished_sim))

outf.close()
