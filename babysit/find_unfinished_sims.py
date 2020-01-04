import os
import glob

workdir = '../../simulations/production'

sim_loc_list = glob.glob('%s/*'%workdir)
unfinished_sim_list = []

for sim_loc in sim_loc_list:
    if not os.path.isdir('%s/DD0100'%sim_loc):
        unfinished_sim_list.append(os.path.basename(sim_loc))
                                   

outf = open('unfinished_sims.dat', 'w')
outf.write('%s\n'%workdir)
for unfinished_sim in unfinished_sim_list:
    outf.write('%s\n'%unfinished_sim)

outf.close()
