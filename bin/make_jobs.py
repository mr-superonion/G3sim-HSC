#!/usr/bin/env python
import galsim
import os

# Settings go up here.
njobs = 1 # Submit this many PBS scripts.
script_dir = 'pbs_scripts'
script_pref = 'job_'
bash_file = 'submit_all.x'


# Generate the PBS files in a subdirectory, for both Wide and Deep.
if not os.path.exists(script_dir):
    os.mkdir(script_dir)

f_bash = open(bash_file, "w")

for jobnum in range(njobs):
    jobname = '%s%03d'%(script_pref,jobnum)
    script_file = os.path.join(script_dir,jobname+'.sh')
    f_bash.write("qsub %s\n"%script_file)
    with open(script_file, "w") as f:
        f.write("#!/bin/sh\n")
        f.write("#PBS -l nodes=2:ppn=20\n")
        f.write("#PBS -q small\n")
        f.write("#PBS -j oe\n")
        f.write("#PBS -N "+jobname+"\n")
        f.write("#PBS -l walltime=100:00:00\n")
        f.write("\n")
        f.write("source ~/.setupHSC6.sh\n")
        f.write("cd %s\n"%os.getcwd())
    f.close()

# Finally, write a bash script that just submits all the PBS files.
f_bash.close()

