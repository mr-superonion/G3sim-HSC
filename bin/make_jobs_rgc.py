"""
This script generates a load of PBS scripts that can be submitted on the Coma cluster at CMU, for
the purpose of mass simulation production.  It also writes a bash script that can be run to submit
all the PBS scripts.
"""
import galsim
import os

# Settings go up here.
njobs = 10 # Submit this many PBS scripts.
output_dir = '/gpfs02/work/xiangchong.li/s18SimHSC/s16Repeat/output'
script_dir = 'pbs_scripts'
script_pref = 'job_'
bash_file = 'submit_all.x'

# Make sure we have the right version of GalSim installed.  If not, throw an exception.
# Note Add 1.5.1
if galsim.__version__ not in ['1.4.3']:
    raise RuntimeError("Wrong version of GalSim installed: %s"%(galsim.__version__))

# Check for the files needed to analyze the sims:
# - the COSMOS 25.2 sample
# - the variance file
# - the PSF file

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
        f.write("#PBS -l nodes=1:ppn=20\n")
        f.write("#PBS -q small\n")
        f.write("#PBS -j oe\n")
        f.write("#PBS -N "+jobname+"\n")
        f.write("#PBS -l walltime=100:00:00\n")
        f.write("\n")
        f.write("source ~/.setupHSC6.sh\n")
        f.write("cd %s\n"%os.getcwd())
        f.write("galsim -v 1 -n %d -j %d rgc_hsc.yaml output.noclobber=True output.dir=%s\n"%(njobs,jobnum+1,output_dir))
    f.close()

# Finally, write a bash script that just submits all the PBS files.
f_bash.close()

