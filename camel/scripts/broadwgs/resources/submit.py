"""
Script to launch a snakemake run on slurm.
Run as: snakemake -s /scratch/thdelcourt/testing_slurm/camel_3.0_TD/camel/scripts/broadwgs/snakefile/main.smk --configfile /scratch/thdelcourt/testing_slurm/working/config.yml --cores 1 --cluster "python3 /scratch/thdelcourt/testing_slurm/submit.py {dependencies}" --is --notemp --jobs 1000 --local-cores 1

 - Adds the dependencies for immediate submission of jobs
 - Wraps around slurm's sbatch command
 - Adds a delay before starting the job to avoid a bug in snakemake due to immediate submission (cf github; ongoing)
 - Assigns the correct memory, cpu and temp disk usage if provided in the resources dict.
"""

#!/usr/bin/env python3
import os
import sys

from snakemake.utils import read_job_properties

# last command-line argument is the job script
jobscript = sys.argv[-1]

# all other command-line arguments are the dependencies
dependencies = list(sys.argv[1:-1])

# parse the job script for the job properties that are encoded by snakemake within
job_properties = read_job_properties(jobscript)
print(job_properties)

# collect all command-line options in an array
cmdline = ["sbatch"]

# set all the slurm submit options as before
print(job_properties["cluster"])
# slurm_args = " -p {partition} -N {nodes} -n {ntasks} -c {ncpus} -t {time} -J {job-name} -o {output} -e {error} ".format(**job_properties["cluster"])

# setting memory. setting to 0 grants all memory on the machine, which is the exact opposite.
if 'mem_mb' in job_properties["resources"]:
    mem_mb = job_properties["resources"]["mem_mb"]
else:
    mem_mb = "1"

threads = job_properties["threads"]

# setting tmpdsk. Empty if not set in resources.
grestmpdsk  = ""
if 'tmpdsk' in job_properties["resources"]:
    tmpdsk = job_properties["resources"]["tmpdsk"]
    grestmpdsk = " --gres=tmpdsk:{tmpdsk}"

# adding 15s to begin time to ensure jobs don't get jumbled up before running.
#time_to_begin = dt.datetime.fromtimestamp(time.time()+20).strftime("%H:%M:%S")
#slurm_args = f" --mem {mem_mb}M -c {threads} -b {time_to_begin}"
#slurm_args = f" --mem {mem_mb}M -c {threads} --begin=now+15 --gres=tmpdsk:100000 "
slurm_args = f" --mem {mem_mb}M -c {threads} --begin=now+15 {grestmpdsk} "
#slurm_args = f" --mem {mem_mb}M -c {threads} --gres=tmpdsk:600000 "

#'threads': 2, 'resources': {'mem_mb': 50000},
cmdline.append(slurm_args)

if dependencies:
    cmdline.append("--dependency")
    # only keep numbers in dependencies list
    # with open("/scratch/slurm/dependencies.txt", "a") as out_str:
    #     out_str.write(str(dependencies)+"\n")
    dependencies = [ x for x in dependencies if x.isdigit() ]
    cmdline.append("afterok:" + ",".join(dependencies))

cmdline.append(jobscript)
# with open("/scratch/slurm/dependencies.txt", "a") as out_str:
#         out_str.write(str(cmdline)+"\n")

os.system(" ".join(cmdline))

