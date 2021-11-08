"""
Script to launch a snakemake run on slurm.
Run as: snakemake -s /scratch/thdelcourt/testing_slurm/camel_3.0_TD/camel/scripts/broadwgs/snakefile/main.smk --configfile \
/scratch/thdelcourt/testing_slurm/working/config.yml --cores 1 --cluster "python3 /scratch/thdelcourt/testing_slurm/submit.py \
{dependencies}" --is --notemp --jobs 1000 --local-cores 1

 - Adds the dependencies for immediate submission of jobs
 - Wraps around slurm's sbatch command
 - Adds a delay before starting the job to avoid a bug in snakemake due to immediate submission (cf github; ongoing)
 - Assigns the correct memory, cpu and temp disk usage if provided in the resources dict.
"""

import os
import sys

from snakemake.utils import read_job_properties

# last command-line argument is the job script
jobscript = sys.argv[-1]

# all other command-line arguments are the dependencies
dependencies = list(sys.argv[1:-1])

# parse the job script for the job properties that are encoded by snakemake within
job_properties = read_job_properties(jobscript)

# collect all command-line options in an array
cmdline = ["sbatch"]

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
slurm_args = f" --mem {mem_mb}M -c {threads} --begin=now+60 {grestmpdsk} "
cmdline.append(slurm_args)

if dependencies:
    cmdline.append("--dependency")
    dependencies = [ x for x in dependencies if x.isdigit() ]
    cmdline.append("afterok:" + ",".join(dependencies))

cmdline.append(jobscript)
os.system(" ".join(cmdline))
