#!/usr/bin/bash

#SBATCH --time 24:00:00
#SBATCH --partition batch
#SBATCH --job-name workflow_submission
#SBATCH --output=logs/workflow_submission_%j.log

# Make sure to do the following before running this script:
# - conda activate snakemake environment
# - check that slurm profile for snakemake is set up

# set the number of jobs to run at a time (no spaces)
num_jobs=50

# if using cluster.json, run this command:
# snakemake -j $num_jobs --use-conda --rerun-incomplete --latency-wait 20 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -N {cluster.N}  -t {cluster.t} -o {cluster.o} -e {cluster.e} -J {cluster.J} -c {cluster.c} --mem {cluster.mem}" -s Snakefile 

# if using cluster.yaml (highly recommended), run this command:
snakemake -j $num_jobs \
--verbose \
--use-conda \
--conda-prefix $CONDA_PREFIX_1/envs \
--profile slurm \
--cluster-config cluster.yaml

exit
