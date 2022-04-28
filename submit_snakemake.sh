#!/usr/bin/bash
#SBATCH --time 35:00:00
#SBATCH --partition exacloud
#SBATCH --job-name workflow_submission
#SBATCH --output=logs/workflow_submission_%j.log

# source activate omic-qc-wf
# if using cluster.json, run this command:
# snakemake -j 100 --use-conda --rerun-incomplete --latency-wait 20 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -N {cluster.N}  -t {cluster.t} -o {cluster.o} -e {cluster.e} -J {cluster.J} -c {cluster.c} --mem {cluster.mem}" -s Snakefile 

# if using cluster.yaml (highly recommended), run this command instead:
snakemake -j 100 --use-conda --profile slurm --cluster-config cluster.yaml

exit