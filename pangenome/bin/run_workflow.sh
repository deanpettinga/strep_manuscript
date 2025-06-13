#!/bin/bash

#SBATCH --time=24:00:00           # walltime limit (HH:MM:SS)
#SBATCH --nodes=1                    # number of nodes
#SBATCH --ntasks-per-node=1          # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=8G
#SBATCH --partition=short
#SBATCH --qos=pgec
#SBATCH --job-name=streps

# load snakemake
module load miniconda
source ~/.bashrc
source activate snakemake

# save DAG job file with time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")

# make logs dir if it does not exist already. Without this, logs/ is automatically generate only after the first run of the pipeline
logs_dir="logs/runs"
[[ -d $logs_dir ]] || mkdir -p $logs_dir

# log the jobs to be run
snakemake -s Snakefile --use-conda -n > logs/runs/workflow_${TIME}.txt
# log a png of the DAG
snakemake -s Snakefile --use-conda --dag | dot -Tpng > logs/runs/workflow_${TIME}.png

snakemake \
-s Snakefile \
--rerun-incomplete \
--cores 1 \
--use-conda \
--use-envmodules \
--jobs 48 \
--executor cluster-generic \
--cluster-generic-submit-cmd "sbatch \
    --time={resources.time} \
    -N {resources.nodes} \
    -n {resources.threads} \
    --mem={resources.mem_gb}G \
    --partition={resources.partition} \
    --qos={resources.qos} \
    --job-name={resources.name} \
    --error=logs/slurm/{resources.name}_%j.log \
    --out=logs/slurm/{resources.name}_%j.log"
