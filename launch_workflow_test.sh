#!/bin/bash
#SBATCH --account awesolo2
#SBATCH --partition shared
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -o ./slurm/slurm-%j.out

module load java

nextflow run main.nf \
    --mode test \
    -profile singularity,rockfish \
    --full_sim_dir /data/awesolo2/projects/GenE8_rel_sims/simulations \
    --resdir /data/awesolo2/projects/GenE8_rel_sims/results

module unload java
