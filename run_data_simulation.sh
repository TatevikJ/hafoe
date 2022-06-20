#!/bin/bash
#SBATCH --mem 40gb
#SBATCH --output=slurm-%j.log
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=tatevik.jalatyan@abi.am

Rscript data_simulation/data_simulation.R
