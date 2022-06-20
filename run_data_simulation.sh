#!/bin/bash
#SBATCH --mem 40gb
#SBATCH --output=slurm-%j.log

Rscript data_simulation/data_simulation.R
