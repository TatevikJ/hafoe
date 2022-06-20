#!/bin/bash
#SBATCH --mem 20gb
#SBATCH --output=slurm-%j.log
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=tatevik.jalatyan@abi.am


./hafoe.sh \
    -parentlib input_files/AAV_all16_new.fasta \
    -chimericlib input_files/Chimeric_lib_simulated.csv \
    -enrichedlib1 input_files/Enriched_lib_simulated.fastq.gz  \
    -o hafoe_out \
    -cdhitest /storage2/proj/kodikaz/softwares/cdhit/cd-hit-est \
    -cdhitest2d /storage2/proj/kodikaz/softwares/cdhit/cd-hit-est-2d \
    --explore \
    --identify \
    --overlap \
    -readlength 100 \
    -stepsize 15




