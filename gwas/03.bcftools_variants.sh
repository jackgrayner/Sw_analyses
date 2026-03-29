#!/bin/bash
#SBATCH --job-name=freebayes_auto  
#SBATCH --cpus-per-task=16   
#SBATCH --mem=16gb                
#SBATCH --partition=long
#SBATCH --array=1-14
#SBATCH --output=bcftools_array_%A_%a.log

chr=${SLURM_ARRAY_TASK_ID}
source activate bcftools
bcftools mpileup -Ov -f TOC.asm.scaffold.fasta -d 100000 --annotate FORMAT/AD,FORMAT/DP -b bamlist.txt -r scaffold_${chr}
