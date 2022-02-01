#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH -c 10
#SBATCH --job-name="test"
#SBATCH --output=test.txt
#SBATCH --mail-user=amina.lemsara@uni-heidelberg.de
#SBATCH --partition=general

source activate snakemake
module load java
module load samtools
srun snakemake --cores --unlock

#For 3way analysis
srun snakemake --cores all   analysis_aggregate

#For downsampling analysis
# srun snakemake --cores all -s /home/alemsara/DirectRNA/rRNA_by_ONT/Computational\ Analysis/snakefile  downsampling_aggregate

#For WTvsIVT analysis
# srun snakemake --cores all -s /home/alemsara/DirectRNA/rRNA_by_ONT/Computational\ Analysis/snakefile IVTvsWT_analysis

#For mixing analysis
# srun snakemake --cores all agregate_mixing --rerun-incomplete

#Heart Biopsies
# srun snakemake --cores all bivariate_analysis 
