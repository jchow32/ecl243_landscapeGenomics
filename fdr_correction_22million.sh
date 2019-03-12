#!/bin/bash
#SBATCH -p med
#SBATCH -t 7-10:00:00
#SBATCH --mem=10000
#SBATCH -n 1

for i in {1..20}; do Rscript fdr_correction.R "cut_transposed_VCF_WIFL_input_s${i}.4.zscore_pvalues" "cut_transposed_VCF_WIFL_input_s${i}.4.zscore_pvalues_adj_22"; done
