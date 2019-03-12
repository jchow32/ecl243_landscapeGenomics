#!/bin/bash
#SBATCH -p med
#SBATCH -t 7-10:00:00
#SBATCH --mem=10000
#SBATCH -n 1

/home/thiagoms/WIFL_rad/flmm/LFMM_CL_v1.5/bin/LFMM -x cut_transposed_VCF_WIFL_input -v LFMM_environ_input -K 4
