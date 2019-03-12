#!/bin/bash
#SBATCH -p med
#SBATCH -t 7-10:00:00
#SBATCH --mem=25000
#SBATCH -n 1

#Rscript /home/thiagoms/WIFL_rad/tutorial_GI.R

#mv Rplots.pdf 21popUsed.pdf

Rscript /home/thiagoms/WIFL_rad/tutorial_GI_24pop.R

mv Rplots.pdf 24popUsed.pdf
