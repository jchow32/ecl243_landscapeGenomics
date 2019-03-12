#!/bin/bash
#SBATCH -p med
#SBATCH -t 24:00:00
#SBATCH --mem=4056
#SBATCH -n 1

for i in /home/thiagoms/WIFL_rad/*.sra
do
	/share/apps/sra-toolkit-2.8.2/bin/fastq-dump.2.8.2 --split-3 -O /home/thiagoms/WIFL_rad $i
done
