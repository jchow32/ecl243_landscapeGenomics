#!/bin/bash
#SBATCH -p med
#SBATCH -t 7-10:00:00
#SBATCH --mem=4056
#SBATCH -n 1

/share/apps/stacks-2.2/bin/process_radtags -p /home/thiagoms/WIFL_rad/WIFL_fastq -o /home/thiagoms/WIFL_rad/WIFL_processed -e sbfI -i fastq -r -c -q

/share/apps/stacks-2.2/bin/clone_filter -p /home/thiagoms/WIFL_rad/WIFL_processed -i fastq -o /home/thiagoms/WIFL_rad/WIFL_processed_clone
