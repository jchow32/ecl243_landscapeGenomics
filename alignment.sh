#!/bin/bash
#SBATCH -p ecl243
#SBATCH -t 7-10:00:00
#SBATCH --mem=10000
#SBATCH -n 8

bowtie2-build --threads 8 /home/thiagoms/WIFL_rad/GCF_003031625.1_ASM303162v1_genomic.fna.gz flycatcher

#bowtie2 -t -x flycatcher -U [list of files with reads] -S aligned_flycatcher.sam

#samtools view -bS aligned_flycatcher.sam > aligned_flycatcher.bam

#samtools sort aligned_flycatcher.bam -o aligned_flycatcher.sorted.bam

#/share/apps/gatk-4.0.11.0/gatk HaplotypeCaller -R /home/thiagoms/WIFL_rad/GCF_003031625.1_ASM303162v1_genomic.fna.gz -I aligned_flycatcher.sorted.bam -O flycatcher.g.vcf.gz -ERC GVCF

# vcftools --gzvcf flycatcher.g.vcf.gz --minGQ 30 --minDP 8  --maf 0.01 --remove-indels --min-alleles 2 --max-alleles 2
