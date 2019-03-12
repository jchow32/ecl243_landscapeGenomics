#### WIFL ECL 243 ####

#### Downloading datasets and scripts ####

# Download scripts and subset of genomic vulnerability data to /home/thiagoms/'genomic vulnerability subset and scripts', genome assembly

scp -P 2022 -r genomic\ vulnerability\ subset\ and\ scripts/ thiagoms@farm.cse.ucdavis.edu:/home/thiagoms

# Download willow flycatcher SRA's start from SRR7064646 to SRR7064864

# import to /home/thiagoms/WIFL_rad 219 sequences of flycatcher
for (( i = 646; i <= 864; i++ ))
  do
  wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR706/SRR7064$i/SRR7064$i.sra
done

# SRA to fastq

sbatch sra_to_fastq.sh

# Stacks used to demultiplex and discard low quality with process_radtags ; remove duplicate read pairs with clone_filter: *stacks_demulti.sh*

sbatch stacks_demulti.sh

# Get the genome assembly, see https://www.ncbi.nlm.nih.gov/assembly/GCF_003031625.1/

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/031/625/GCF_003031625.1_ASM303162v1//GCF_003031625.1_ASM303162v1_genomic.fna.gz

#### SNP association with environment ####

# Build index, maps reads to genome assembly using bowtie2, convert SAM to BAM and sort, use Haplotype Caller from GATK to get SNPs 
# Note that the uncommented lines were not run for the sake of time. VCF file of 105,000 SNPs was provided

sbatch alignment.sh

# Filtering occurs in alignment.sh: remove genotype quality < 30, depth < 8, minor allele frequency < 0.01, indels and non-biallelic SNPs via vcftools

######################

# After receiving the VCF.

# Get the packages necessary to run gradientForest (extendedForest)

wget http://download.r-forge.r-project.org/src/contrib/extendedForest_1.6.1.tar.gz

tar -xzvf extendedForest_1.6.1.tar.gz

# Have to grep some R packages dependencies
# Uncomment them inside the tutorial_GI.R if needed

# tutorial_GI.R: The resulting .pdf contains some of the figures from the paper... but was made from only using 10,000 SNPs rather than all 37855 variants

# So we need to format /home/thiagoms/genomic\ vulnerability\ subset\ and\ scripts/wiflforest.allfreq.sample.csv as 

# meaning of fields in VCF per individual
# GT:AD:DP:GQ:PL
# Genotype
# Allelic depths for the ref and alt alleles in the order listed
# Approximate read depth (reads with MQ=255 or with bad mates are filtered)
# Genotype Quality
# Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification

# Get the IDs of individuals and associate them with their respective populations (1 population per line)
for i in {1..24}; do awk -v var="$i" -F$'\t' '{if ($11 == var) {print $1}}' WIFL_175inds_meta.txt | tr '\n' ',' | sed  's/.$/\n/'; done > cats

# append vcf-subset command
while IFS='' read -r line || [[ -n "$line" ]]; do echo "vcf-subset -c $line WIFL_final.vcf > subset_"; done < cats > vcf_subset.sh

# attach SBATCH headers 
cat <(head -n5 sra_to_fastq.sh) <(paste vcf_subset.sh <(for i in {1..24}; do echo $i; done) | sed 's/\t//g') > vcf_subset_head.sh
sbatch vcf_subset_head.sh

# get vcf from individuals within each population separately
for i in {1..24}; do vcftools --vcf "/home/thiagoms/WIFL_rad/subset_${i}" --freq --out "subset_${i}"; done

# get allele frequencies for a population (1 population per file) as a row, for all 105,000 variants (in columns)
for i in {1..24}; do cut -f5 "subset_${i}.frq" | cut -d: -f2 | tr '\n' ',' | sed 's/FREQ},//g' | sed '$ s/.$/\n/' > "subset_${i}_row"; done

# cat together all the populations (variants in columns, populations in rows)
for i in {1..24}; do cat "subset_${i}_row" >> subset_all_rows; done

# run gradient Forest for 24 populations using 20 ClimWorld variables

###### !!!!INSERT HOW DOWNLOADED CLIMWORLD VARIABLES!!!! ######

sbatch gradientForestTutorial.sh

######################

# Format files as input for LFMM

# Get genotype matrix for individuals, where in 0/0 in VCF means homozygous for major allele, 0/1 heterozygous, 1/1 homozygous alternate. Individuals (columns), variants (rows)
cut -d$'\t' -f10- WIFL_final.vcf | grep -v "##" | sed 's/:[0-9]*,[0-9]*:[0-9]*:[0-9]*:[0-9]*,[0-9]*,[0-9]*//g' | sed 's/0\/0/0/g' | sed 's/1\/0/1/g' | sed 's/1\/1/2/g' | sed 's/0\/1/1/g' | sed 's/\.\/\.:\.:\.:\.:\./9/g' | sed 's/\.\/\./9/g' > genotypes_individuals_vcf

###### !!!!INSERT HOW TRANSPOSED genotypes_individuals_vcf!!!! ######

# remove individual who do not belong to any population, take only the fields containing genotype information
grep -v 1590-97279 transposed_VCF_WIFL | grep -v 1590-97493 | grep -v 1710-20525 | grep -v 1710-20526 | cut -f2- | tail -n+2 > cut_transposed_VCF_WIFL_input


# To get environmental matrix for individuals, copy same environmental data

# make WIFL_BIO environmental variables easier to read
sed 's/"//g' WIFL_Bio.csv | sed 's/,/\t/g' > cleaned_WIFL_Bio.tsv

# cut -f11 WIFL_175inds_meta.txt | tail -n+2 | tr '\n' ',' | sed '$ s/.$/\n/'
# above gives the list, which is pasted below
for i in {3,3,3,2,2,2,1,3,3,11,11,NA,13,13,24,3,3,3,3,15,15,15,15,NA,21,21,21,12,4,4,4,4,2,2,2,2,2,8,8,11,11,11,11,11,11,11,11,NA,NA,13,13,13,13,24,24,24,1,1,3,3,3,1,17,17,17,17,17,17,17,12,4,4,1,15,22,22,22,22,23,23,23,23,15,15,15,1,1,1,16,16,16,16,16,16,1,1,2,2,2,2,2,3,2,10,10,10,10,18,18,18,18,12,12,12,12,15,15,15,9,9,9,9,9,9,8,8,8,8,8,7,7,7,7,7,7,7,7,7,7,7,7,7,5,5,19,5,5,6,6,19,19,19,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,20,20}; do awk -v var="$i" '{if (var == $1) {print $0}}' cleaned_WIFL_Bio.tsv; done > LFMM_environ_input

# Download and install LFMM command line

cd /home/thiagoms/WIFL_rad/flmm
wget http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/LFMM_CL_v1.5.zip
unzip LFMM_CL_v1.5.zip
./install.command

# Run LFMM using genotypic and environmental matrices
sbatch lfmm.sh

# Create a tab-separated file which contains the chromosome, position, and adjusted p-value for each variant for all 20 environmental variables (So there is one variant per row, and columns are the 20 environmental variables)
# Get the position of each variant including its scaffold
cut -d$'\t' -f1,2 /home/thiagoms/WIFL_rad/WIFL_final.vcf | grep -v "##" | tail -n+2  > /home/thiagoms/WIFL_rad/flmm/LFMM_adj_pvalues_position 

# Get the uncorrected p-values from LFMM
for i in {1..20}; do cut -d' ' -f2 "cut_transposed_VCF_WIFL_input_s${i}.4.zscore" > "cut_transposed_VCF_WIFL_input_s${i}.4.zscore_pvalues"; done

# Some scaffolds don't have LFMM results for certain environmental variables, so have to read in each vector of p-values individually (column) for FDR correction
sbatch fdr_correction.sh

# for i in {1..20}; do echo "cut_transposed_VCF_WIFL_input_s${i}.4.zscore_pvalues_adj"; done | tr '\n' ' '
# above is pasted to create command below

paste /home/thiagoms/WIFL_rad/flmm/LFMM_adj_pvalues_position cut_transposed_VCF_WIFL_input_s1.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s2.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s3.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s4.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s5.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s6.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s7.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s8.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s9.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s10.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s11.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s12.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s13.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s14.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s15.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s16.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s17.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s18.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s19.4.zscore_pvalues_adj cut_transposed_VCF_WIFL_input_s20.4.zscore_pvalues_adj > results_adj_nlog10_allVar

# Correction for 2.2 million comparisons

sbatch fdr_correction_22million.sh

# for i in {1..20}; do echo "cut_transposed_VCF_WIFL_input_s${i}.4.zscore_pvalues_adj_22"; done | tr '\n' ' '
# above is pasted to create command below
paste /home/thiagoms/WIFL_rad/flmm/LFMM_adj_pvalues_position cut_transposed_VCF_WIFL_input_s1.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s2.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s3.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s4.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s5.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s6.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s7.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s8.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s9.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s10.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s11.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s12.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s13.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s14.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s15.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s16.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s17.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s18.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s19.4.zscore_pvalues_adj_22 cut_transposed_VCF_WIFL_input_s20.4.zscore_pvalues_adj_22 > results_adj_nlog10_allVar_22
