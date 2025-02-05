#!/bin/bash

# author: Susanne Voigt

wdir=/mnt/c/Users/ennas/Desktop/mito/mito_seq		#set working directory  
ref_folder=$wdir/mito_mt_seq_27/mito_mt_seq_ref	#directory to reference
ref=$ref_folder/Dmel_6.13.fasta						#reference
bam_folder=$wdir/mito_mt_seq_27/mito_mt_seq_bam	#directory to .bam data

#----------------
#- SNP calling -
#---------------

# make directory for output
snp_folder=$wdir/mito_mt_seq_27/mito_mt_seq_snps
mkdir $snp_folder

# get sample names
cd $bam_folder
sample_name=($(ls -1 | grep bam$ | cut -d'_' -f 1)) # sample names
cd $wdir

# .bam files with only mt genome
for i in "${!sample_name[@]}"; do 
samtools view $bam_folder/${sample_name[$i]}_V1_dup.bam mitochondrion_genome -b > $bam_folder/${sample_name[$i]}_mt.bam
done

# list all bam files 
bam_mt=($(ls -d -1 $bam_folder/* | grep mt.bam$))

# merge all .bam files in single mpileup file (only retaining nucleotides with BQ >20 and reads with MQ >20)
samtools mpileup -B -q 20 -Q 20 -f $ref ${bam_mt[@]} | awk '/mitochondrion_genome/' > $snp_folder/mito_mt_seq_27.mpileup

# remove mt .bam files
rm $bam_folder/*_mt.bam

# mt only ref
samtools faidx $ref mitochondrion_genome > $ref_folder/dmel_ref_mt.fasta

# call SNPs with PoolSNP (https://github.com/capoony/PoolSNP - cite doi:10.1093/molbev/msaa120)
./scripts/PoolSNP-master/PoolSNP.sh mpileup=$snp_folder/mito_mt_seq_27.mpileup output=$snp_folder/mito_mt_seq_27 reference=$ref_folder/dmel_ref_mt.fasta \
names=${!sample_name[@]} \
max-cov=0.9999 min-cov=10 min-count=10 min-freq=0.001 miss-frac=0.2 jobs=6 BS=1 
	#max-cov: maximum coverage threshold based on a percentile cutoff from the coverage distribution which is calculated for each chromosomal arm/scaffold and all sample separately
	#min-cov: minimum coverage threshold which is tested for each sample separately
	#min-count: minimum allele count across all libraries pooled
	#min-freq: minimum allele frequency across all libraries pooled
	#miss-frac: a consistency parameter that defines how many libraries need to fulfill all of the above-mentioned threshold parameters so that a site is considered
	#jobs: number of parallel jobs/cores used for the SNP calling
	#
	#PoolSNP creates multiple output files:
	# - A gzipped VCF file (v.4.2) containing allele counts and frequencies for every position and library
	# - A Max-coverage file containing the maximum coverage thresholds for all chromosomal arms and libraries in the mpileup file (separated by a column)
	# - Optionally a Bad-sites file (by setting the parameter BS=1), which contains a list of (variable and invariable) sites that did not pass the SNP calling criteria. This file can be used to weight windows for the calculation of population genetic estimators with PoolGEN



#-------------------------------------------------------------------------------------
#- genetic differentiation: allele frequency differences (Fisher's exact test) & FST -
#-------------------------------------------------------------------------------------

# convert .mpilup to .sync
java -jar scripts/popoolation2_master/mpileup2sync.jar \
--min-qual 20 --threads 6 --input $snp_folder/mito_mt_seq_27.mpileup --output $snp_folder/mito_mt_seq_27_tmp.sync 
	#.sync format:
	#-------------
	#col1: reference contig (chromosome)
	#col2: position in the reference contig
	#col3: reference character
	#col4: population 1
	#col5: population 2
	#coln: population n
	#
	#population data are in the form
	#A:T:C:G:N:*
	#A: count of character A
	#T: count of character T
	#C: count of character C
	#G: count of character G
	#N: count of character N
	#*: deletion, count of deletion

#remove del positions in .sync file, ie. set them  to 0 so that sites are not removed from popoolation analysis
awk 'BEGIN {OFS="\t"}{for(i=4;i<=NF;i++) sub(/[0-9]+$/,0,$i)}1' $snp_folder/mito_mt_seq_27_tmp.sync \
> $snp_folder/mito_mt_seq_27.sync
rm $snp_folder/mito_mt_seq_27_tmp.sync 


# calculate allele frequency differences/Fisher's exact tests
perl scripts/popoolation2_master/fisher-test.pl \
-window-size 1 --step-size 1  --min-covered-fraction 1.0 --min-coverage 10 \
--max-coverage 0.0001% --min-count 1%  --suppress-noninformative \
--input $snp_folder/mito_mt_seq_27.sync \
--output $snp_folder/mito_mt_seq_27.fet
	#--minimum coverage: the coverage in each of the populations has to be higher or equal to this threshold
	#--min-count: % minimum count of the minor allele across all populations
	#--max-coverage: % highest coverages will be ignored, this value is independently estimated for every population
	#--window-size 1 --step-size 1  --min-covered-fraction 1.0 => each site is considered
	#for further details see perl script
	#
	#output format:
	#col1: reference contig (chromosome)
	#col2: mean position of the sliding window
	#col3: number of SNPs found in the window (not considering sites with a deletion) 
	#col4: fraction of the window which has a sufficient coverage (min. coverage <= cov <= max. coverage) in every population;
	#col5: average minimum coverage in all populations
	# ....
	#1:2 -log10(product of p-value) for comparision of population 1 and 2
	#1:3 -log10(product of p-value) for comparision of population 1 and 3
	#1:4 -log10(product of p-value) for comparision of population 1 and 4
	#
	#for further details see perl script

# convert output into .csv - all SNPs 
sed 's/[0-9]\{1,\}:[0-9]\{1,\}=//g' $snp_folder/mito_mt_seq_27.fet | 
sed 's/\t/,/g' > $snp_folder/mito_mt_seq_27_fet.csv
# convert output into .csv - only sig. SNPs (ie. p-value < 1e-12 or -log10(p-value)>12)
awk '/\=(1[2-9]|[2-9][0-9]|[0-9]{3,})\./' $snp_folder/mito_mt_seq_27.fet | sed 's/[0-9]\{1,\}:[0-9]\{1,\}=//g' |
sed 's/\t/,/g' > $snp_folder/mito_mt_seq_27_fet_sig.csv

# calculate FST (according to Hudson et al 1992 https://doi.org/10.1093/genetics/132.2.583)
perl scripts/popoolation2_master/fst-sliding.pl \
--window-size 1 --step-size 1  --min-covered-fraction 1.0 --min-coverage 10 \
--max-coverage 0.0001% --min-count 1% --pool-size 150 --suppress-noninformative \
--input $snp_folder/mito_mt_seq_27.sync \
--output $snp_folder/mito_mt_seq_27.fst
	#--minimum coverage: the coverage in each of the populations has to be higher or equal to this threshold
	#--min-count: % minimum count of the minor allele across all populations
	#--max-coverage: % highest coverages will be ignored, this value is independently estimated for every population
	#--window-size 1 --step-size 1  --min-covered-fraction 1.0 => each site is considered
	#--pool-size: size of the population pools; May be provided for each population individually (number of sampled chroms per pop)
	#
	#output format:
	#col1: reference contig (chromosome)
	#col2: mean position of the sliding window
	#col3: number of SNPs found in the window (not considering sites with a deletion) 
	#col4: fraction of the window which has a sufficient coverage (min. coverage <= cov <= max. coverage) in every population;
	#col5: average minimum coverage in all populations
	# ....
	#1:2 the pairwise Fst for population 1 and 2
	#1:3 the pairwise Fst for population 1 and 3
	#
	#for further details see perl script

# convert FST output to .csv
sed 's/[0-9]\{1,\}:[0-9]\{1,\}=//g' $snp_folder/mito_mt_seq_27.fst | 
sed 's/\t/,/g' > $snp_folder/mito_mt_seq_27_fst.csv
# extracts sig diff SNPs (as given by Fisher's exact tests; p-value < 1e-12) from FST output 
awk -F, '{print $2}' $snp_folder/mito_mt_seq_27_fet_sig.csv > $snp_folder/sig_pos_tmp
mapfile -t a < $snp_folder/sig_pos_tmp
for i in ${a[@]}; do awk -F, -v i="$i" '($2==i)' $snp_folder/mito_mt_seq_27_fst.csv; done \
> $snp_folder/mito_mt_seq_27_sig_fst.csv


#---------------------------------------
#- allele frequencies & major alleles -
#---------------------------------------

# calculate allele frequencies
perl scripts/popoolation2_master/snp-frequency-diff.pl \
--min-coverage 10 --max-coverage 0.0001%  --min-count 1% \
--input $snp_folder/mito_mt_seq_27.sync \
--output $snp_folder/mito_mt_seq_27_freq
	##2 output files:
	#"output_prefix"_pwc: contains the differences in allele frequencies for all pairwise comparisions
	#"output_prefix"_rc: shows the major and minor allele in a succinct format for all populations

python2 scripts/DrosEU_pipeline-master/scripts/sync2AF.py \
--inp $snp_folder/mito_mt_seq_27.sync --out $snp_folder/mito_mt_seq_27
#using python script from DrosEU pipeline - https://github.com/capoony/DrosEU_pipeline - cite doi: 10.1093/molbev/msaa120

# extract major alleles - sig. SNPs
awk 'NR>1' $snp_folder/mito_mt_seq_27_freq_rc |awk 'BEGIN {OFS=","}{print $2,$8}' | 
sed 's/A/A,/g' | sed 's/T/T,/g' | sed 's/C/C,/g' | sed 's/G/G,/g'| sed 's/,$//g' \
> $snp_folder/mito_mt_seq_27_major_alleles_tmp
#
header=$( IFS=$','; echo "${sample_name[*]}" )
mapfile -t a < $snp_folder/sig_pos_tmp
for i in ${a[@]}; do awk -F, -v i="$i" '($1==i)' $snp_folder/mito_mt_seq_27_major_alleles_tmp; done |
awk -v var="$header" 'BEGIN {print ","var}{print}'> $snp_folder/mito_mt_seq_27_major_alleles_sig.csv
#


# create snp table - .csv major alleles (consensus)
awk -F, '
{
    for (i = 1; i <= NF; i++) {
        if(NR == 1) {
            s[i] = $i;
        } else {
            s[i] = s[i] " " $i;
        }
    }
}
END {
    for (i = 1; s[i] != ""; i++) {
        print s[i];
    }
}' $snp_folder/mito_mt_seq_27_major_alleles_sig.csv |
sed 's/ /,/g'> $snp_folder/mito_mt_seq_27_major_alleles_sig_snptable.csv
#
# create .fasta major alleles (consensus)
awk 'NR>1' $snp_folder/mito_mt_seq_27_major_alleles_sig_snptable.csv |
awk -F, '{fas=">"$1"\n"; for(i=2;i<=NF;i++){fas=fas$i}; print fas}' \
> $snp_folder/mito_mt_seq_27_major_alleles_sig.fas


#------------------
#- SNP annotation -
#------------------

gunzip $snp_folder/mito_mt_seq_27.vcf.gz
#extract sig SNPs
awk '/^#/ ' $snp_folder/mito_mt_seq_27.vcf > $snp_folder/mito_mt_seq_27_sig.vcf
mapfile -t a < $snp_folder/sig_pos_tmp
for i in ${a[@]}; do awk -v i="$i" '($2==i)' $snp_folder/mito_mt_seq_27.vcf; done \
>> $snp_folder/mito_mt_seq_27_sig.vcf

# annotate SNPs
java -jar scripts/snpEff/snpEff.jar download -v Drosophila_melanogaster
java -Xmx8g -jar scripts/snpEff/snpEff.jar -ud 0 -no-downstream -no-upstream -no-utr Drosophila_melanogaster -stats \
$snp_folder/mito_mt_seq_sig.html $snp_folder/mito_mt_seq_27_sig.vcf > $snp_folder/mito_mt_seq_27_sig_annotated.vcf

rm $snp_folder/*tmp*
rm -r $snp_folder/mito_mt_seq_27