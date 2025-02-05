#!/bin/bash

# author: Susanne Voigt

wdir=/mnt/c/Users/ennas/Desktop/mito/mito_seq	#set working directory  
ref_folder=$wdir/mito_seq_ref                   #directory to reference
ref=$ref_folder/dmel628shortheaders.fasta		#reference
bam_folder=$wdir/mito_seq_bam                  #directory to .bam data


#-------------------------------
#- SNP calling and annotation  -
#-------------------------------

# make directory for output
snp_folder=$wdir/mito_seq_snps
mkdir $snp_folder

# get sample names
cd $bam_folder
sample_name=($(ls -1 | grep bam$ | cut -d'm' -f 1)) # sample names - use -d'.' for more general usage
cd ..

# list all bam files 
bam=($(ls -d -1 $bam_folder/* | grep .bam$))

# merge all .bam files in single mpileup file (only retaining nucleotides with BQ >20 and reads with MQ >20) 
samtools mpileup -B -q 20 -Q 20 -f $ref ${bam[@]}  | gzip > $snp_folder/mito_seq.mpileup.gz

# call SNPs with PoolSNP (https://github.com/capoony/PoolSNP -  doi:10.1093/molbev/msaa120)
#-------------------------------------------------------------------------------------------
./scripts/PoolSNP-master/PoolSNP.sh mpileup=$snp_folder/mito_seq.mpileup.gz output=$snp_folder/mito_seq_poolsnp reference=$ref \
names=${!sample_name[@]} \
max-cov=0.99 min-cov=10 min-count=10 min-freq=0.001 miss-frac=0.2 jobs=6 BS=1 
rm -r $snp_folder/mito_seq_poolsnp

# remove SNPs around InDels and in TEs
#-------------------------------------
# identify sites in proximity of InDels with a minimum count of 10 across all samples pooled and mask sites 5bp up- and downstream of InDels 
# (python scripts in the following from https://github.com/capoony/DrosEU_pipeline - doi:10.1093/molbev/msaa120)
python2 scripts/DrosEU_pipeline-master/scripts/DetectIndels.py --mpileup $snp_folder/mito_seq.mpileup.gz --minimum-count 10 --mask 5 | 
gzip > $snp_folder/mito_seq_indel_pos.txt.gz
#
# use Repeatmasker to generate a GFF with location of known TEs
# obtain TE libraries:
cd $ref_folder
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster//dmel_r6.28_FB2019_03/fasta/dmel-all-transposon-r6.28.fasta.gz
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster//dmel_r6.28_FB2019_03/fasta/dmel-all-chromosome-r6.28.fasta.gz
gunzip dmel-all-transposon-r6.28.fasta.gz
gunzip dmel-all-chromosome-r6.28.fasta.gz
# only keep contig name in headers (no spaces):
awk '{print $1}' dmel-all-transposon-r6.28.fasta > dmel-all-transposon-r6.28_fixed-id.fasta
cd ..
# repeat mask D. melanogaster genome using Repeatmasker:
scripts/RepeatMasker/RepeatMasker -pa 1 -e rmblast --lib $ref_folder/dmel-all-transposon-r6.28_fixed-id.fasta \
--gff --qq --no_is --nolow $ref_folder/dmel-all-chromosome-r6.28.fasta
#
# filter SNPs around indels and in TEs from the original VCF obtained from PoolSNP 
python2 scripts/DrosEU_pipeline-master/scripts/FilterPosFromVCF.py --indel $snp_folder/mito_seq_indel_pos.txt.gz \
--te $ref_folder/dmel-all-chromosome-r6.28.fasta.out.gff --vcf $snp_folder/mito_seq_poolsnp.vcf.gz |
gzip > $snp_folder/mito_seq_poolsnp_clean.vcf.gz
rm $ref_folder/dmel-all-transposon-r6.28.fasta
rm $ref_folder/dmel-all-chromosome-r6.28.fasta


#----------------------
#- genetic variation  -
#----------------------

# make directory for output
var_folder=$wdir/mito_seq_var
mkdir $var_folder

# annotate SNPs with SNPeff (doi: 10.4161/fly.19695)
#---------------------------------------------------m
java -jar scripts/snpEff/snpEff.jar download -v Drosophila_melanogaster
java -Xmx8g -jar scripts/snpEff/snpEff.jar -ud 2000 Drosophila_melanogaster -stats $snp_folder/mito_seq_poolsnp.html \
$snp_folder/mito_seq_poolsnp_clean.vcf.gz | gzip > $var_folder/mito_seq_poolsnp_clean_annotated.vcf.gz

 
# convert the VCF to SYNC file format (python scripts in the following from https://github.com/capoony/DrosEU_pipeline - doi:10.1093/molbev/msaa120)
python2 scripts/DrosEU_pipeline-master/scripts/VCF2sync.py --vcf $var_folder/mito_seq_poolsnp_clean_annotated.vcf.gz |
gzip > $snp_folder/mito_seq_snps.sync.gz

# resample SNPS to a 40x coverage for autosomes and 20x for X by random subsampling 
# (script written for pools of males - chr X sampled by half)
# (Watterson’s θ and Tajima’s D are sensitive to coverage variation)
gunzip $snp_folder/mito_seq_snps.sync.gz
python2 scripts/DrosEU_pipeline-master/scripts/SubsampleSync.py --sync $snp_folder/mito_seq_snps.sync --target-cov 40 --min-cov 10 \
> $snp_folder/mito_seq_snps_40x.sync
gzip $snp_folder/mito_seq_snps.sync


# calculate window-wise population genetics statistics Tajima's pi, Watterson's Theta and Tajima's D for 40kB-overlapping 50kB windows
#-------------------------------------------------------------------------------------------------------------------------------------
# ("true" window-sizes based on the number of sites that passed the coverage criteria as calculated from PoolSNP) 
for chromlen in {2L:23513712,2R:25286936,3L:28110227,3R:32079331,X:23542271,4:1348131,Y:3667352}
do
    chrom=$(echo "$chromlen" | cut -d ':' -f1)
    awk -v chrom="$chrom" '$1==chrom' $snp_folder/mito_seq_snps_40x.sync | gzip > $snp_folder/mito_seq_snps_40x_$chrom.sync.gz
    python2 scripts/DrosEU_pipeline-master/scripts/TrueWindows.py --badcov $snp_folder/mito_seq_poolsnp_BS.txt.gz \
    --indel $snp_folder/mito_seq_indel_pos.txt.gz --te $ref_folder/dmel-all-chromosome-r6.28.fasta.out.gff \
    --window 50000 --step 10000 --output truewindows --output $snp_folder/mito_seq_truewin_$chrom \
    -- chromosomes $chromlen
    python2 scripts/DrosEU_pipeline-master/scripts/PoolGen-var.py --input $snp_folder/mito_seq_snps_40x_$chrom.sync.gz --min-count 2  \
    --window 50000 --step 10000 --sitecount $snp_folder/mito_seq_truewin_$chrom-50000-10000.txt --min-sites-frac 0.75 \
    --output $var_folder/mito_seq_var_tmp_$chrom \
    --pool-size 100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100 
done
cat $var_folder/*_50000_10000.pi  > $var_folder/mito_seq_var_50000_10000.pi
cat $var_folder/*_50000_10000.th > $var_folder/mito_seq_var_50000_10000.th
cat $var_folder/*_50000_10000.D  > $var_folder/mito_seq_var_50000_10000.D
rm $var_folder/*tmp*
gzip $snp_folder/mito_seq_snps_40x.sync


# recombination rate per site 
#----------------------------
# back-convert genomic coordinates from R6 (dm6) to R5 (dm3) with liftOver tool from the UCSC Genome Browser (doi: 10.1093/nar/gkac1072)
# for SNPs: 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
wget http://hgdownload.soe.ucsc.edu/goldenPath/dm6/liftOver/dm6ToDm3.over.chain.gz
gunzip $snp_folder/mito_seq_snps.sync.gz 
awk -v OFS='\t' '{print "chr"$1, $2, $2}' $snp_folder/mito_seq_snps.sync | awk '$1 != "chrY"' > $snp_folder/mito_seq_snps_pos_R6.bed # exclude chrom Y 
gzip $snp_folder/mito_seq_snps.sync 
./scripts/liftOver $snp_folder/mito_seq_snps_pos_R6.bed scripts/dm6ToDm3.over.chain.gz \
$snp_folder/mito_seq_snps_pos_R5_tmp.bed $snp_folder/mito_seq_snps_pos_R6toR5_unmapped.txt
grep -v -f $snp_folder/mito_seq_snps_pos_R6toR5_unmapped.txt $snp_folder/mito_seq_snps_pos_R6.bed > $snp_folder/mito_seq_snps_pos_R6.tmp
paste $snp_folder/mito_seq_snps_pos_R5_tmp.bed $snp_folder/mito_seq_snps_pos_R6.tmp > $snp_folder/mito_seq_snps_pos_R5_R6.bed
# for midpoints of window-wise popgen analysis:
cat $var_folder/*_50000_10000.pi | cut -d "." -f 1 | awk -v OFS='\t' '{print "chr"$1, $2, $2}' > $snp_folder/mito_seq_var_mipoints_R6_tmp.bed
./scripts/liftOver $snp_folder/mito_seq_var_mipoints_R6_tmp.bed scripts/dm6ToDm3.over.chain.gz \
$snp_folder/mito_seq_var_mipoints_R5_tmp.bed $snp_folder/mito_seq_var_mipoints_R6_tmp_R6toR5_unmapped.txt
grep -v -f $snp_folder/mito_seq_var_mipoints_R6_tmp_R6toR5_unmapped.txt $snp_folder/mito_seq_var_mipoints_R6_tmp.bed \
> $snp_folder/mito_seq_var_mipoints_R6_tmp.tmp
paste $snp_folder/mito_seq_var_mipoints_R5_tmp.bed $snp_folder/mito_seq_var_mipoints_R6_tmp.tmp \
> $snp_folder/mito_seq_var_midpoints_R5_R6.bed


# recombination rate per site with Drosophila melanogaster Recombination Rate Calculator (http://petrov.stanford.edu/cgi-bin/recombination-rates_updateR5.pl)
# (doi: 10.1016/j.gene.2010.04.015)
cd scripts/RRCv2.3/
# for SNPs:
awk '{print $1":"$2".."$2}' $snp_folder/mito_seq_snps_pos_R5_tmp.bed | cut -d "r" -f 2 > mito_seq_snps_pos_R5_tmp
perl RRC-open-v2.3.pl -M mito_seq_snps_pos_R5_tmp
cut -d "." -f 1  mito_seq_snps_pos_R5_tmp.rrc | sed 's/:/\t/g' > mito_seq_snps_rrc_R5_pos_tmp
grep -f mito_seq_snps_rrc_R5_pos_tmp $snp_folder/mito_seq_snps_pos_R5_R6.bed | cut -f 4,5 | cut -d "r" -f 2,3 \
> mito_seq_snps_rrc_R6_pos_tmp
cut -f 3,6 mito_seq_snps_pos_R5_tmp.rrc > mito_seq_snps_rrc_rates_tmp
paste mito_seq_snps_rrc_R6_pos_tmp mito_seq_snps_rrc_rates_tmp | 
awk 'BEGIN {print "chr\tposR6\tFistonLavier2010\tComeron2012"}{print}' \
> $var_folder/mito_seq_snps_recombrates.txt
# for midpoints of window-wise popgen analysis:
awk '{print $1":"$2".."$2}' $snp_folder/mito_seq_var_mipoints_R5_tmp.bed | cut -d "r" -f 2 > mito_seq_var_mipoints_R5_tmp
perl RRC-open-v2.3.pl -M mito_seq_var_mipoints_R5_tmp
cut -d "." -f 1  mito_seq_var_mipoints_R5_tmp.rrc | sed 's/:/\t/g' > mito_seq_var_mipoints_rrc_R5_pos_tmp
grep -f mito_seq_var_mipoints_rrc_R5_pos_tmp $snp_folder/mito_seq_var_midpoints_R5_R6.bed | cut -f 4,5 | cut -d "r" -f 2,3 \
> mito_seq_var_mipoints_rrc_R6_pos_tmp
cut -f 3,6 mito_seq_var_mipoints_R5_tmp.rrc > mito_seq_var_mipoints_rrc_rates_tmp
paste mito_seq_var_mipoints_rrc_R6_pos_tmp mito_seq_var_mipoints_rrc_rates_tmp | 
awk 'BEGIN {print "chr\tposR6\tFistonLavier2010\tComeron2012"}{print}' \
> $var_folder/mito_var_mipoints_recombrates.txt

rm *tmp*
cd $snp_folder/*tmp*
cd ..


# calculate mean population genetics statistics Tajima's pi, Watterson's Theta and Tajima's D per chrom arm
#-----------------------------------------------------------------------------------------------------------------

# identify SNPs inside introns of < 60bp length (small bug in original python script - line 41 should be input instead of sync)
python2 scripts/DrosEU_pipeline-master/scripts/IntronicSnps_revised.py --gff $ref_folder/dmel-all-filtered-r6.28.gff.gz \
--sync $snp_folder/mito_seq_snps_40x.sync.gz --target-length 60 > $snp_folder/mito_seq_snps_40x_intron60_tmp.sync

# remove SNPs within and in 1Mb distance to chromosomal inversions and with recombination rates <3
python2 scripts/DrosEU_pipeline-master/scripts/FilterByRecomRateNInversion_adapted.py \
--inv scripts/DrosEU_pipeline-master/data/inversions_breakpoints_v5v6.txt \
--RecRates $var_folder/mito_seq_snps_recombrates.txt --input $snp_folder/mito_seq_snps_40x_intron60_tmp.sync \
--D 500000 --r 1 > $snp_folder/mito_seq_snps_40x_intron60.sync

# calculate mean Tajima's pi, Watterson's Theta and Tajima's D per chrom arm
for chromlen in {2L:23513712,2R:25286936,3L:28110227,3R:32079331,X:23542271}
do
    chrom=$(echo "$chromlen" | cut -d ':' -f1) 
    awk -v chrom="$chrom" '$1==chrom' $snp_folder/mito_seq_snps_40x_intron60.sync > $snp_folder/mito_seq_snps_40x_intron60_${chrom}_tmp.sync
    python2 scripts/DrosEU_pipeline-master/scripts/TrueWindows.py --badcov $snp_folder/mito_seq_poolsnp_BS.txt.gz \
    --indel $snp_folder/mito_seq_indel_pos.txt.gz --te $ref_folder/dmel-all-chromosome-r6.28.fasta.out.gff \
    --window 50000000 --step 50000000 --output $snp_folder/mito_seq_intron60_truewin_tmp_$chrom \
    -- chromosomes $chromlen
    python2 scripts/DrosEU_pipeline-master/scripts/PoolGen-var.py --input $snp_folder/mito_seq_snps_40x_intron60_${chrom}_tmp.sync --min-count 2  \
    --window 50000000 --step 50000000 --sitecount $snp_folder/mito_seq_intron60_truewin_tmp_$chrom-50000000-50000000.txt --min-sites-frac 0.75 \
    --output $var_folder/mito_seq_var_intron60_mean_tmp_$chrom \
    --pool-size 100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100 
done

cat $var_folder/*_50000000_50000000.pi  > $var_folder/mito_seq_var_intron60_mean.pi
cat $var_folder/*_50000000_50000000.th  > $var_folder/mito_seq_var_intron60_mean.th
cat $var_folder/*_50000000_50000000.D  > $var_folder/mito_seq_var_intron60_mean.D
rm $var_folder/*tmp*
rm $snp_folder/*tmp*


# estimate major allele frequencies
#----------------------------------
python2 scripts/DrosEU_pipeline-master/scripts/sync2AF.py --inp $snp_folder/mito_seq_snps.sync.gz --out $var_folder/mito_seq_snps_major_allele



# calculate Fst (Weir&Cockerham)
#-------------------------------
# use only SNPs at least 1kb apart (to minimize effects of LD)
gunzip $snp_folder/mito_seq_snps.sync.gz
awk '$1=="2L"' $snp_folder/mito_seq_snps.sync | awk 'f && $2-f<1000 {next} {print $0}; {f=$2}'  > $snp_folder/mito_seq_snps_1kb_apart.tmp
awk '$1=="2R"' $snp_folder/mito_seq_snps.sync | awk 'f && $2-f<1000 {next} {print $0}; {f=$2}'  >> $snp_folder/mito_seq_snps_1kb_apart.tmp
awk '$1=="3L"' $snp_folder/mito_seq_snps.sync | awk 'f && $2-f<1000 {next} {print $0}; {f=$2}'  >> $snp_folder/mito_seq_snps_1kb_apart.tmp
awk '$1=="3R"' $snp_folder/mito_seq_snps.sync | awk 'f && $2-f<1000 {next} {print $0}; {f=$2}'  >> $snp_folder/mito_seq_snps_1kb_apart.tmp
awk '$1=="4"' $snp_folder/mito_seq_snps.sync | awk 'f && $2-f<1000 {next} {print $0}; {f=$2}'  >> $snp_folder/mito_seq_snps_1kb_apart.tmp
awk '$1=="X"' $snp_folder/mito_seq_snps.sync | awk 'f && $2-f<1000 {next} {print $0}; {f=$2}'  >> $snp_folder/mito_seq_snps_1kb_apart.tmp
gzip $snp_folder/mito_seq_snps.sync
# excluding SNPs in regions of no (and NA) recombination rate (Comeron2012) (to minimize effects of LD)
awk '$4>0' $var_folder/mito_seq_snps_recombrates.txt | cut -f 1,2 > $snp_folder/nonzero_recomb_pos.tmp
grep -f $snp_folder/nonzero_recomb_pos.tmp $snp_folder/mito_seq_snps_1kb_apart.tmp > $snp_folder/mito_seq_snps_1kb_apart_nonzero_recomb.sync
rm $snp_folder/*tmp*

# (using python script from https://github.com/capoony/DrosEU_pipeline - doi:10.1093/molbev/msaa120)
# Fst per SNP:
python2 scripts/DrosEU_pipeline-master/scripts/FST.py \
--pool-size 100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100 \
--input $snp_folder/mito_seq_snps.sync.gz --minimum-count 2 --minimum-cov 10 | 
gzip > $var_folder/mito_seq_snps.fst.gz

# average Fst across all SNPs at least 1kb apart & nonzero recombination:
python2 scripts/DrosEU_pipeline-master/scripts/FST.py \
--pool-size 100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100 \
--input $snp_folder/mito_seq_snps_1kb_apart_nonzero_recomb.sync.gz --minimum-count 2 --minimum-cov 10 | 
gzip > $var_folder/mito_seq_snps_1kb_apart_nonzero_recomb.fst.gz
python2 scripts/DrosEU_pipeline-master/scripts/CombineFST.py --diff $var_folder/mito_seq_snps_1kb_apart_nonzero_recomb.fst.gz \
--stat 0 > $var_folder/mito_seq_snps_1kb_apart_nonzero_recomb_average.fst 
rm $var_folder/mito_seq_snps_1kb_apart_nonzero_recomb.fst.gz
