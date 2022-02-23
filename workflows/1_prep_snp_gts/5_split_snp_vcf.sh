#!/bin/bash
#PBS -q hotel
#PBS -N my_script
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:45:00
#PBS -o job.o
#PBS -e job.e 
#PBS -V
#PBS -M mikmaksi@gmail.com
#PBS -m abe
#PBS -t 1-19
#PBS -A gymreklab-group

# about: split the full snp vcf into per-chromosome files, while subsetting for snps of 
# interest (see $regions_file) and fixing * alt alleles to <*> for PLINK

# cd into main analysis directory
cd /projects/ps-gymreklab/mikhail/072321_bxd_mutator_paper/workflows/1_prep_snp_gts

# debug
# PBS_ARRAYID=1

# config
vcf_dir='/projects/ps-gymreklab/resources/datasets/BXD/david_dropbox'
vcf_file='Merged_gvcf_files_all_chr_screen_recalibrated_INDEL_variants_99.9_PASSED_variants.recode_in_at_least_20.vcf.gz'
regions_file='../../data/snp_gts/nohetfou_nomiss_noident_nosegdup.tsv.gz'

# set up outputs
out_dir='../../data/snp_gts/vcf' && mkdir -p $out_dir
chrom=chr${PBS_ARRAYID}

# get nlines
nlines=`zcat $regions_file | awk -v chr=$chrom '$1==chr {print $1,$2,$3}' OFS='\t' | wc -l`
# echo $nlines

# convert the * alt allele to <*> so that PLINK can understand it
bcftools view ${vcf_dir}/${vcf_file} \
    --include 'TYPE="snp"' \
    -R <(zcat $regions_file | awk -v chr=$chrom '$1==chr {print $1,$2,$3}' OFS='\t') | \
    sed 's/*/<*>/g' | \
    pv -pleIt -s $nlines | \
    bcftools sort - -o ${out_dir}/${chrom}.vcf.gz -O z
tabix -p vcf ${out_dir}/${chrom}.vcf.gz
