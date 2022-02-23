#!/bin/bash
#PBS -q hotel
#PBS -N my_script
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:30:00
#PBS -o job.o
#PBS -e job.e 
#PBS -V
#PBS -M mikmaksi@gmail.com
#PBS -m abe
#PBS -t 1-19
#PBS -A gymreklab-group

# about: pull fou and ri genotypes for a list snps (see $regions_file)
# output are locus-by-strain matrix; tsv separated; header is commented and needs formatting
# one file per chromosome

# cd into main analysis directory
cd /projects/ps-gymreklab/mikhail/072321_bxd_mutator_paper/workflows/1_prep_snp_gts

# for running locally
# PBS_ARRAYID=$1

# debug
# PBS_ARRAYID=1

# config
vcf_dir='/projects/ps-gymreklab/resources/datasets/BXD/david_dropbox'
vcf_file='Merged_gvcf_files_all_chr_screen_recalibrated_INDEL_variants_99.9_PASSED_variants.recode_in_at_least_20.vcf.gz'
regions_file='../../data/snp_gts/nohetfou_nomiss_noident_nosegdup.tsv.gz'

# set outputs
out_dir='../../data/snp_gts/raw' && mkdir -p $out_dir
out_file="${out_dir}/chr${PBS_ARRAYID}.tsv.gz"

# get nlines
nlines=`zcat $regions_file | awk -v chr=chr${PBS_ARRAYID} '$1==chr {print $1,$2,$3}' OFS='\t' | wc -l`
# echo $nlines

# run bcftools
bcftools query -f '%CHROM\t%POS\t%POS[\t%GT]\n' \
    --include 'TYPE="snp"' \
    -R <(zcat $regions_file | awk -v chr=chr${PBS_ARRAYID} '$1==chr {print $1,$2,$3}' OFS='\t') \
    -H \
    ${vcf_dir}/${vcf_file} | pv -pleIt -s $nlines | gzip -c > $out_file
    # ${vcf_dir}/${vcf_file} | head -n10000 | pv -pl -s 10000 | gzip -c > $out_file
