#!/bin/bash

# about: get all founder and ri genotypes into a raw tsv formatted file
# header commented directly from bcftools

# config
vcf_dir='/projects/ps-gymreklab/mikhail/051820_unified_workflow/data/call_sets/strs'
vcf_file='BXD_variants_MM_strs_highQ_calls.vcf.gz'
out_dir='../../data/str_gts/' && mkdir -p $out_dir
out_file="${out_dir}/all_repcn_raw.tsv.gz"

# bcftools view -h ${vcf_dir}/${vcf_file} | tail -n1 | datamash transpose

# get nlines (do this one time)
# nlines=`bcftools query -f '%CHROM\t%POS\n' ${vcf_dir}/${vcf_file} | wc -l`
nlines=1176016
# echo $nlines

# print everything
bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%REPCN]\n' \
    -H \
    ${vcf_dir}/${vcf_file} | \
    pv -pleIt -s $nlines | gzip -c > $out_file
