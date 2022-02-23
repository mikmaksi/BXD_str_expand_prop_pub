#!/bin/bash

# about: get founder genotypes for a list of all snps where 
#   - founders are not identical
#   - neither founder is missing
#   - neither founder is het
# output is headerless file with chr, pos, end, C57 genotype, DBA genotype

# config
vcf_dir='/projects/ps-gymreklab/resources/datasets/BXD/david_dropbox'
vcf_file='Merged_gvcf_files_all_chr_screen_recalibrated_INDEL_variants_99.9_PASSED_variants.recode_in_at_least_20.vcf.gz'
out_dir='../../data/snp_gts/' && mkdir -p $out_dir
out_file="${out_dir}/nohetfou_nomiss_noident.tsv.gz"

# bcftools view -h ${vcf_dir}/${vcf_file} | tail -n1 | datamash transpose

# run bcftools
# NOTE: bcftools allows for filtering of hets, but it is awkward
bcftools query -f '%CHROM\t%POS\t%POS[\t%GT]\n' \
    --include 'TYPE="snp"' \
    -S samples.txt \
    -H \
    ${vcf_dir}/${vcf_file} | \
    awk '{if ($4==$5) next;
          if ($4=="./.") next;
	  if ($5=="./.") next;
	  split($4, a, "/"); 
	  split($5, b, "/");
	  if (a[1] != a[2]) next;
	  if (b[1] != b[2]) next;
	  print $0
    }' OFS='\t' | gzip -c > $out_file
