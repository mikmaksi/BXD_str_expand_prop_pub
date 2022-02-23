#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: smaller manageable vcfs and bams, which can be downloaded locally for use with IGV/JBrowse

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(cowplot)
library(fs)
library(qtl2)
devtools::load_all('../../BXDstrs_package/BXDstrs') # to get vcf files

# config
out_dir = '../../data/vep_annot_new_wind'; dir_create(out_dir)

# regions of interest
peak_chr = 'chr13'
roi_st = 83.78112
roi_en = 93.41913
region = sprintf('%s:%s-%s', peak_chr, roi_st*1e6, roi_en*1e6) 

# loop over vcf files
to_proc = BXDstrs::vcf_files[c('bxd_strs', 'bxd_snp_indel_unfilt')]
to_proc$bxd_svs = path('/projects/ps-gymreklab/resources/datasets/BXD/pangenome/symlinks/chr13.pang.10mb.large.vcf.gz')
to_proc = tibble(
	vcf_name = names(to_proc), 
	vcf_file = as.character(to_proc), 
	out_file = c('bxd_strs', 'bxd_snp_indel', 'bxd_svs'), 
	loc_type = c('str', 'snp', 'sv')
)
pwalk(to_proc, function(vcf_file, vcf_name, out_file, loc_type) {
	out_file = path(out_dir, out_file, ext = 'annot.vcf.gz')
	if (loc_type == 'str') {
		cmd = sprintf('bcftools view %s -r %s | bcftools annotate - --set-id \'%%CHROM\\_%%POS\\_%%INFO/END\' -o %s -O z; tabix -f %s', 
				  vcf_file,
				  region,
				  out_file,
				  out_file)
	} else if (loc_type == 'snp') {
		cmd = sprintf('bcftools view %s -r %s | bcftools annotate - --set-id \'%%CHROM\\_%%POS\\_%%POS\' -o %s -O z; tabix -f %s', 
				  vcf_file,
				  region,
				  out_file,
				  out_file)
	} else if (loc_type == 'sv') {
		cmd = sprintf('bcftools view %s -r %s | bcftools annotate - --set-id \'%%CHROM\\_%%POS\\_%%END\' -o %s -O z; tabix -f %s', 
				  vcf_file,
				  region,
				  out_file,
				  out_file)
	}
	print(cmd)
	system(cmd)
})
