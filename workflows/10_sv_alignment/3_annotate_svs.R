#!/home/momaksimov/anaconda3/envs/r4/bin/Rscript

# about: annotated structural variants using vep
# use previse INS/DEL coordinates instead of sequences in vcf

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(cowplot)
library(fs)

# config
out_dir = '../../data/sv_alig'
sv_data_file = path(out_dir, 'sv_data.rds')
out_dir = '../../data/sv_alig'; dir_create(out_dir)

# read the structural variants data
sv_data = readRDS(sv_data_file)

# make a temporary directory
tmp_dir = path(out_dir, 'vep'); dir_create(tmp_dir)

# write a vep formatted file
vep_files = list(
	in_file      = path(tmp_dir, 'vep_in.vcf'),
	out_file     = path(tmp_dir, 'vep_out.tsv'),
	stats_file   = path(tmp_dir, 'stats_out.tsv'),
	run_file     = 'run_vep.sh'
)

# format the data into default vep format and write to temp file
sv_regions = sv_data$locs_data %>%
	mutate(strand = '+') %>%
	unite('locus', c('chr', 'pos', 'end', 'sv_type'), remove = FALSE) %>%
	select(chr, pos, end, sv_type, strand, locus) %>%
	distinct %>%
	arrange(chr, pos, end)
# sv_regions %>% distinct(chr, pos, end)
write_tsv(sv_regions, vep_files$in_file, col_names = FALSE)

# set up vep run
vep_cache_dir = '/projects/ps-gymreklab/resources/dbase/mouse/vep/'
cmd = sprintf('vep -v --format ensembl --cache --offline --dir %s --species mus_musculus --assembly GRCm38 --cache_version 102 -i %s -o %s --stats_file %s --tab --everything --force_overwrite', 
		  vep_cache_dir, vep_files$in_file, vep_files$out_file, vep_files$stats_file)
write_lines(c('#!/bin/bash', '', cmd), vep_files$run_file)
# NOTE: have to run manually because it is in a different conda environment

# remove temporary files
file_delete(vep_files$in_file)

# load annotation after running vep
vep_annot = read_tsv(vep_files$out_file, comment = '##')

# check set membership
to_comp = list(
	req = sv_regions %>% select(locus) %>% unique,
	ret = vep_annot %>% select(locus = '#Uploaded_variation') %>% distinct
)
all_equal(intersect(to_comp$req, to_comp$ret), union(to_comp$req, to_comp$ret))
