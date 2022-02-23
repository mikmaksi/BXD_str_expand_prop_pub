#!/home/momaksimov/anaconda3/envs/r4/bin/Rscript

# about:

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(cowplot)
library(fs)

# config
cc_info_file = '../../info/cc_sample_name_to_bamfile'
out_dir = '../../data/colab_cross/coverage'; dir_create(out_dir)

# read info on CC bam files
strain_info = read_tsv(cc_info_file, 
					   col_names = c('bam_file', 'strain', 'bam_dir'), 
					   col_types = cols(bam_file = 'c', strain = 'c', bam_dir = 'c'))

# add C57 to the list of founders
strain_info = strain_info %>%
	mutate(base_bam_dir = '/projects/ps-gymreklab/resources/datasets') %>%
	add_row(bam_file = '4512-JFI-0333_C57BL_6J_two_lanes_phased_possorted_bam.bam',
			base_bam_dir = '/projects/ps-gymreklab/resources/datasets/BXD',
			bam_dir = 'drive3',
			strain = 'C57BL') %>%
	mutate(strain_type = if_else(str_detect(strain, '^CC0'), 'cc_strain', 'founder')) %>%
	mutate(path = path(base_bam_dir, bam_dir, bam_file))

# check
# strain_info %>% count(strain_type)

# define function for querying coverage
# step_size = 100
# step_size = 1
get_read_depth = function(bam_files, step_size, out_file, coor) {
	# make windows
	# chr13:92,340,715-92,358,532
	region = sprintf('%s:%d-%d', coor$chr, coor$st, coor$en)
	winds = tibble(x = c(coor$st, coor$en) %>% 
				   cut_width(width = step_size, dig.lab = 10, boundary = coor$st) %>% 
				   levels) %>%
		separate('x', c('from', 'to'), ',') %>%
		mutate(across(c(from, to), ~str_replace(.x, '[\\(\\]\\[]', ''))) %>%
		mutate(chr = coor$chr)
		write_tsv(winds %>% select(chr, from, to), '_tmp_regions.bed', col_names = FALSE)

	# query coverage per region window
	# .x = bam_files %>% slice(1) %>% pull(path)
	cov_data = bam_files %>%
		mutate(read_cov = map(path, function(.x) {
			# user out
			print(.x)

			# run bedtools coverage
			# much faster to samtools view a small subregion first
			# .x = 'ccri_sym_links/CC001_Unc_M4363_NYGC.bam' 
			cmd = sprintf("samtools view -b %s %s | bedtools coverage -a _tmp_regions.bed -b -", .x, region)
			cov_data = pipe(cmd) %>% 
				read_tsv(col_names = c('chr', 'pos', 'end', 'n', 'nonzero_bp', 'len', 'nonzero_frac'),
						 col_types = cols(pos = 'i', end = 'i', n = 'i', nonzero_bp = 'i', len = 'i', nonzero_frac = 'd'))
			# debug
			# bedtools coverage -a _tmp_regions.bed -b ccri_sym_links/CC001_Unc_M4363_NYGC.bam
			# bedtools coverage -a _tmp_regions.bed -b /projects/ps-gymreklab/resources/datasets/ColabCross/JAX_bam/CC001_Unc_M4363_NYGC.sorted.bam
			# samtools view /projects/ps-gymreklab/resources/datasets/ColabCross/JAX_bam/CC001_Unc_M4363_NYGC.sorted.bam chr13:92340715-92358532 | less -S
			# samtools view -b /projects/ps-gymreklab/resources/datasets/ColabCross/JAX_bam/CC001_Unc_M4363_NYGC.sorted.bam chr13:92340715-92358532 | bedtools coverage -a _tmp_regions.bed -b - | head

			# output
			cov_data %>% select(chr, pos, end, n)
		}))
	if (file_exists('_tmp_regions.bed')) file_delete('_tmp_regions.bed')

	# collate
	cov_data = cov_data %>%
	select(strain, read_cov) %>%
	unnest(read_cov)

	# save coverage
	saveRDS(cov_data, out_file)
}

### CC strains

# high res (single bp resolution), low res would be 10 bp
coor = list(chr = 'chr13', st = 92340715-1e3, en = 92358532+1e3)
get_read_depth(strain_info %>% filter(strain_type == 'cc_strain'), 
			   step_size = 1, 
			   out_file = path(out_dir, 'coverage_highres.rds'), 
			   coor = coor)

# read the results
cov_data = readRDS(path(out_dir, 'coverage_highres.rds'))

# remove extra strains and format
strains_exclude = c(
	'CC071_TauUnc_M3425_NYGC', 
	'CC035_Unc_M1489_NYGC', 
	'C040_TauUncJ_F1000_NYGC', 
	'CC041_TauUncJ_F1000_NYGC',
	'CC018_UncJ_F1000_NYGC',
	'CC040_TauUncJ_F1000_NYGC',
	'CC019_TauUncJ_F1000_NYGC'
)
cov_data_filt = cov_data %>%
    filter(!strain %in% strains_exclude) %>%
    # mutate(across(c(pos, end), ~.x/1e6)) %>%
    # mutate(Mb = (pos+end)/2) %>%
    separate('strain', 'strain', sep = '_', extra = 'drop')

# check object size
# cov_data %>% object.size %>% print(unit = 'Mb')
# cov_data_filt %>% object.size %>% print(unit = 'Mb')

### CC founders
coor = list(chr = '13', st = 92340715-1e3, en = 92358532+1e3)
get_read_depth(strain_info %>% filter(strain_type == 'founder') %>% filter(strain != 'C57BL'), 
			   step_size = 1, 
			   out_file = path(out_dir, 'coverage_highres_fou.rds'), 
			   coor = coor)

### C57BL
coor = list(chr = 'chr13', st = 92340715-1e3, en = 92358532+1e3)
get_read_depth(strain_info %>% filter(strain_type == 'founder') %>% filter(strain == 'C57BL'), 
			   step_size = 1, 
			   out_file = path(out_dir, 'coverage_highres_c57.rds'), 
			   coor = coor)

# read the results
cov_data_fou = list(fou = path(out_dir, 'coverage_highres_fou.rds'), 
	 c57 = path(out_dir, 'coverage_highres_c57.rds')) %>% map(readRDS)
cov_data_fou$fou = cov_data_fou$fou %>% mutate(chr = str_c('chr', chr))
cov_data_fou = map_df(cov_data_fou, ~.x)

# combine founder and cc strain data
cov_data_all = bind_rows(cov_data_filt, cov_data_fou)

# save
saveRDS(cov_data_all, path(out_dir, 'coverage_highres_filt.rds'))
