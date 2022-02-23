#!/home/momaksimov/anaconda3/envs/r4/bin/Rscript

# about: pull actual probe sequences and whether the probe is a single probe or a ProbeSet

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(cowplot)
library(fs)
library(progress)

# get list of GN dataset files
gn_files = dir_info('/projects/ps-gymreklab/mikhail/070521_expanded_GN_dbase/data/GN_express/no_crlf_line_fix', glob = '*.txt')
gn_files = gn_files %>%
	select(fpath = path) %>%
	mutate(fname = fpath %>% path_file) %>%
	filter(str_detect(fname, '^GN\\d+_')) %>%
	mutate(GN = fname %>% str_replace('^(GN\\d+)[_-].*', '\\1'))

# define columns of interest
# cols_of_int = c('ProbeID', 'ProbeSetID', 'ProbeSet', 
# 				'Gene.Symbol', 'symbol', 'Gene Symbol',
# 				'BlatSeq', 'Blat Sequence', 
# 				'Blat.Mb.Start', 'Blat.Mb.End', 'Blat Mb Start', 'Blat Mb End')
cols_of_int = c('ProbeId' = 'ProbeId', 
				'ProbeId' = 'ProbeID', 
				'ProbeSet' = 'ProbeSetID', 
				'ProbeSet' = 'ProbeSet', 

				'Symbol' = 'Gene.Symbol', 
				'Symbol' = 'symbol', 
				'Symbol' = 'Gene Symbol',

				'BlatSeq' = 'BlatSeq', 
				'BlatSeq' = 'Blat Sequence', 

				'pos' = 'Blat.Mb.Start', 
				'pos' = 'Blat Mb Start', 
				'pos' = 'Probe_set_Blat_Mb_start',

				'end' = 'Blat Mb End',
				'end' = 'Probe_set_Blat_Mb_end', 
				'end' = 'Blat.Mb.End')

# read every file and get just one row to get columns
pb <- progress_bar$new(total = nrow(gn_files))
gn_files = gn_files %>%
	mutate(data = map(fpath, function(.x) {
		pb$tick()
		.x %>% path_file %>% print

		# for checking column names
		.x = read_tsv(.x, skip = 33, n_max = 1, col_types = cols(.default = 'c'))
		# print(.x)

		# detect ProbeId (these need to be concatenated with dset)
		if ('ProbeId' %in% names(.x)) {
			.x = .x %>% mutate(ProbeId = str_c('gn_', ProbeId))
		}

		# subselect columns of interest
		.x = .x %>% select(any_of(cols_of_int))

		# add a label if probeset
		if ('ProbeSet' %in% names(.x)) {
			.x = .x %>% mutate(is_probeset = TRUE) %>% rename(ProbeId = ProbeSet)
		} else {
			.x = .x %>% mutate(is_probeset = FALSE)
		}
		# print(.x)

		# output
		return(.x)
	}))


# check if every dataset has every columns name
gn_files %>%
	mutate(cols = map(data, ~names(.x))) %>%
	select(GN, cols) %>%
	unnest(cols) %>%
	mutate(present = TRUE) %>%
	complete(GN, cols, fill = list(present = FALSE)) %>%
	filter(!present) %>%
	pull(GN) %>% unique

# set an output directory
out_dir = '../../data/probe_blat_seqs'; dir_create(out_dir)
by_dset_dir = path(out_dir, 'by_dset'); dir_create(by_dset_dir)

# loop over each file and write the data to file
# redo = c("GN152", "GN157", "GN298", "GN282", "GN284", "GN285", "GN208", "GN209", "GN260", "GN149", "GN145", "GN178", "GN179", "GN204", "GN205", "GN146", "GN150", "GN148", "GN144", "GN151", "GN147") 
pb <- progress_bar$new(total = min(redo %>% length, nrow(gn_files)))
pwalk(gn_files %>% select(fpath, GN), function(fpath, GN) {
	pb$tick()
	fpath %>% path_file %>% print

	# for checking column names
	.x = read_tsv(fpath, skip = 33, n_max = Inf, col_types = cols(.default = 'c'))
		
	# detect ProbeId (these need to be concatenated with dset)
	if ('ProbeId' %in% names(.x)) {
		.x = .x %>% mutate(ProbeId = str_c(1:n(), '_', GN))
	}

	# subselect columns of interest
	.x = .x %>% select(any_of(cols_of_int))

	# add a label if probeset
	if ('ProbeSet' %in% names(.x)) {
		.x = .x %>% mutate(is_probeset = TRUE) %>% rename(ProbeId = ProbeSet)
	} else {
		.x = .x %>% mutate(is_probeset = FALSE)
	}

	# save
	saveRDS(.x, path(by_dset_dir, str_c(GN, '.rds')))
})

# check output files
# dir_ls(by_dset_dir) %>% length
# readRDS(path(by_dset_dir, 'GN129.rds'))

# try loading everything
all_probe_blats = dir_ls(by_dset_dir) %>% as.character %>% set_names(. %>% path_file %>% path_ext_remove) %>% map_df(~readRDS(.x), .id = 'GN')

# check on GN152
# all_probe_blats %>% filter(GN == 'GN152')

# make list of unqiue dataset/probe ids
dset_probe = all_probe_blats %>% distinct(GN, ProbeId)

# take distinct ... becomes much smaller
all_probe_blats = all_probe_blats %>% select(!GN) %>% distinct

# convert pos and end to integer
all_probe_blats = all_probe_blats %>%
	mutate(across(c(pos, end), ~as.numeric(.x)*1e6 %>% as.integer))

# nearly all can be converted so this is perfect
all_probe_blats %>% count(is.na(pos)) %>% mutate(prop = n/sum(n))
# A tibble: 2 x 3
# `is.na(pos)`       n   prop
# FALSE        3030923 0.968 
# TRUE           98868 0.0316

# save the combined object
saveRDS(all_probe_blats, path(out_dir, 'probe_blat_seqs.rds'))

# save list of probes/datasets
dset_probe %>% mutate(across(everything(), as.factor)) %>%
	# object.size %>% print(unit = 'Mb') %>%
	saveRDS(path(out_dir, 'dset_probe.rds'))
