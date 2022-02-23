#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: raw QTL traces are too cumbersome to deal with when working with many datasets and many genes
# especially since many probes may be associated with each gene
# it often it makes sense to work with a reduced representation
#   - best QTL trace per gene/dataset combination
#   - peak QTL (single value) per gene/dataset combination

# options
options(stringsAsFactors = FALSE)

# libraries
library(tidyverse)
# library(cowplot)
library(fs)
library(jsonlite)

# config
ci_lo = 83.78112
ci_hi = 93.41913

# load probe to gene information
{
	# load probe to gene information (have gene position and fixed missing gene names for probes)
	probe_info = readRDS('../../../data/gene_expr/probe_to_gene.rds')

	# check
	# probe_info %>% count(is.na(gene_name))
	# probe_info %>% count(probe_chr, gene_chr)
}

# load QTL mapping results for each gene in each dataset
{
	# qtl coloc results: %expanded and gene expr QTL signals in matched strains bundeled together with coloc analysis
	qtl_res_dir = '../../../data/gene_expr/qtlres_dense'
	qtl_coloc_res = dir_info(qtl_res_dir, glob = '*.rds') %>% 
		select(path) %>%
		mutate(f_name = path_file(path)) %>%
		separate('f_name', c('GN'), '_', extra = 'drop') %>%
		pull(path) %>%
		map_df(readRDS)

	# get results for %expanded qtl mapping
	perc_expand = qtl_coloc_res %>%
		select(GN, perc_expand_qtl) %>%
		unnest(cols = perc_expand_qtl) %>%
		distinct() 

	# check size
	# object.size(perc_expand) %>% print(unit = 'Mb')

	# get results for gene expression qtl mapping
	# NOTE: did not calculate permutation-based thresholds with high density snp data
	gene_expr = qtl_coloc_res %>%
		select(GN, gene_expr_qtl) %>%
		unnest(cols = gene_expr_qtl) %>%
		rename(probe = ProbeSet) %>%
		inner_join(probe_info, by = 'probe') %>%
		# group_by(GN, gene_name) %>%
		# mutate(thresh = mean(thresh)) %>%
		# ungroup %>%
		select(!pval)
		# nest(data = c(marker, LOD))

	# get results for colocalization b/w gene eQTL signal and %expanded phenotype QTL
	# NOTE: clean up the upstream script to only save the `summmary` element of the list
	qtl_eqtl_coloc = qtl_coloc_res %>%
		select(GN, coloc_res) %>%
		unnest(cols = coloc_res) %>%
		rename(probe = ProbeSet) %>%
		mutate(data_type = names(coloc_res)) %>%
		filter(data_type == 'result') %>%
		mutate(smry = map(coloc_res, 'summary'),
			   details = map(coloc_res, 'results'),
			   priors = map(coloc_res, 'priors')
		) %>%
		select(GN, probe, smry, details, priors)

	# check
	# gene_expr %>% select(matches('gene_'), matches('probe_')) %>% distinct %>% skimr::skim()

	# check size
	# object.size(gene_expr) %>% print(unit = 'Mb')
	# object.size(perc_expand) %>% print(unit = 'Mb')

	# make sure genes are only on chr13
	gene_expr = gene_expr %>% filter(gene_chr == 'chr13')

	# extract position from marker
	marker_pos = gene_expr %>%
		distinct(marker) %>%
		separate('marker', c('mark_chr', 'mark_pos', NA), '_', convert = TRUE, remove = FALSE) %>%
		mutate(across(mark_pos, ~.x/1e6))

	# check on datasets
	# gene_expr %>% distinct(GN)
}

# reduce data
{
	# shouldn't have more than 1 probe per GN/gene combination; check on this
	if((gene_expr %>% distinct(GN, gene_name, probe) %>% count(GN, gene_name, probe) %>% pull(n) %>% unique) != 1) stop('More than one probe per gene')

	# join marker position and sort gene_expr
	gene_expr = gene_expr %>% 
		left_join(marker_pos, by = 'marker') %>%
		arrange(desc(LOD))

	# best point
	# NOTE: only chr13 analyzed
	# make sure to restrict best point to region of interest and not all of chr13
	best_point = gene_expr %>%
		filter(mark_pos >= ci_lo, mark_pos <= ci_hi) %>%
		distinct(GN, gene_name, mark_chr, .keep_all = TRUE)
		# best_point %>% object.size %>% print(unit = 'Mb')
		# best_point %>% skimr::skim(mark_pos)

	# NOTE: in this case best_trace and gene_expr should be the same
	if (0) {
		# best trace
		# this is a single trace per GN/gene name
		# select top value per GN, gene_name, then use the associate probe to select all value for that probe
		best_trace = gene_expr %>%
			semi_join(best_point, by = c('GN', 'gene_name', 'probe'))

		# transform best_trace into a kind-of database optimized for storage and retrieval
		best_trace_db = list(
			# convert to factors to save space
			signal = best_trace %>% 
			mutate(across(c(GN, marker, probe, gene_name, probe_chr), as.factor)) %>%
			select(GN, marker, probe, gene_name, LOD, thresh),
			markers = best_trace %>% 
			mutate_at('marker', as.factor) %>%
			select(marker, contains('mark_')) %>%
			distinct,
			probes = best_trace %>% 
			mutate_at('probe', as.factor) %>%
			select(probe, contains('probe_')) %>%
			distinct,
			genes = best_trace %>% 
			mutate_at('gene_name', as.factor) %>%
			select(contains('gene_')) %>%
			distinct
		)
		# walk(best_trace_db, ~object.size(.x) %>% print(unit = 'Mb'))
	}

	# transform gene_expr into a kind-of database optimized for storage and retrieval
	gene_expr_db = list(
		# convert to factors to save space
		signal = gene_expr %>% 
			mutate(across(c(GN, marker, probe, gene_name, probe_chr), as.factor)) %>%
			# select(GN, marker, probe, gene_name, LOD, thresh),
			select(GN, marker, probe, gene_name, LOD),
		markers = gene_expr %>% 
			mutate_at('marker', as.factor) %>%
			select(marker, contains('mark_')) %>%
			distinct,
		probes = gene_expr %>% 
			mutate_at('probe', as.factor) %>%
			select(probe, contains('probe_')) %>%
			distinct,
		genes = gene_expr %>% 
			mutate_at('gene_name', as.factor) %>%
			select(contains('gene_')) %>%
			distinct
	)
	# walk(gene_expr_db, ~object.size(.x) %>% print(unit = 'Mb'))

	# transform perc_expand into a kind-of database optimized for storage and retrieval 
	perc_expand_db = list(
		# convert to factors to save space
		signal = perc_expand %>% 
			mutate(across(c(GN, marker), as.factor)) %>%
			select(GN, marker, LOD, thresh),
		markers = perc_expand %>% 
			distinct(marker) %>%
			mutate_at('marker', as.factor) %>%
			separate('marker', c('mark_chr', 'mark_pos'), convert = TRUE, remove = FALSE, extra = 'drop') %>%
			mutate_at('mark_pos', ~.x/1e6)
	)

	# transfor qtl_eqtl_coloc into a kind-of database optimized for storage and retrieval
	coloc_db = list(
		smry = qtl_eqtl_coloc %>% 
			select(GN, probe, smry) %>% 
			mutate(smry = map(smry, ~tibble(var = names(.x), val = .x) %>% spread(var, val))) %>%
			unnest(smry) %>%
			rename_at(vars(matches('^PP')), ~str_replace(.x, 'PP.(.+).abf', '\\1')) %>%
			select(!nsnps),
		details = qtl_eqtl_coloc %>% 
			select(GN, probe, details) %>% 
			unnest(details),
		priors = qtl_eqtl_coloc %>% 
			select(GN, probe, priors) %>% 
			mutate(priors = map(priors, ~tibble(var = names(.x), val = .x) %>% spread(var, val))) %>%
			unnest(priors)
	)

	# save objects
	out_dir = '../../../data/gene_expr/qtl_agg_dense'; dir_create(out_dir)
	# saveRDS(best_point   ,  path(out_dir , 'best_point.rds'))
	# saveRDS(best_trace_db,  path(out_dir , 'best_trace_db.rds'))
	saveRDS(gene_expr_db ,  path(out_dir , 'gene_expr_db.rds'))
	saveRDS(perc_expand_db, path(out_dir , 'perc_expand_db.rds'))
	saveRDS(coloc_db, path(out_dir , 'coloc_db.rds'))
}
