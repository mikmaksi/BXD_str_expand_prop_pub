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
	qtl_res_dir = '../../../data/gene_expr/qtlres'
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
	gene_expr = qtl_coloc_res %>%
		select(GN, gene_expr_qtl) %>%
		unnest(cols = gene_expr_qtl) %>%
		rename(probe = ProbeSet) %>%
		inner_join(probe_info, by = 'probe') %>%
		group_by(GN, gene_name) %>%
		mutate(thresh = mean(thresh)) %>%
		ungroup %>%
		select(!pval)
		# nest(data = c(marker, LOD))

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

# load expression values
{
	# read expression values
	expr_vals = readRDS('../../../data/gene_expr/expr_vals.rds')

	# subset expr vals to those in gene expr QTL results
	expr_vals = expr_vals %>%
		select(GN, probe = ProbeSet, strain, expr_val) %>%
		inner_join(probe_info, by = 'probe')

	# check
	# expr_vals %>% select(matches('gene_'), matches('probe_'), in_range) %>% distinct %>% skimr::skim()

	# filter these for genes of interest as well
	expr_vals = expr_vals %>%
		semi_join(gene_expr %>% distinct(gene_name), by = 'gene_name')
}

# reduce data
{
	# join marker position and sort gene_expr
	gene_expr = gene_expr %>% 
		left_join(marker_pos, by = 'marker') %>%
		arrange(desc(LOD))

	# NOTE: only chr13 analyzed
	# gene_expr %>% distinct(mark_chr)

	if (0) {
		### NOTE: best point and trace calculation have now moved into analysis section ###
		# easier to trace source of data and gene_expr_db is not that large

		# best point
		if (0) {
			# best point within QTL region
			# this unfairly favors genes that are more centrally located within QTL window
			best_point = gene_expr %>%
				filter(mark_pos >= ci_lo, mark_pos <= ci_hi) %>%
				distinct(GN, gene_name, mark_chr, .keep_all = TRUE)
				# best_point %>% object.size %>% print(unit = 'Mb')
				# best_point %>% skimr::skim(mark_pos)
		} else {
			# best point within QTL sized region around gene
			ci_hwind = (ci_hi - ci_lo)/2
			best_point = gene_expr %>%
				filter(mark_pos >= (gene_end+gene_pos)/2 - ci_hwind, mark_pos <= (gene_end+gene_pos)/2 + ci_hwind) %>%
				distinct(GN, gene_name, mark_chr, .keep_all = TRUE)
		}

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
		best_trace_db$expr_vals = expr_vals %>%
			mutate(GN = fct_relevel(GN, best_trace_db$signal %>% pull(GN) %>% levels),
				   probe = fct_relevel(probe, best_trace_db$probes %>% pull(probe) %>% levels),
				   gene_name = fct_relevel(gene_name, best_trace_db$genes %>% pull(gene_name) %>% levels)
			) %>%
			select(GN, probe, gene_name, strain, expr_val) %>%
			semi_join(best_trace_db$signal %>% distinct(GN, probe, gene_name), by = c('GN', 'probe', 'gene_name'))
			
		# walk(best_trace_db, ~object.size(.x) %>% print(unit = 'Mb'))

		### END OF BEST TRACE ###
	}

	# transform gene_expr into a kind-of database optimized for storage and retrieval
	gene_expr_db = list(
		# convert to factors to save space
		signal = gene_expr %>% 
			mutate(across(c(GN, marker, probe, gene_name, probe_chr), as.factor)) %>%
			select(GN, marker, probe, gene_name, LOD, thresh),
		markers = gene_expr %>% 
			mutate_at('marker', as.factor) %>%
			select(marker, contains('mark_')) %>%
			distinct,
		probes = gene_expr %>% 
			mutate_at('probe', as.factor) %>%
			select(probe, contains('probe_'), BlatSeqLen, is_probeset, n_probes, n_var_per_probe, alig_data) %>%
			distinct,
		genes = gene_expr %>% 
			mutate_at('gene_name', as.factor) %>%
			select(contains('gene_')) %>%
			distinct
	)
	gene_expr_db$expr_vals = expr_vals %>%
		mutate(GN = fct_relevel(GN, gene_expr_db$signal %>% pull(GN) %>% levels),
			   probe = fct_relevel(probe, gene_expr_db$probes %>% pull(probe) %>% levels),
			   gene_name = fct_relevel(gene_name, gene_expr_db$genes %>% pull(gene_name) %>% levels)
		) %>%
		select(GN, probe, gene_name, strain, expr_val)
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

	# check
	# best_trace_db$expr_vals %>% nrow
	# gene_expr_db$expr_vals %>% nrow

	# save objects
	out_dir = '../../../data/gene_expr/qtl_agg'; dir_create(out_dir)
	# saveRDS(best_point   ,  path(out_dir , 'best_point.rds'))
	# saveRDS(best_trace_db,  path(out_dir , 'best_trace_db.rds'))
	saveRDS(gene_expr_db ,  path(out_dir , 'gene_expr_db.rds'))
	saveRDS(perc_expand_db, path(out_dir , 'perc_expand_db.rds'))
}
