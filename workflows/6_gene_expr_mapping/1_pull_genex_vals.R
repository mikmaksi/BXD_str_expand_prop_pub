#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: query gene expression data for probes with a define region

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(fs)
library(DBI)
library(dbplyr)

# function for querying gene expr from databases
query_expr = function(expr_db_file, gene_names = NULL, gene_range = NULL) {
	# query either by gene name or by a genomic range
	# gene_names is a vector of gene names
	# gene_range is a named list with chr, pos and end elements

	# connect to expression database
	conn = dbConnect(RSQLite::SQLite(), expr_db_file)

	# query probes per dataset
	probes_per_dset = tbl(conn, sql('select dset_id, probe_id from probes_per_dset')) %>% collect()
	dset_names = tbl(conn, sql('select id as dset_id, GN from dset_names')) %>% collect()
	strains = tbl(conn, 'strains') %>% collect()

	# join the dataset name
	probes_per_dset = probes_per_dset %>% left_join(dset_names, by = 'dset_id')

	# user out
	print('Getting probe information')
	if (!is.null(gene_names) & is.null(gene_range)) {
		# query probe information for a list of genes
		probe_info = tbl(conn, 
				 sql(sprintf("select * from probe_info
						 where \"Gene.Symbol\" in ('%s')", 
						 str_c(gene_names, collapse = "', '"))
				 )) %>% collect()
	} else if (is.null(gene_names) & !is.null(gene_range)) {
	# query probe information for all probes within a certain window
		probe_info = tbl(conn, 
				 sql(sprintf("select * from probe_info
						 where chr = '%s' and
						 pos_mm10 >= %d and
						 end_mm10 <= %d", 
						 gene_range$chr,
						 gene_range$pos,
						 gene_range$end)
				 )) %>% collect()
	}

	# subselect relevant columns
	probe_info = probe_info %>%
		select(probe_id, ProbeSet, Gene.Symbol, chr_mm10, pos_mm10, end_mm10, strand = Strand.Gene, 
			   Description, Aliases, UniGeneId, OMIM, HomoloGeneID, TargetId)

	# join dataset info to probe_info
	probe_info = probe_info %>% left_join(probes_per_dset, by = 'probe_id')

	# gene expression data for separate anaysis 
	print('Getting expression values')
	tmp_probes = probe_info %>% distinct(probe_id, dset_id)
	dbWriteTable(conn = conn, name = 'tmp_probes', value = tmp_probes, overwrite = TRUE, temporary = TRUE) 
	expr_vals = tbl(conn, sql("select * from expr_vals inner join tmp_probes using ('probe_id', 'dset_id')")) %>% collect()

	# join other information
	expr_vals = expr_vals %>%
		left_join(strains, by = 'strain_id') %>%
		left_join(probe_info %>% select(probe_id, dset_id, ProbeSet, GN), by = c('probe_id', 'dset_id'))
	 
	# close connection to database
	dbDisconnect(conn)

	# output
	return(list(probe_info = probe_info, expr_vals = expr_vals))
}

# config
data_dir = '../../data/gene_expr'; dir_create(data_dir)
expr_db_dir = '../../../070521_expanded_GN_dbase'
out_pref = ''

# GN206 is in a separate database due to its size
expr_dbs = list(GN206 = path(expr_db_dir, 'data/db_files/gn206_db.sqlite'),
				all_other = path(expr_db_dir, 'data/db_files/bxd_expr_db.sqlite'))

# define the region of interest based on 1.5 LOD drop-off
ci_lo = 83.78112
ci_hi = 93.41913
gene_range = list(
	chr = 'chr13',
	pos = ci_lo*1e6,
	end = ci_hi*1e6)

# query probe information and expression values by gene name and by gene range
probe_data = map(expr_dbs, ~query_expr(.x, gene_range = gene_range))

# consolidate ids across databases
probe_data$GN206 = map(probe_data$GN206, function(.x) {
	.x %>%
		left_join(probe_data$all_other$probe_info %>% select(probe_id, ProbeSet), 
			  by = 'ProbeSet', suffix = c('.rem', '.keep')) %>%
		select(!probe_id.rem) %>%
		rename(probe_id = probe_id.keep)
})

# combine
expr_vals  = map_df(probe_data, 'expr_vals') %>% distinct
probe_info = map_df(probe_data, 'probe_info') %>% distinct

# save expression values
saveRDS(expr_vals, path(data_dir, str_c(out_pref, 'expr_vals'), ext = 'rds'))

# save probe_info
probe_info %>%
	rename(chr = chr_mm10, pos = pos_mm10, end = end_mm10) %>%
	saveRDS(path(data_dir, str_c(out_pref, 'probe_info'), ext = 'rds'))
