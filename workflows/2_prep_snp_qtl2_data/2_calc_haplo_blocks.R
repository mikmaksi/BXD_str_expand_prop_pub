#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: group founder haplotypes into haploblocks
# similar to going from bed to bedgraph file

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(fs)
library(qtl2)
library(tidygraph)

# for parallel
library(furrr)
library(progressr)

# clean variables
rm(list = ls())

# function for determining inheritance blocks using tidygraph
# inherit_smry = inherit_smry %>% filter(lab %in% c('B', 'D', 'H')); connect_depth = 10
calc_inherit_blocks = function(inherit_smry, connect_depth = 1, parallel = FALSE) {

    # loop over strains and chromosomes
    # .x = inherit_smry %>% filter(chr == 'chr1', strain == 'BXD001')

    if (parallel) { plan(multisession, workers = 10) } else { plan(sequential) }
    with_progress({
	p <- progressor(steps = inherit_smry %>% distinct(strain, chr) %>% nrow)
	res = inherit_smry %>% 
	    nest(data = !c(strain, chr)) %>%
	    mutate(smry = future_map(data, function(.x) {
		# progress
		p()

		# add id
		.x = .x %>% 
		    arrange(pos, end) %>%
		    mutate(loc_id = 1:n())

		# make pairs of loci at a certain maximum connection distance
		a = rep(1:max(.x$loc_id), each = connect_depth)
		b = rep(1:connect_depth, times = max(.x$loc_id))
		loc_pairs = tibble(varA = a, varB = a+b) %>%
		    left_join(.x %>% select(loc_id, A = lab), by = c('varA' = 'loc_id')) %>%
		    left_join(.x %>% select(loc_id, B = lab), by = c('varB' = 'loc_id')) %>%
		    mutate(connected = A == B) %>%
		    mutate(connected = replace_na(connected, FALSE))

		# create network graph
		nodes = .x %>% distinct(loc_id)
		edges = loc_pairs %>%
		    filter(connected) %>%
		    select(from = varA, to = varB)    
		ld_graph = tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

		# cluster based on connected components
		ld_graph = ld_graph %>% mutate(clst = group_components(type = 'strong'))
		    
		# check that one locus can not be in multiple clusters
		# ld_graph %>%
		#     as_tibble %>%
		#     group_by(loc_id) %>%
		#     summarise(n_clust = length(unique(clst))) %>%
		#     pull(n_clust) %>% unique

		# get first last variant from every cluster
		blocks = ld_graph %>% 
		    as_tibble %>% 
		    group_by(clst) %>% 
		    summarise(from = min(loc_id), to = max(loc_id),
			      n_loci = length(unique(loc_id))) %>%
		    ungroup

		# join the actual genomic positions
		blocks = blocks %>%
		    left_join(.x %>% select(loc_id, pos, lab), by = c('from' = 'loc_id')) %>%
		    left_join(.x %>% select(loc_id, end), by = c('to' = 'loc_id'))

		# find every pair of overlapping blocks
		# set the cross up so that pos_b > pos_a
		# means that overlapping blocks are always like this:
		#   --------------
		#          -------------
		blocks = blocks %>% arrange(from, to) %>% mutate(clst = 1:n())

		# check that cblocks are arranged left to right
		if (0) {
		    # should have equal number has_olap_ab=has_olap_ba (dups)
		    crossing(
			blocks %>% select(clst_a = clst, pos_a = pos, end_a = end),
			blocks %>% select(clst_b = clst, pos_b = pos, end_b = end) 
		    ) %>%
			filter(clst_a != clst_b) %>%
			mutate(has_olap_ab = pos_b < end_a & pos_b > pos_a) %>%
			mutate(has_olap_ba = pos_a < end_b & pos_a > pos_b) %>%
			count(has_olap_ab, has_olap_ba)

		    # should be all has_olap_ab
		    crossing(
			blocks %>% select(clst_a = clst, pos_a = pos, end_a = end),
			blocks %>% select(clst_b = clst, pos_b = pos, end_b = end) 
		    ) %>%
			filter(clst_a < clst_b) %>%
			mutate(has_olap_ab = pos_b < end_a & pos_b > pos_a) %>%
			mutate(has_olap_ba = pos_a < end_b & pos_a > pos_b) %>%
			count(has_olap_ab, has_olap_ba)

		    # should be all has_olap_ba
		    crossing(
			blocks %>% select(clst_a = clst, pos_a = pos, end_a = end),
			blocks %>% select(clst_b = clst, pos_b = pos, end_b = end) 
		    ) %>%
			filter(clst_a > clst_b) %>%
			mutate(has_olap_ab = pos_b < end_a & pos_b > pos_a) %>%
			mutate(has_olap_ba = pos_a < end_b & pos_a > pos_b) %>%
			count(has_olap_ab, has_olap_ba)
		}

		#
		block_olaps = crossing(
		    blocks %>% select(clst_a = clst, pos_a = pos, end_a = end),
		    blocks %>% select(clst_b = clst, pos_b = pos, end_b = end) 
		) %>%
		    filter(clst_a < clst_b) %>%
		    mutate(has_olap = pos_b < end_a & pos_b > pos_a) %>%
		    filter(has_olap) %>%
		    mutate(olap_bp = end_a - pos_b + 1,
			    size_b = end_b - pos_b + 1) %>%
		    mutate(olap_bp   = if_else(size_b <= olap_bp, size_b, olap_bp)) %>%
		    mutate(olap_type = if_else(olap_bp == size_b, 'full', 'partial')) %>%
		    select(matches('clst'), pos_a, end_a, pos_b, end_b, olap_bp, olap_type)

		# check
		# block_olaps %>% count(olap_type)

		# output
		return(list(blocks = blocks, block_olaps = block_olaps))
	    }))
    })

    # extract results from list
    res = res %>%
	mutate(blocks = map(smry, 'blocks'),
	       block_olaps = map(smry, 'block_olaps'))

    # separate data.frames from output
    blocks      = res %>% select(chr, strain, blocks) %>% unnest(blocks)
    block_olaps = res %>% select(chr, strain, block_olaps) %>% unnest(block_olaps)

    # output
    return(list(blocks = blocks, block_olaps = block_olaps))
}

# load the cross
# bxd_cross = read_cross2('bxd.json')

# config
data_dir = '../../data/snp_qtl2/gw'

# load physical map
phys_map = readRDS(path(data_dir, 'pmap.rds')) %>% 
    map_df(~tibble(marker = names(.x), pos = .x*1e6 %>% as.integer), .id = 'chr') %>%
    mutate(chr = str_c('chr', chr))

# load genotype probabilities
snp_probs = readRDS(path(data_dir, 'probs.rds'))

# set to minprob to 0.5 to get return a haplotype for every locus/strain combination 
snp_haplo = maxmarg(snp_probs, minprob = 0.5)
snp_haplo = map_df(snp_haplo, 
   ~.x %>% 
       as.data.frame %>% 
       rownames_to_column(var = 'strain') %>%
       pivot_longer(cols = !strain, names_to = 'marker', values_to = 'gt')) %>%
    mutate(gt = recode(gt, `1` = 'B', `2` = 'D'))

# make haplotype blocks
inherit_smry   = snp_haplo %>% left_join(phys_map, by = 'marker') %>% rename(lab = gt) %>% mutate(end = pos)
inherit_blocks = calc_inherit_blocks(inherit_smry, connect_depth = 1, parallel = TRUE)
blocks         = inherit_blocks$blocks

# save blocks
blocks %>% 
    select(strain, chr, pos, end, lab, n_loci) %>%
    write_tsv('../../data/inherit_smry/snp_haplo_blocks.tsv')
