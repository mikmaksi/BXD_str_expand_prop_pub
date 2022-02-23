#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript
#PBS -q hotel
#PBS -N qtl_map_gene_expr
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:30:00
#PBS -o job.o
#PBS -e job.e 
#PBS -V
#PBS -M mikmaksi@gmail.com
#PBS -m abe
#PBS -t 11-241
#PBS -A gymreklab-group

# about: run association testing between snp founder genotypes and gene expression
# run qtl2 in multi-phenotype mode (probe expression = phenotype)

# clean variables
# rm(list = ls())

# options
options(stringsAsFactors = FALSE)

# cd into analysis directory
setwd('/projects/ps-gymreklab/mikhail/072321_bxd_mutator_paper/workflows/6_gene_expr_mapping/map_eqtl_snps_covar')

# get job id
PBS_ARRAYID = as.integer(Sys.getenv('PBS_ARRAYID'))
# PBS_ARRAYID = 1
print(sprintf('Processing job_id %s', PBS_ARRAYID))

# libraries
library(tidyverse)
library(cowplot)
library(fs)
library(qtl2)
library(jsonlite)

# load strain info
strain_info = readRDS('../../../info/strain_info.rds')

# config
# min_strains = 5 # number of strains required to run association testing
n_perm = 100
    
# set up output directory
out_dir = '../../../data/gene_expr/qtlres_snps_covar'; dir_create(out_dir)

# covert a qtl2 type matrix to a tibble
qtl2Totbl = function(data) {
	as.data.frame(data) %>%
		rownames_to_column(var = 'marker') %>%
		as_tibble
}

# estimate pvalues from LOD
lodToPval = function(x) {
	pchisq(x*(2*log(10)), df=1, lower.tail = FALSE)/2
} 

# run scan1 and scan1 coef and merge the data
# gp = this_probs_drop
# ph = this_expr_vals
# ks = this_kinship
# cv = this_covar
# a_thresh = 0.05
# n_perm = 50
run_scan1 = function(gp, ph, ks, cv, a_thresh = 0.05, n_perm = 100) {
	# get LOD values
	lod_res = scan1(
		genoprobs = gp, 
		pheno = ph,
		kinship = ks,
		addcovar = cv
	) %>% qtl2Totbl

	# permutation testing
	perm_res = scan1perm(
		genoprobs = gp, 
		pheno = ph,
		kinship = ks,
		addcovar = cv,
		n_perm = n_perm
	) %>% qtl2Totbl

	# aggregate permutation threshold
	lod_thresh = perm_res %>% pull(pheno) %>% 
		quantile(probs = 1-a_thresh, names = FALSE)

	# convert LOD to p-value
	lod_res = lod_res %>% mutate(pval = lodToPval(pheno)) %>% rename(LOD = pheno)

	# add threshold
	lod_res$thresh = lod_thresh

	# output
	return(lod_res)
}

# run scan1 and scan1 coef and merge the data (multi-phenotype matrix version)
# gp = this_probs_drop
# ph = this_expr_vals
# ks = this_kinship
# cv = this_covar
# a_thresh = 0.05
# n_perm = 50
run_scan1_multi_pheno = function(gp, ph, ks, cv, a_thresh = 0.05, n_perm = 100) {
	# get LOD values
	lod_res = scan1(
		genoprobs = gp, 
		pheno = ph,
		kinship = ks,
		addcovar = cv
	) %>% qtl2Totbl

	# permutation testing
	perm_res = scan1perm(
		genoprobs = gp, 
		pheno = ph,
		kinship = ks,
		addcovar = cv,
		n_perm = n_perm
	) %>% qtl2Totbl

	# aggregate permutation threshold (one per phenotype)
	lod_thresh = perm_res %>% 
		summarise(across(!marker, ~quantile(.x, probs = 1-a_thresh))) %>%
		pivot_longer(cols = everything(), names_to = 'ProbeSet', values_to = 'thresh')

	# convert LOD to p-value
	lod_res = lod_res %>% 
		pivot_longer(cols = !marker, names_to = 'ProbeSet', values_to = 'LOD') %>%
		mutate(pval = lodToPval(LOD))

	# add threshold
	lod_res = lod_res %>% left_join(lod_thresh, by = 'ProbeSet')

	# output
	return(lod_res)
}
		    
# output template
tmplt_res = tibble(
	gene_expr_qtl = NA,
	perc_expand_qtl = NA,
	maf = NA,
	N = NA,
	gene_expr_sdY = NA,
	perc_expand_sdY = NA
)

# load probe information from analysis cache
{
	# load both probe info and variants per snp/probe
	probe_info = readRDS('../../../data/gene_expr/probe_to_gene.rds')
	vars_per_probe_strain = readRDS('../../../data/gene_expr/vars_per_probe_strain.rds')

	# fill strains with 0 that don't explicitly appear
	vars_per_probe_strain = vars_per_probe_strain %>% 
		complete(probe, strain = strain_info %>% filter(is_seq_str) %>% pull(bxd_id), fill = list(n_var = 0))
}

# Load expression data
{
	# read expression values
	expr_vals = readRDS('../../../data/gene_expr/expr_vals.rds')

	# some probes were filtered out so filtered these out from expression values
	expr_vals = expr_vals %>%
		semi_join(probe_info %>% distinct(probe), by = c('ProbeSet' = 'probe'))

	# make list of datasets to process
	gns = expr_vals %>% pull(GN) %>% unique

	# check completed
	if (0) {
		completed_gns = dir_ls(out_dir, glob = '*.rds') %>% 
			path_file %>% path_ext_remove %>%
			str_replace('_qtl_res', '')
		incomplete_gns = gns %>% discard(~.x %in% completed_gns)
		which(gns %in% incomplete_gns)
	}
}

# Calculate %expanded phenotype
{
	# load denovo str list
	denovo_strs_file = '../../../data/denovo_info/denovo_ri_gts_hom.tsv'
	denovo_strs = read_tsv(denovo_strs_file,
		col_types = cols(
		chr = col_character(),
		pos = col_integer(),
		end = col_integer(),
		strain = col_character(),
		RN_A = col_integer(),
		RN_B = col_integer(),
		founder = col_character(),
		fou_het = col_logical(),
		fou_gt = col_character(),
		fou_rn = col_integer(),
		delta_fou = col_integer(),
		expand_sign = col_integer(),
		expand_type = col_character()
		)
	)
	denovo_strs = denovo_strs %>% mutate(RN_T = RN_A + RN_B) %>% rename(founder_rn = fou_rn)

	# load motif info
	motif_info_file = path('../../../../090520_unified_workflow/data/ref/motif_comp', 'str_regions_mm10_filt_w_hom.tsv.gz')

	# load motif info
	motif_info = read_tsv(motif_info_file, 
		col_types = cols(
			chr = col_character(),
			pos = 'i',
			end = 'i',
			motif_len = 'i',
			motif = 'c',
			unq_motif = 'c',
			unq_motif_len = 'i',
			A = 'i',
			C = 'i',
			G = 'i',
			T = 'i',
			canon_motif = 'c',
			canon_unq_motif = 'c'
		)
	) 

	# join motif length to denovo STRs (note: these are homozygous)
	denovo_strs = denovo_strs %>%
		left_join(motif_info %>% select(chr, pos, end, motif_len), by = c('chr', 'pos', 'end'))

	# load function for calculating phenotypes
	source('../../../BXDstrs_package/BXDstrs/R/calc_pheno_vals.R')

	# calculate %expanded phenotypes
	# use only 4-mer STRs
	mls = c(4) %>% set_names(.)
	pheno_vals = map_df(mls, function(motif_len) {
		# print(motif_len)
		# calculate phenotypes
		pheno_vals = calc_pheno_vals(denovo_strs, 'proportion_expanded', 
						 min_pts_per_phe = 0, 
						 max_denovo_strains_per_loc = 10, ml = motif_len)
		pheno_vals = pheno_vals %>% mutate_at('motif_len', as.character)
	}, .id = 'motif_len')
}

# Load qtl data
{
	# config
	qtl_data_dir = '../../../data/snp_qtl2/gw'

	# load qtl2 formatted objects
	# genotype probabilities, physical map and strain kinship
	snp_probs   = readRDS(path(qtl_data_dir, 'probs.rds'))
	snp_pmap    = readRDS(path(qtl_data_dir, 'pmap.rds'))
	snp_kinship = readRDS(path(qtl_data_dir, 'kinship.rds'))

	# collate into one list
	qtl_data = list(
		snp_probs   = snp_probs, 
		snp_pmap    = snp_pmap, 
		snp_kinship = snp_kinship)
	# names(qtl_data)
}

# QTL mapping for gene expression; only chr13
{
	# NOTE: it makes sense to redo %expanded QTL mapping with the same strains which are available for gene expression

	# set strains
	strains = strain_info %>% filter(is_seq_str) %>% pull(bxd_id)

	# select dataset
	GN = gns[PBS_ARRAYID]

	# user out
	print(sprintf("Processing %s", GN))

	# select probes that have some variants and for which eqtls need to calculte with # snps as covariate
	probes = expr_vals %>% 
		filter(GN == !!GN) %>%
		semi_join(probe_info %>% filter(n_var_per_probe != 0), by = c('ProbeSet' = 'probe')) %>%
		pull(ProbeSet) %>% unique

	# subset probs and kinship
	expr_strains = expr_vals %>% filter(GN == !!GN) %>% pull(strain) %>% unique
	this_probs = qtl_data$snp_probs[,'13']
	this_kinship = qtl_data$snp_kinship['13']

	# calculate MAF
	# minprob = 0 makes it so that all probabilities are converted to gentype and no NAs returned
	maf = maxmarg(this_probs, minprob = 0)[['13']] %>% 
		t() %>% qtl2Totbl %>%
		rowwise(id = 'marker') %>% 
		# genotype counts
		summarise(
			table(c_across(everything())) %>% 
			as.list %>% as_tibble, 
			.groups = 'drop') %>%
		rename(c(B = '1', D = '2')) %>%
		mutate(across(c(B, D), ~replace_na(.x, 0))) %>%
		# find smaller of the two
		mutate(minor = pmin(B, D),
			   maf = minor/(B + D)) %>%
		select(marker, maf)

	# remove any loci with 0 maf because coloc will not accept this
	marks_to_keep = map(this_probs, function(.x) {
		marks = dimnames(.x)[[3]]
		marks[!marks %in% (maf %>% filter(maf == 0) %>% pull(marker))]
	})

	# NOTE: need to do this this way becase class of this_probs need to be preseved and not make it a generic list (otherwise use map2)
	this_probs_drop = this_probs
	for (chrom in names(this_probs_drop)) {
		this_probs_drop[[chrom]] = this_probs_drop[[chrom]][,,marks_to_keep[[chrom]]] 
	}

	# check
	# class(this_probs) == class(this_probs_drop)
	# map(this_probs, dim)
	# map(this_probs_drop, dim)

	# subset strains and chromosomes
	this_perc_expanded = pheno_vals %>% 
		filter(motif_len == 4) %>% 
		filter(strain %in% expr_strains) %>%
		select(strain, pheno) %>% 
		column_to_rownames(var = 'strain') %>% as.matrix

	# subset covariates (for %expanded we don't need number of snps as covariate, because mutator phenotype isn't connected to array probes
	this_covar = strain_info %>% 
		filter(is_seq_str) %>%
		filter(bxd_id %in% expr_strains) %>%
		select(bxd_id, gen_inbreeding) %>% 
		column_to_rownames(var = 'bxd_id') %>% as.matrix
	
	# run qtl analysis on %expanded
	# perc_expand_qtl = run_scan1(this_probs_drop, this_perc_expanded, this_kinship, this_covar, 0.05, n_perm)
	# subselect for common strains again to be certain this is correct, but checked for a single one and it works even without realigning
	use_strains = intersect(this_perc_expanded %>% row.names, expr_strains)
	perc_expand_qtl = run_scan1(
		this_probs_drop[use_strains], 
		this_perc_expanded[use_strains,,drop = FALSE], 
		this_kinship, 
		this_covar[use_strains,,drop = FALSE],
		0.05, n_perm)

	# loop over probes because we need a separate set of covariates for every probe, namely the number of snps per probe
	# library(progress)
	# pb = progress_bar$new(total = length(probes))
	# probe = probes[1]
	gene_expr_qtl = map_df(probes, function(probe) {
		# user out
		print(sprintf("Processing %s: %s", GN, probe))
		# pb$tick()

		# subset strains and chromosomes for expression values
		this_expr_vals = expr_vals %>% 
			filter(GN == !!GN, ProbeSet == !!probe) %>% 
			filter(strain %in% strains)

		# filter out zero-variance expression values
		this_expr_vals = this_expr_vals %>%
			group_by(ProbeSet) %>%
			mutate(stdev = sd(expr_val, na.rm = TRUE)) %>%
			filter(stdev != 0) %>%
			ungroup

		# NOTE: that for expression datasets there may be a lot fewer strains for which phenotype are known then genotypes
		this_expr_vals = this_expr_vals %>%
			select(strain, pheno = expr_val, pheno_id = ProbeSet) %>%
			pivot_wider(id_cols = strain, names_from = pheno_id, values_from = pheno) %>%
			column_to_rownames(var = 'strain') %>% as.matrix

		# subset covariates
		this_covar = strain_info %>% 
			filter(is_seq_str) %>%
			filter(bxd_id %in% expr_strains) %>%
			select(bxd_id, gen_inbreeding) %>% 
			# add number of variants per probe per strain
			left_join(vars_per_probe_strain %>% 
						filter(probe == !!probe) %>% 
						select(strain, n_var), by = c('bxd_id' = 'strain')) %>%
			column_to_rownames(var = 'bxd_id') %>% as.matrix

		# run eqtl analysis on gene expression
		# gene_expr_qtl = run_scan1_multi_pheno(this_probs_drop, this_expr_vals, this_kinship, this_covar, 0.05, n_perm)
		use_strains = intersect(this_expr_vals %>% row.names, expr_strains)
		gene_expr_qtl = run_scan1_multi_pheno(
			this_probs_drop[use_strains], 
			this_expr_vals[use_strains,,drop = FALSE], 
			this_kinship, 
			this_covar[use_strains,,drop = FALSE],
			0.05, n_perm)

		# # form a tibble suitable for output
		# # keep this tibbles separate as list elements, because we'll need to pass them separately to `coloc` package anyway
		# res = tibble(
		# 	gene_expr_qtl = list(gene_expr_qtl),
		# 	perc_expand_qtl = list(perc_expand_qtl),
		# 	maf = list(maf),
		# 	N = length(expr_strains), # number of samples 
		# 	GN = GN
		# ) %>% relocate(GN, .before = 1)
		# 
		# # output
		# return(res)
		return(gene_expr_qtl)
	})

	# combine
	res = tibble(
		gene_expr_qtl = list(gene_expr_qtl), 
		perc_expand_qtl = list(perc_expand_qtl),
		maf = list(maf),
		N = length(expr_strains), # number of samples 
		GN = GN,
	)

	# save result individually
	saveRDS(res, path(out_dir, str_c(GN, '_qtl_res'), ext = 'rds'))
}

