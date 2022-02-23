#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript
#PBS -q hotel
#PBS -N qtl_map_gene_expr
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
#PBS -o job.o
#PBS -e job.e 
#PBS -V
#PBS -M mikmaksi@gmail.com
#PBS -m abe
#PBS -t 1-54
#PBS -A gymreklab-group

# about: run association testing between snp founder genotypes and gene expression
# run qtl2 in multi-phenotype mode (probe expression = phenotype)
# NOTE: this requires first running ../../analysis/4_eqtl_mapping.Rmd

# clean variables
# rm(list = ls())

# options
options(stringsAsFactors = FALSE)

# cd into analysis directory
setwd('/projects/ps-gymreklab/mikhail/072321_bxd_mutator_paper/workflows/6b_gene_expr_mapping_new_wind')

# get job id
PBS_ARRAYID = as.integer(Sys.getenv('PBS_ARRAYID'))
# PBS_ARRAYID = 2
print(sprintf('Processing job_id %s', PBS_ARRAYID))

# libraries
library(tidyverse)
library(cowplot)
library(fs)
library(qtl2)
library(coloc)
library(jsonlite)

# load strain info
strain_info = readRDS('../../info/strain_info.rds')

# config
# min_strains = 5 # number of strains required to run association testing
# NOTE: turning off permutation testing, b/c we care more about the resolution of the peak and less about the significance, since we already demonstrated significance using the sparse dataset
n_perm = 100 
do_perm = FALSE
get_coef = TRUE
    
# set up output directory
out_dir = '../../data/gene_expr_new_wind/chr13qtl_expr_qtlres_dense'; dir_create(out_dir)

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
# ph = this_perc_expanded
# ks = this_kinship
# cv = this_covar
# a_thresh = 0.05
# n_perm = 50
# get_coef = TRUE
run_scan1 = function(gp, ph, ks, cv, a_thresh = 0.05, n_perm = 100, do_perm = TRUE, get_coef = FALSE) {
    # get LOD values
    lod_res = scan1(
	genoprobs = gp, 
	pheno = ph,
	kinship = ks,
	addcovar = cv
    ) %>% qtl2Totbl
    
    # convert LOD to p-value
    lod_res = lod_res %>% mutate(pval = lodToPval(pheno)) %>% rename(LOD = pheno)
    
    # permutation testing
    if (do_perm) {
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
    } else {
	lod_thresh = NA
    }

    # add threshold
    lod_res$thresh = lod_thresh

    # calculate effect sizes
    if (get_coef) {
	coef_res = scan1coef(
	    genoprobs = gp, 
	    pheno = ph,
	    kinship = ks,
	    addcovar = cv,
	    se = TRUE
	)

	# extract the standard errors
	se_res = attr(coef_res, 'SE') %>% qtl2Totbl
	
	# for coefficients
	coef_res = coef_res %>% qtl2Totbl
	
	# join errors to effect sizes
	# calculate the variance
	coef_res = coef_res %>%
	    select(marker, Beta = BB) %>%
	    left_join(se_res %>% select(marker, se = BB), by = 'marker') %>%
	    mutate(var = se^2) %>%
	    select(!se)

	# join effect sizes to p-values
	lod_res = lod_res %>%
	    left_join(coef_res, by = 'marker')
    }

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
# do_perm = FALSE
# get_coef = TRUE
run_scan1_multi_pheno = function(gp, ph, ks, cv, a_thresh = 0.05, n_perm = 100, do_perm = TRUE, get_coef = FALSE) {
    # get LOD values
    lod_res = scan1(
	genoprobs = gp, 
	pheno = ph,
	kinship = ks,
	addcovar = cv
    ) %>% qtl2Totbl
    
    # convert LOD to p-value
    lod_res = lod_res %>% 
	pivot_longer(cols = !marker, names_to = 'ProbeSet', values_to = 'LOD') %>%
	mutate(pval = lodToPval(LOD))

    # permutation testing
    if (do_perm) {
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
    } else {
	lod_thresh = tibble(ProbeSet = lod_res %>% pull(ProbeSet) %>% unique, lod_thresh = NA)
    }
	
    # add threshold
    lod_res = lod_res %>% left_join(lod_thresh, by = 'ProbeSet')

    # calculate effect sizes
    if (get_coef) {
	coef_res = map_df(colnames(ph) %>% set_names(.), function(.x) {
	    # user out
	    # print(.x)

	    coef_res = scan1coef(
		genoprobs = gp, 
		pheno = ph[,.x],
		kinship = ks,
		addcovar = cv,
		se = TRUE
	    )

	    # extract the standard errors
	    se_res = attr(coef_res, 'SE') %>% qtl2Totbl
	    
	    # for coefficients
	    coef_res = coef_res %>% qtl2Totbl
	    
	    # join errors to effect sizes
	    # calculate the variance
	    coef_res = coef_res %>%
		select(marker, Beta = BB) %>%
		left_join(se_res %>% select(marker, se = BB), by = 'marker') %>%
		mutate(var = se^2) %>%
		select(!se)
	    return(coef_res)
	}, .id = 'ProbeSet')

	# join effect sizes to p-values
	lod_res = lod_res %>%
	    left_join(coef_res, by = c('ProbeSet', 'marker'))
    }
    # output
    return(lod_res)
}

# quiet version of coloc.abf
coloc.abf.quiet = quietly(coloc.abf)
		    
# output template
tmplt_res = tibble(
    gene_expr_qtl = NA,
    perc_expand_qtl = NA,
    maf = NA,
    N = NA,
    gene_expr_sdY = NA,
    perc_expand_sdY = NA
)

# Read list of representative GN datasets
{
    # read sets
    repr_dsets = readRDS('../../data/analysis_cache/final_gn_table.rds')
}

# Load expression data
{
    # read expression values
    expr_vals = readRDS('../../data/gene_expr_new_wind/chr13qtl_expr_vals.rds')
    
    # filter in representative datasets
    expr_vals = expr_vals %>% filter(GN %in% (repr_dsets %>% pull(sel_dset)))
    
    # make list of datasets to process
    gns = expr_vals %>% pull(GN) %>% unique
    # gns = gns[c(1, 39)]
    
    # subset to probes (especially important because coloc is also being run)
    if (0) { # instead of just representative datasets also subset to best probe per gene based on sparse qtl mapping data
	best_point = readRDS('../../data/gene_expr/chr13qtl_agg/best_point.rds')
	expr_vals = expr_vals %>% semi_join(best_point, by = c('GN', 'ProbeSet' = 'probe'))
    } else {
	# can reduce computational load even further if only GN/gene pairs with eQTL signals are evaluated (from sparse snp analysis)
	# NOTE: filtering this way means that script may fail for certain datasets if no genes have significant eQTL in that dataset
	eqtl_data = readRDS('../../data/analysis_cache/eqtl_data.rds')
	probe_genes_to_keep = eqtl_data$eqtl_dsets_genes %>%
	    left_join(eqtl_data$unq_expr_vals$top_qtl_probe %>% 
		      distinct(GN, gene_name, probe), by = c('GN', 'gene_name'))
	expr_vals = expr_vals %>% semi_join(probe_genes_to_keep, by = c('GN', 'ProbeSet' = 'probe'))
    }
    
    # check completed
    if (0) {
	completed_gns = dir_ls(out_dir, glob = '*.rds') %>% 
	    path_file %>% path_ext_remove %>%
	    str_replace('_qtl_res', '')
	incomplete_gns = gns %>% discard(~.x %in% completed_gns)
	which(gns %in% incomplete_gns)

	# some of the incomplete gns may not have had any genes with eQTL in filter set
	# if this frame has no rows then analysis is done
	probe_genes_to_keep %>% filter(GN %in% incomplete_gns)
    }
}

# Calculate %expanded phenotype
{
    # load denovo str list
    denovo_strs_file = '../../data/denovo_info/denovo_ri_gts_hom.tsv'
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
    motif_info_file = path('../../../090520_unified_workflow/data/ref/motif_comp', 'str_regions_mm10_filt_w_hom.tsv.gz')

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
    source('../../BXDstrs_package/BXDstrs/R/calc_pheno_vals.R')

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
    qtl_data_dir = '../../data/snp_qtl2/chr13'

    # load qtl2 formatted objects
    # genotype probabilities, physical map and strain kinship
    snp_probs   = readRDS(path(qtl_data_dir, 'probs.rds'))
    snp_pmap    = readRDS(path(qtl_data_dir, 'pmap.rds'))
    snp_kinship = readRDS(path('../../data/snp_qtl2/gw', 'kinship.rds')) # use gw version of kinship

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

    # subset strains and chromosomes for expression values
    this_expr_vals = expr_vals %>% 
	filter(GN == !!GN) %>% 
	filter(strain %in% strains)

    # filter out zero-variance expression values
    this_expr_vals = this_expr_vals %>%
	group_by(ProbeSet) %>%
	mutate(stdev = sd(expr_val, na.rm = TRUE)) %>%
	filter(stdev != 0) %>%
	ungroup
    
    # NOTE: that for expression datasets there may be a lot fewer strains for which phenotype are known then genotypes
    expr_strains = this_expr_vals %>% pull(strain) %>% unique
    this_expr_vals = this_expr_vals %>%
	select(strain, pheno = expr_val, pheno_id = ProbeSet) %>%
	pivot_wider(id_cols = strain, names_from = pheno_id, values_from = pheno) %>%
	column_to_rownames(var = 'strain') %>% as.matrix

    # check if there are enough strains
    # NOTE: remove this for next reprocessing, b/c there aren't actually any datasets like this
    # if (length(expr_strains) < min_strains) {
    #     print(sprintf('Only %d strains ... skipping', length(expr_strains)))
    #     quit(save = 'no')
    # }

    # subset probs, kinship and covariates
    this_probs = qtl_data$snp_probs[expr_strains,'13']
    this_kinship = qtl_data$snp_kinship['13']
    this_covar = strain_info %>% 
	select(bxd_id, gen_inbreeding) %>% 
	filter(bxd_id %in% expr_strains) %>%
	column_to_rownames(var = 'bxd_id') %>% as.matrix

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

    # run qtl analysis on gene expression
    gene_expr_qtl = run_scan1_multi_pheno(this_probs_drop, this_expr_vals, this_kinship, this_covar, 0.05, n_perm, do_perm, get_coef)

    # subset strains and chromosomes
    this_perc_expanded = pheno_vals %>% 
	filter(motif_len == 4) %>% 
	filter(strain %in% expr_strains) %>%
	select(strain, pheno) %>% 
	column_to_rownames(var = 'strain') %>% as.matrix
    
    # run qtl analysis on %expanded
    perc_expand_qtl = run_scan1(this_probs_drop, this_perc_expanded, this_kinship, this_covar, 0.05, n_perm, do_perm, get_coef)

    # .x = '1418913_at'; .y = gene_expr_qtl %>% filter(ProbeSet == .x)
    res = gene_expr_qtl %>%
	nest(gene_expr_qtl = !ProbeSet) %>%
	mutate(coloc_res = map2(ProbeSet, gene_expr_qtl, function(.x, .y) {
	    # user out
	    print(.x) 

	    # calculate standard deviations
	    gene_expr_sdY   = sd(this_expr_vals[,.x])
	    perc_expand_sdY = sd(this_perc_expanded[,'pheno'])

	    # assemble list for coloc
	    signals = list(perc_expand_qtl = perc_expand_qtl, gene_expr_qtl = .y) %>%
		map(function(.x) {
			.x = .x %>% left_join(maf, by = 'marker')
			list(beta = .x$Beta, 
			     varbeta = .x$var,
			     N = length(expr_strains),
			     MAF = .x$maf,
			     snp = .x$marker,
			     type = 'quant')
		})
	    # check
	    # names(signals)

	    # add sdY to both
	    signals$gene_expr_qtl$sdY   = gene_expr_sdY
	    signals$perc_expand_qtl$sdY = perc_expand_sdY

	    # run coloc
	    coloc_res = coloc.abf.quiet(signals$gene_expr_qtl, signals$perc_expand_qtl)
	    # coloc_res = coloc.abf(signals$gene_expr_qtl, signals$perc_expand_qtl)
	    return(coloc_res)
	}))

    # assembly output frame (renest)
    res = tibble(
	gene_expr_qtl = list(res %>% select(ProbeSet, gene_expr_qtl) %>% unnest(gene_expr_qtl)),
	perc_expand_qtl = list(perc_expand_qtl),
	coloc_res = list(res %>% select(ProbeSet, coloc_res) %>% unnest(coloc_res)),
	GN = GN,
	maf = list(maf),
	N = length(expr_strains) # number of samples 
    ) %>% relocate(GN, .before = 1)

    # save result individually
    saveRDS(res, path(out_dir, str_c(GN, '_qtl_res'), ext = 'rds'))
}
