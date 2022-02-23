#!/home/momaksimov/anaconda3/envs/r4/bin/Rscript

# about: prep qtl2 compatible data (probabilties and maps) for individual chromosomes using plink-pruned snps
# all > plink-pruned > ld-indep
# these are dense marker maps that we didn't end up using

# clean variables
rm(list = ls())

# options
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

# libraries
library(tidyverse)
library(fs)
library(qtl2)
library(jsonlite)

# for parallel
library(furrr)
library(progressr)

# config
snp_geno_dir = '../../data/snp_gts/fou_lab'
out_dir = '../../data/snp_qtl2'; dir_create(out_dir)
chroms = str_c('chr', c(4, 13, 17))
parallel = FALSE
prune_dir = '../../data/plink'
	
# for coordinate interpolation
gw_dir = '../../data/snp_qtl2/gw'
gw_pmap = readRDS(path(gw_dir, 'pmap.rds')) %>% map_df(~tibble(marker = names(.x), pos = .x), .id = 'chr')
gw_gmap = readRDS(path(gw_dir, 'gmap.rds')) %>% map_df(~tibble(marker = names(.x), cM = .x), .id = 'chr')
interp_fns = gw_gmap %>%
    left_join(gw_pmap, by = c('chr', 'marker')) %>%
    nest(data = !chr) %>%
    group_by(chr) %>%
    summarise(fn = map(data, ~approxfun(x = .x$pos, y = .x$cM, rule = 2)))
    
if (parallel) { plan(multicore, workers = 10) } else { plan(sequential) }
with_progress({
    p <- progressor(steps = length(chroms))
    # chrom = 'chr17'
    future_walk(chroms, function(chrom) {
	# read genotypes
	geno = read_tsv(path(snp_geno_dir, str_c(chrom, '.tsv.gz')), col_types = cols(.default = 'c'))

	# subset to list of pruned loci
	pruned_in = read_tsv(path(prune_dir, str_c(chrom, '.prune.in')),
			     col_names = c('chr', 'pos', 'end'), col_types = cols(chr = 'c', pos = 'i', end = 'i')) %>%
	    unite('locus', c('chr', 'pos', 'end'))
	geno = geno %>% semi_join(pruned_in, by = 'locus')

	# separate out 
	geno = geno %>% separate(locus, c('chr', 'pos', 'end'), convert = TRUE, remove = FALSE)

	# keep only those markers of interest in the region
	# geno = geno %>% filter(pos >= 86e6, end <= 96e6)

	# add interpolated cM
	interp_fn = (interp_fns %>% filter(chr == str_replace(chrom, 'chr', '')) %>% pull(fn))[[1]]
	geno = geno %>% mutate(cM = interp_fn(geno$pos/1e6))

	# make physical and genetic maps
	geno_map = geno %>% select(chr, marker = locus, cM)
	phys_map = geno %>% select(chr, marker = locus, pos)

	# write outputs in qtl2 format
	sub_dir = path(out_dir, chrom); dir_create(sub_dir)
	phys_map %>% 
	    select(marker, chr, pos) %>%
	    mutate(chr = str_replace(chr, 'chr', ''),
		   pos = pos/1e6) %>% 
	    write_csv(path(sub_dir, 'pmap.csv'))
	geno_map %>% 
	    select(marker, chr, pos = cM) %>%
	    mutate(chr = str_replace(chr, 'chr', '')) %>% 
	    write_csv(path(sub_dir, 'gmap.csv'))
	geno %>%
	    rename(marker = locus) %>%
	    select(!c(chr, pos, end, cM)) %>%
	    relocate(marker, .before = 1) %>%
	    write_csv(path(sub_dir, 'geno.csv'))

	# make a config file
	config_file = path(sub_dir, 'config.json')
	write_control_file(config_file, 
	    description = str_c('BXD snps on ', chrom),
	    crosstype = 'risib',
	    sep = ",",
	    na.strings = c("-", "NA"),
	    comment.char = "#",
	    geno_file = sprintf('../../../data/snp_qtl2/%s/geno.csv', chrom),
	    gmap_file = sprintf('../../../data/snp_qtl2/%s/gmap.csv', chrom),
	    pmap_file = sprintf('../../../data/snp_qtl2/%s/pmap.csv', chrom),
	    alleles = c("B", "D"),
	    geno_codes = c(B = 1, D = 2),
	    geno_transposed = TRUE, 
	    overwrite = TRUE)
	# load the cross object
	bxd_cross = read_cross2(config_file, quiet = FALSE)
	# file_delete(config_file)

	# calculate genotype probabilities
	snp_probs = calc_genoprob(bxd_cross, cores = 10, quiet = FALSE)

	# save rds objects
	saveRDS(snp_probs,      file = path(sub_dir, 'probs.rds'))
	saveRDS(bxd_cross$gmap, file = path(sub_dir, 'gmap.rds'))
	saveRDS(bxd_cross$pmap, file = path(sub_dir, 'pmap.rds'))
    })
})
