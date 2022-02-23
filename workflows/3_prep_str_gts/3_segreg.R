#!/home/momaksimov/anaconda3/envs/r4.2/bin/Rscript

# about: find segregating str loci within BXD

# options
options(stringsAsFactors = FALSE)

# libraries
library(tidyverse)
library(fs)
library(DBI)
library(dbplyr)
library(progress)

# for parallel
# library(furrr)
# library(progressr)

# config
out_dir = '../../data/str_gts'; dir_create(out_dir)
geno_file = 'all_repcn_proc_nosegdup_nolowcr.rds'

# read processed version of all the genotypes
geno = readRDS(path(out_dir, geno_file))

# blank out hets (this is fast column-wise)
# NAs are already blanked out
geno_het_clean = geno %>%
    # select(chr, pos, end, BXD001, BXD002, C57BL, DBA) %>%
    mutate(across(c(matches('BXD'), C57BL, DBA), function(.x) {
	.y = .x %>% str_split(',', simplify = TRUE)
	.x[.y[,1] != .y[,2]] = NA
	.x
    }))

# count number of unique genotypes per row
# CHECK: BXD_variants_MM_strs_highQ_calls_strain_filt.vcf.gz with 154 strains (BXD210 filtered)
# CHECK: founders are counted in number of unique genotypes

cache_file = path(out_dir, 'unq_gt_per_loc.tsv')
if (!file_exists(cache_file)) {
    pb <- progress_bar$new(total = dim(geno_het_clean)[1])
    n_unq_gt = geno_het_clean %>% 
	select(matches('BXD'), C57BL, DBA) %>%
	apply(MARGIN = 1, FUN = function(.x) {pb$tick(); .x %>% discard(is.na) %>% unique %>% length})
    unq_gt_per_loc = geno %>% select(chr, pos, end) %>% mutate(n_unq_gt = !!n_unq_gt)

    # more complicated way with parallel furrr
    # plan(multisession, workers = 10)
    # n_batch = 10
    # with_progress({
    #     p <- progressor(steps = n_batch)
    #     gt_per_loc = geno_het_clean %>%
    # 	# slice(1:50000) %>%
    # 	mutate(batch = cut_interval(row_number(), n = n_batch, label = FALSE)) %>%
    # 	nest(data = !batch) %>%
    # 	mutate(res = future_map(data, function(.x) {
    # 	    p()
    # 	    .x$n_gt_per_loc = .x %>% 
    # 		select(matches('BXD'), C57BL, DBA) %>% 
    # 		apply(MARGIN = 1, FUN = function(.x) .x %>% discard(is.na) %>% unique %>% length)
    # 	    .x %>% select(chr, pos, end, n_gt_per_loc)
    # 	}))
    # })
    # gt_per_loc = gt_per_loc %>% select(batch, res) %>% unnest(res) %>% select(!batch)

    # save unq_gt_per_loc
    write_tsv(unq_gt_per_loc, cache_file)
} else { unq_gt_per_loc = read_tsv(cache_file, col_types = cols(chr = 'c', pos = 'i', end = 'i', n_unq_gt = 'i')) }

# spot check some of these
geno_het_clean %>%
    semi_join(unq_gt_per_loc %>% filter(n_unq_gt == 2) %>% slice_sample(n = 1)) %>%
    pivot_longer(matches('BXD'), 'var', 'val') %>%
    count(value)

# make list of segregating strs
segreg = geno %>%
    semi_join(unq_gt_per_loc %>% filter(n_unq_gt >= 2), by = c("chr", "pos", "end"))

# make a list of monoallelic loci
monoal_loci = unq_gt_per_loc %>% filter(n_unq_gt < 2) %>% select(chr, pos, end)

{ # DEBUGGING
    # check with previous lists
    ori_regions = read_tsv('/projects/ps-gymreklab/mikhail/090520_unified_workflow/data/ref/filt/str_regions_wo_zv_wo_hetonly_wo_segdup_wo_lowcr.bed',
	col_names = c('chr', 'pos', 'end'),
	col_types = cols(chr = 'c', pos = 'i', end = 'i')
    )
    list(ori_regions = ori_regions, segreg = segreg) %>% map(nrow)
    setdiff(ori_regions, segreg %>% select(chr, pos, end)) # note the extra 93 new excluded lowcr loci
    setdiff(segreg %>% select(chr, pos, end), ori_regions)
}
   
# sort list by locus
segreg = segreg %>% 
    mutate(chr = fct_relevel(chr, str_c('chr', 1:19))) %>%
    arrange(chr, pos, end)

# write segregating genotypes and locus list
saveRDS(segreg, path(out_dir, 'all_repcn_proc_nosegdup_nolowcr_segreg.rds'))
write_tsv(segreg %>% select(chr, pos, end), '../../data/ref/str_segreg.tsv')
write_tsv(monoal_loci %>% select(chr, pos, end), '../../data/ref/str_zv_hetonly.tsv')
